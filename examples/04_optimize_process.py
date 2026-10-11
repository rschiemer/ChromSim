# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Tutorial 04 — Optimise the elution gradient under a purity constraint
#
# With the calibrated model from Tutorial 03 we now design the process. The mAb must be pooled at
# **≥ 99 % purity**, with HMW as the impurity (by mass). We will:
#
# 1. evaluate the current method: the best pool, its yield (recovery) and productivity;
# 2. find the gradient length and end salt that give the highest recovery at 99 % purity;
# 3. confirm the optimum on a finer grid;
# 4. map the trade-off between recovery and productivity (Pareto front).
#
# The cut points are not guessed. For every candidate method, ChromSim chooses the pool exactly:
# the start and end that give the highest recovery while meeting the purity.

# %%
from pathlib import Path

import numpy as np

from chromsim import api, units

ROOT = Path.cwd() if (Path.cwd() / "examples").is_dir() else Path.cwd().parent
project = api.load_project(ROOT / "examples" / "projects" / "cex_mab_hmw_lge.yaml")
project.isotherm.parameters["mAb"]["k_eq"] = 1.5e-2  # calibrated in Tutorial 03
project.solver.n_col, project.solver.n_par = 20, 4  # coarse grid for the search; checked below

PURITY = 0.99
TARGET = "mAb"
IMPURITIES = api.default_impurities(project, TARGET)  # every component except salt and mAb
print("impurities:", IMPURITIES)


def report(label, performance):
    """One line of KPIs in lab units."""
    productivity = units.from_si(performance.productivity_kg_per_m3_s, "g/L/h")
    buffer = units.from_si(performance.buffer_consumption_m3_per_kg, "L/g")
    print(
        f"{label:<22} recovery {performance.recovery:6.1%}  purity {performance.purity:6.2%}  "
        f"productivity {productivity:5.2f} g/(L·h)  buffer {buffer:6.1f} L/g  "
        f"pool {units.from_si(performance.volume_m3, 'mL'):5.1f} mL"
    )


# %% [markdown]
# ## 1. The current method
#
# The pool may only be collected during the elution step (`pool_steps`). Otherwise, the highest
# recovery could come from collecting traces of product during the wash, which dilutes the pool.

# %%
baseline_result = api.simulate(project)
baseline_pool = api.optimal_pool(baseline_result, TARGET, IMPURITIES, PURITY, pool_steps=("elute",))
baseline = api.pool_performance(baseline_result, baseline_pool, IMPURITIES)
report("current method", baseline)
api.plot_chromatogram(baseline_result, x="CV", pools=[baseline_pool])

# %% [markdown]
# ## 2. Maximise recovery at 99 % purity
#
# The decision variables are entries of the project, addressed by path, as in fitting:
# the gradient length (step 3, in CV) and the salt concentration at its end (in mM).

# %%
variables = [
    api.FitParameter(path="method.3.duration", lower=5, upper=40, unit="CV", name="gradient_CV"),
    api.FitParameter(
        path="method.3.inlet_end.salt", lower=250, upper=1000, unit="mM", name="end_salt_mM"
    ),
]


def progress(evaluation, objectives):
    if evaluation % 10 == 0:
        print(f"  {evaluation} methods simulated")


best = api.optimize_process(
    project,
    variables,
    TARGET,
    PURITY,
    objectives=("recovery",),
    pool_steps=("elute",),
    max_evaluations=30,
    callback=progress,
)
gradient_cv, end_salt = best.optimization.x
print(f"best method: {gradient_cv:.1f} CV gradient to {end_salt:.0f} mM")
report("optimised method", best.best.performance)
api.plot_chromatogram(best.best.result, x="CV", pools=[best.best.pool])

# %% [markdown]
# A shallower gradient that ends at lower salt separates mAb from HMW better, so far more of the
# mAb can be pooled at the same purity. The best gradient length sits at the upper bound we
# allowed (40 CV): recovery keeps rising with longer gradients. What stops us is the cost in
# cycle time, which is the trade-off of section 4.
#
# ## 3. Check the optimum on a finer grid
#
# The search used a coarse grid (20 × 4 cells) to stay fast. Before trusting the result,
# simulate the chosen method on a finer grid and pool it again.

# %%
fine = best.best.project.model_copy(deep=True)
fine.solver.n_col, fine.solver.n_par = 50, 12
fine_result = api.simulate(fine)
fine_pool = api.optimal_pool(fine_result, TARGET, IMPURITIES, PURITY, pool_steps=("elute",))
fine_performance = api.pool_performance(fine_result, fine_pool, IMPURITIES)
report("optimised, fine grid", fine_performance)
assert fine_performance.purity >= PURITY - 1e-9
assert fine_performance.recovery > baseline.recovery + 0.05, (
    "the optimised method should be clearly better"
)
assert abs(fine_performance.recovery - best.best.performance.recovery) < 0.05, (
    "grid too coarse for the search"
)

# %% [markdown]
# ## 4. Recovery or productivity?
#
# Shallow gradients raise recovery but lengthen the cycle, which lowers productivity. A
# multi-objective optimizer (NSGA-II) maps the best compromises, the Pareto front: no point on
# it can be improved in one objective without losing in the other.
#
# This short run (40 simulations) is enough to see the shape. For a real study, use more
# generations and `parallel=N` to run simulations side by side. In scripts, put the call under
# `if __name__ == "__main__":` when using `parallel`.

# %%
pareto = api.optimize_process(
    project,
    variables,
    TARGET,
    PURITY,
    objectives=("recovery", "productivity"),
    pool_steps=("elute",),
    pop_size=10,
    n_gen=4,
)
molar_mass = project.components[project.component_names.index(TARGET)].molar_mass_kg_per_mol
front = pareto.table.copy()
front["productivity_g_per_L_h"] = units.from_si(
    front["productivity"].to_numpy() * molar_mass, "g/L/h"
)
columns = ["gradient_CV", "end_salt_mM", "recovery", "productivity_g_per_L_h"]
print(front.sort_values("recovery")[columns].round(3).to_string(index=False))
api.plot_pareto(
    front,
    x="recovery",
    y="productivity_g_per_L_h",
    color="gradient_CV",
    labels={
        "recovery": "Recovery",
        "productivity_g_per_L_h": "Productivity / g/(L·h)",
        "gradient_CV": "Gradient / CV",
    },
)

# %% [markdown]
# Along the front, longer gradients buy recovery at the cost of productivity. Which point to run
# is a process decision: product value, column time and buffer cost decide. Every point meets
# the 99 % purity, because the pool is cut exactly for each candidate.

# %%
productivity = front["productivity_g_per_L_h"].to_numpy()
recovery = front["recovery"].to_numpy()
order = np.argsort(recovery)
assert np.all(np.diff(productivity[order]) <= 1e-9), (
    "a Pareto front trades one objective for the other"
)
