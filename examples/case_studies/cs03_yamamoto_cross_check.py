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
# # Case study 03 — Yamamoto method: ChromSim against CADET-Process
#
# The CADET-Workshop (lesson 99.1, CC BY-NC-SA 4.0) ships simulated linear-gradient elutions
# of three proteins on a 10 cm × 7.7 mm column at four gradient volumes (18.8–75.2 mL). Each
# file has three columns per protein set, without header: time (s), outlet salt (mM) and the
# protein concentrations (mM). Fetch them first:
#
#     pixi run python examples/case_studies/fetch_data.py cadet_workshop_yamamoto
#
# (The files are model output, but they are not perfectly consistent, as we will see.)
#
# CADET-Process has its own implementation of the Yamamoto method
# (`CADETProcess.tools.yamamoto`). Both should give the same ν and k_eq from the same data:
# this is an independent check of `api.yamamoto_fit`, `api.elution_salt_at_peak` and of the
# importer for headerless files.

# %%
from pathlib import Path

import numpy as np

from chromsim import api

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
DATA = HERE / "data" / "cadet_workshop_yamamoto"
if not (DATA / "18.8mL.csv").is_file():
    raise SystemExit(f"Data not found in {DATA}; run examples/case_studies/fetch_data.py first.")

ML = 1e-6
GRADIENT_VOLUMES_ML = (18.8, 37.6, 56.4, 75.2)
PROTEINS = ("A", "B", "C")
# Column and resin as in the workshop lesson.
COLUMN = api.Column.model_validate({"model": "LRMP", "length": "10 cm", "diameter": "7.7 mm"})
COL_POROSITY, PAR_POROSITY = 0.36, 0.85
TOTAL_POROSITY = COL_POROSITY + (1.0 - COL_POROSITY) * PAR_POROSITY
IONIC_CAPACITY = 4.7 * 175.0  # mM

MAPPING = {
    "header": False,
    "x_column": 0,
    "x_unit": "s",
    "signals": [
        {"column": 1, "name": "salt", "unit": "mM"},
        *({"column": 2 + i, "name": name, "unit": "mM"} for i, name in enumerate(PROTEINS)),
    ],
}

# %% [markdown]
# ## ChromSim
#
# Salt and proteins are recorded by the same (simulated) detector, so no delay correction is
# needed, and both signals can stay on the time axis. The gradient slope is the salt increase
# over the gradient volume.

# %%
runs = {name: [] for name in PROTEINS}
experiments = {}
for volume_ml in GRADIENT_VOLUMES_ML:
    experiment = api.read_chromatogram(DATA / f"{volume_ml}mL.csv", "generic_csv", mapping=MAPPING)
    experiments[volume_ml] = experiment
    salt = experiment.signals["salt"]
    assert salt.kind == "salt"
    slope = (salt.values.max() - salt.values.min()) / (volume_ml * ML)  # mM per m³
    for name in PROTEINS:
        elution_salt = api.elution_salt_at_peak(experiment.signals[name], salt)
        runs[name].append(api.LgeRun(f"{volume_ml} mL", slope, elution_salt))

chromsim_results = {
    name: api.yamamoto_fit(runs[name], COLUMN.volume_m3, TOTAL_POROSITY, IONIC_CAPACITY)
    for name in PROTEINS
}
for name, result in chromsim_results.items():
    print(
        f"ChromSim       {name}: ν = {result.nu:.4f}, k_eq = {result.k_eq:.4e}, "
        f"R² = {result.r_squared:.5f}"
    )

# %% [markdown]
# ## CADET-Process
#
# Same data, CADET-Process's `GradientExperiment` and `fit_parameters` (it fits all proteins at
# once and needs a column object for volume and porosity).

# %%
from CADETProcess.processModel import (  # noqa: E402  (only needed for this comparison)
    ComponentSystem,
    LumpedRateModelWithPores,
    StericMassAction,
)
from CADETProcess.tools.yamamoto import GradientExperiment, fit_parameters  # noqa: E402

component_system = ComponentSystem(["Salt", *PROTEINS])
binding_model = StericMassAction(component_system)
binding_model.capacity = IONIC_CAPACITY
column = LumpedRateModelWithPores(component_system, "column")
column.binding_model = binding_model
column.length, column.diameter = 0.1, 0.0077
column.bed_porosity, column.particle_porosity = COL_POROSITY, PAR_POROSITY
column.particle_radius = 34e-6

cadet_experiments = []
for volume_ml, experiment in experiments.items():
    salt = experiment.signals["salt"]
    proteins = np.column_stack([experiment.signals[name].values for name in PROTEINS])
    cadet_experiments.append(GradientExperiment(salt.x, salt.values, proteins, volume_ml * ML))
cadet_result = fit_parameters(cadet_experiments, column)
cadet_nu = np.atleast_1d(cadet_result.characteristic_charge)
cadet_k_eq = np.atleast_1d(cadet_result.k_eq)
for i, name in enumerate(PROTEINS):
    print(f"CADET-Process  {name}: ν = {cadet_nu[i]:.4f}, k_eq = {cadet_k_eq[i]:.4e}")

# %% [markdown]
# ## Comparison
#
# ChromSim fits the straight line in log–log space by linear least squares; CADET-Process uses
# a nonlinear curve fit of the same equation. On these noise-free data the results agree
# closely.

# %%
for i, name in enumerate(PROTEINS):
    ours = chromsim_results[name]
    nu_difference = abs(ours.nu - cadet_nu[i]) / cadet_nu[i]
    k_eq_difference = abs(ours.k_eq - cadet_k_eq[i]) / cadet_k_eq[i]
    print(f"{name}: Δν = {nu_difference:.2%}, Δk_eq = {k_eq_difference:.2%}")
    assert nu_difference < 0.01 and k_eq_difference < 0.05, f"protein {name} disagrees"


# %% [markdown]
# ## Why is R² so low for "simulated" data?
#
# Both implementations agree, yet R² is only 0.79–0.96. The residual table shows why: the
# 56.4 mL run sits off the line for **all three proteins**, and in the wrong direction. A
# shallower gradient must elute a protein at lower salt, but protein A elutes at 274 mM in the
# 56.4 mL run and at 265 mM in the steeper 37.6 mL run. The salt programme of that file is
# consistent with the others (each extra 18.8 mL of gradient adds the same ≈2200 s), so the
# protein peaks of that one run deviate. Lesson for real data: always look at the residuals,
# and check that I_R decreases monotonically with the gradient volume before trusting ν.

# %%
table = chromsim_results["A"].table[["name", "gh", "elution_salt", "residual"]]
print(table.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
elution_salt = table["elution_salt"].to_numpy()
assert not np.all(np.diff(elution_salt) < 0), "data set changed: I_R is now monotonic"

without_outlier = {
    name: api.yamamoto_fit(
        [run for run in runs[name] if run.name != "56.4 mL"],
        COLUMN.volume_m3,
        TOTAL_POROSITY,
        IONIC_CAPACITY,
    )
    for name in PROTEINS
}
for name, result in without_outlier.items():
    print(
        f"without 56.4 mL  {name}: ν = {result.nu:.3f}, k_eq = {result.k_eq:.3e}, "
        f"R² = {result.r_squared:.4f}"
    )

api.plot_yamamoto(chromsim_results["A"])

# %% [markdown]
# Without the 56.4 mL run the remaining three runs lie exactly on a line (R² = 1.0000) for every
# protein, which confirms that this one run, not the method, causes the scatter.
