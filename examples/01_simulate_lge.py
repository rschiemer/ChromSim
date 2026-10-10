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
# # Tutorial 01 — Simulate a linear-gradient elution
#
# This tutorial simulates the example project from SPEC §4.3: a monoclonal antibody (mAb) loaded
# onto a cation-exchange column at 50 mM salt, washed, and eluted with a 20 CV linear gradient
# to 500 mM salt. Binding follows the steric mass action (SMA) isotherm, and the column is
# described by the general rate model (GRM).
#
# You will:
#
# 1. load the project file,
# 2. add the transport parameters the GRM needs,
# 3. simulate with CADET and plot the chromatogram,
# 4. check mass balance and grid resolution.

# %%
from pathlib import Path

import numpy as np

from chromsim import api

EXAMPLES_DIR = Path.cwd() if (Path.cwd() / "projects").is_dir() else Path.cwd() / "examples"
project = api.load_project(EXAMPLES_DIR / "projects" / "cex_mab_lge.yaml")
print(f"{project.name}: {project.column.model} column, isotherm {project.isotherm.model}")
print("Available isotherms:", ", ".join(api.list_isotherms()))

# %% [markdown]
# ## 1. What is missing?
#
# The project file describes the column geometry, resin, method and isotherm, but no transport
# parameters (axial dispersion, film and pore diffusion). These normally come from calibration
# (tracer experiments and inverse fitting, Tutorial 03). ChromSim lists everything that is
# missing instead of guessing:

# %%
try:
    api.simulate(project)
except api.SimulationSetupError as error:
    print(error)

# %% [markdown]
# ## 2. Add transport parameters
#
# Until we have calibrated values, we use typical orders of magnitude for a 45 µm CEX resin:
#
# | parameter | salt | mAb | why |
# |---|---|---|---|
# | axial dispersion D_ax | 2·10⁻⁷ m²/s | (same) | particle Péclet number ≈ 0.5 |
# | film transfer k_f | 1·10⁻⁴ m/s | 2·10⁻⁶ m/s | small ions transfer faster than a mAb |
# | pore diffusion D_p | 1·10⁻¹⁰ m²/s | 5·10⁻¹² m²/s | hindered diffusion of the mAb |
#
# These values are illustrative, not measured.

# %%
project.column.axial_dispersion = "2e-7 m^2/s"
project.column.film_diffusion = {"salt": "1e-4 m/s", "mAb": "2e-6 m/s"}
project.column.pore_diffusion = {"salt": "1e-10 m^2/s", "mAb": "5e-12 m^2/s"}

# %% [markdown]
# ## 3. Simulate
#
# `api.simulate` converts the project to a CADET-Process flow sheet:
#
# inlet → mixer → tubing → column → tubing → outlet
#
# It then runs CADET-Core and returns the outlet concentrations on time, volume and CV axes.

# %%
result = api.simulate(project)
print(f"CADET took {result.solver_stats['wall_time_s']:.1f} s")
print(f"{len(result.time_s)} time points, {result.column_volumes[-1]:.1f} CV in total")

# %%
figure = api.plot_chromatogram(result, x="CV")
figure

# %% [markdown]
# ## 4. Look at the peak
#
# The elution salt concentration at the peak maximum is what the Yamamoto method (Tutorial 03)
# uses to estimate ν and k_eq from several gradient slopes.

# %%
mab = result.outlet("mAb", "g/L")
peak = int(np.argmax(mab))
print(f"peak maximum {mab[peak]:.3f} g/L at {result.column_volumes[peak]:.1f} CV")
print(f"salt at the outlet at the peak: {result.outlet('salt')[peak]:.0f} mM")

# %% [markdown]
# ## 5. Mass balance
#
# All of the loaded mAb (2 mL × 1 g/L = 2 mg) should elute. The flow rate is constant, so the
# eluted mass is the integral of concentration over volume.

# %%
eluted_mg = np.trapezoid(mab, result.volume_m3) * 1e6  # g/L × m³ = kg/m³ × m³ → kg; ×1e6 → mg
print(f"eluted {eluted_mg:.4f} mg of 2 mg loaded ({eluted_mg / 2:.4%})")
assert abs(eluted_mg / 2.0 - 1.0) < 1e-3

# %% [markdown]
# ## 6. Is the grid fine enough?
#
# `api.check_discretisation` simulates again on a doubled grid and reports the normalised RMS
# difference of the outlet profile. Below about 1 %, the grid is fine enough. Here we check a
# coarser grid (50 × 10 cells) to keep the tutorial fast.

# %%
coarse = project.model_copy(deep=True)
coarse.solver.n_col, coarse.solver.n_par = 50, 10
check = api.check_discretisation(coarse)
grids = f"{check.n_col}×{check.n_par} vs {check.refined_n_col}×{check.refined_n_par} cells"
for name, value in check.nrmse.items():
    print(f"{name}: NRMSE {grids} = {value:.2%}")

# %% [markdown]
# ## 7. Reproducibility
#
# Every result carries its provenance: the exact project file, package versions, solver
# settings, time and user.

# %%
print({key: value for key, value in result.provenance.items() if key != "project_yaml"})
