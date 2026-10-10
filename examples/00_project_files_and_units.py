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
# # 00 — Project files and units
#
# A ChromSim **project** describes one chromatography run: components, column, resin, system,
# isotherm, method and solver settings. This notebook loads the example project from SPEC §4.3,
# looks at how user units become SI units, edits the method and saves it again.

# %%
import tempfile
from pathlib import Path

from pydantic import ValidationError

from chromsim import api, units

# Works both from the repository root and from inside examples/.
EXAMPLES_DIR = Path.cwd() if (Path.cwd() / "projects").is_dir() else Path.cwd() / "examples"
project = api.load_project(EXAMPLES_DIR / "projects" / "cex_mab_lge.yaml")
print(project.name, "with components", project.component_names)

# %% [markdown]
# ## Column geometry in SI
#
# Users write `20 cm` and `0.5 cm`; the model gives SI values as properties with the unit in
# the name. One column volume is V_c = π·(d/2)²·L.

# %%
column = project.column
print(f"length        {column.length}  ->  {column.length_m} m")
print(f"cross-section {column.cross_section_m2:.4e} m²")
print(f"1 CV          {column.volume_m3:.4e} m³ = {units.from_si(column.volume_m3, 'mL'):.3f} mL")
print(f"total porosity ε_t = {project.resin.total_porosity:.3f}")

# %% [markdown]
# ## Flow: cm/h ↔ mL/min
#
# Linear velocity is the superficial velocity, so Q = u·A.

# %%
for flow in ["150 cm/h", "300 cm/h"]:
    converted = units.convert_flow(api.Quantity.parse(flow), "mL/min", column.cross_section_m2)
    print(f"{flow:>9} on a {column.diameter} column = {converted.magnitude:.3f} mL/min")

# %% [markdown]
# ## The method as an SI timeline
#
# Step durations in CV or mL become seconds using the step's flow rate, and concentrations
# (mM or g/L) become mol/m³ in project component order (salt, mAb).

# %%
print(f"{'step':<12}{'start / s':>10}{'end / s':>10}   inlet start → end / (mol/m³)")
for step in project.resolved_method():
    start = ", ".join(f"{c:.4g}" for c in step.inlet_start_mol_per_m3)
    end = ", ".join(f"{c:.4g}" for c in step.inlet_end_mol_per_m3)
    print(f"{step.name:<12}{step.start_s:>10.1f}{step.end_s:>10.1f}   ({start}) → ({end})")

# %% [markdown]
# ## Edit, validate and save
#
# Changes are validated straight away; a wrong unit is rejected with a message saying how to fix it.

# %%
project.method[3].duration = "30 CV"  # a shallower gradient

try:
    project.column.length = "20 mL"
except ValidationError as error:
    print("Rejected:", error.errors()[0]["msg"])

with tempfile.TemporaryDirectory() as folder:
    path = api.save_project(project, Path(folder) / "cex_mab_lge_30cv.yaml")
    print(path.read_text(encoding="utf-8")[:400], "...")
    assert api.load_project(path) == project  # lossless round trip
