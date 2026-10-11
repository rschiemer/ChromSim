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
# # Tutorial 02 — Compare a simulation with experimental data
#
# A mAb with a high-molecular-weight impurity (HMW) was loaded onto a CEX column and eluted with
# a salt gradient. We have:
#
# - the ÄKTA chromatogram (UNICORN export with UV, conductivity, %B, pH and fraction marks),
# - total protein per fraction by A280: semicolon-separated text with decimal commas,
# - SEC-HPLC per fraction (monomer and HMW): an Excel file in long format.
#
# We import all three, preprocess the signals and compare them with a simulation that uses
# first-guess isotherm parameters. The data in `tests/data/` are synthetic, but they are
# formatted like real exports.

# %%
from pathlib import Path

import numpy as np

from chromsim import api

ROOT = Path.cwd() if (Path.cwd() / "examples").is_dir() else Path.cwd().parent
DATA = ROOT / "tests" / "data"
EXAMPLES = ROOT / "examples"

project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
result = api.simulate(project)
print(f"simulated {project.name} in {result.solver_stats['wall_time_s']:.1f} s")

# %% [markdown]
# ## 1. Import the chromatogram
#
# `read_chromatogram` recognises the UNICORN export by itself. It handles UTF-16 text,
# tab/comma/semicolon separators and decimal commas, and it never modifies the raw file. All
# values are converted to SI units; each signal remembers its original unit for display.

# %%
experiment = api.read_chromatogram(DATA / "akta" / "cex_mab_hmw_run01.csv")
for signal in experiment.signals.values():
    points = len(signal.x)
    unit = signal.original_unit
    print(f"{signal.name:<10} {signal.kind:<16} {points:>5} points, exported in {unit}")
first, last = experiment.fractions[0].label, experiment.fractions[-1].label
print(f"{len(experiment.fractions)} fractions, {first}–{last}")
print(f"injection at {experiment.injections_m3[0] * 1e6:.2f} mL")

# %% [markdown]
# ## 2. Preprocess
#
# - **UV:** remove the baseline offset and drift, using two protein-free windows (equilibration
#   and the end of the strip).
# - **Conductivity → salt**, using the buffer's conductivity calibration (it is not linear).

# %%
ML = 1e-6
uv = api.subtract_baseline(experiment.signal("uv"), (2 * ML, 15 * ML), (140 * ML, 147 * ML))
calibration = api.load_conductivity_calibration(
    EXAMPLES / "calibrations" / "nacl_conductivity.yaml"
)
salt = api.conductivity_to_salt(experiment.signal("conductivity"), calibration)

# %% [markdown]
# ## 3. Overlay measured and simulated detector signals
#
# The simulation is converted to what the detectors would show: UV by Lambert–Beer, and
# conductivity through the same calibration, each behind its delay volume.

# %%
uv_sim = api.predict_uv(result, project)
cond_sim = api.predict_conductivity(result, project, calibration)
corrected_experiment = api.Experiment(
    name="run01 (baseline-corrected UV)",
    signals={"UV (measured)": uv, "Cond (measured)": experiment.signal("conductivity")},
    fractions=experiment.fractions,
)
figure = api.plot_overlay(corrected_experiment, simulated=[uv_sim, cond_sim])
figure

# %% [markdown]
# The salt gradient matches well, which confirms the method, the system volumes and the
# conductivity calibration. The protein peak does not: the measured peak elutes later than
# simulated.


# %%
def peak_ml(signal):
    return signal.x[np.argmax(signal.values)] / ML


elution = (salt.x > 60 * ML) & (salt.x < 120 * ML)
salt_sim = np.interp(
    salt.x[elution],
    result.volume_m3 + project.system.conductivity_delay_volume_m3,
    result.outlet("salt"),
)
salt_error = np.max(np.abs(salt.values[elution] - salt_sim))
print(f"salt during the gradient: max |measured − simulated| = {salt_error:.1f} mM")
print(f"UV peak: measured {peak_ml(uv):.1f} mL, simulated {peak_ml(uv_sim):.1f} mL")

# %% [markdown]
# ## 4. Import the fraction assays
#
# Each lab layout is described once in a *mapping profile* (`examples/assay_profiles/`). The two
# files here differ in separator, decimal mark, layout (wide vs long), file type and in whether
# they contain fraction volumes.

# %%
total = api.read_assay(
    DATA / "assays" / "total_protein_a280.csv",
    EXAMPLES / "assay_profiles" / "total_protein_a280.yaml",
)
sec = api.read_assay(
    DATA / "assays" / "sec_hplc.xlsx", EXAMPLES / "assay_profiles" / "sec_hplc_long.yaml"
)
print(total.results.head(5).to_string(index=False))
print()
print(sec.results.head(6).to_string(index=False))

# %% [markdown]
# ## 5. Fraction-averaged comparison
#
# The simulated outlet is averaged over each fraction (∫c dV / ΔV) and converted to the assay's
# unit. The A280 file contains its own fraction volumes. The SEC file does not, so its sample ids
# are matched to the ÄKTA fraction marks. Ids that don't match are **reported**, never guessed:

# %%
total_comparison = api.compare_with_assay(result, total)
sec_comparison = api.compare_with_assay(result, sec, fractions=experiment.fractions)
print("SEC samples without a collected fraction:", sec_comparison.assay_fractions_not_collected)
print("collected fractions without SEC result:", sec_comparison.collected_fractions_without_assay)

# %%
figure = api.plot_fraction_comparison(
    sec_comparison, "monomer", result=result, components=sec.readouts["monomer"]
)
figure

# %%
figure = api.plot_fraction_comparison(
    total_comparison, "total_protein", result=result, components=total.readouts["total_protein"]
)
figure

# %% [markdown]
# ## 6. How far off is the first guess?

# %%
for name, comparison in (("total protein", total_comparison), ("SEC", sec_comparison)):
    table = comparison.table[~comparison.table["below_loq"]]
    for analyte, rows in table.groupby("analyte"):
        error = (rows["simulated"] - rows["measured"]) / rows["measured"].max()
        rms = np.sqrt(np.mean(error**2))
        print(f"{name:>13} / {analyte:<13} RMS deviation {rms:.1%} of the peak value")

# %% [markdown]
# The HMW readout agrees to within the assay noise (3–4 %). Its first-guess parameters happen to
# be right, which also shows that import, fraction matching and averaging work end to end. The
# mAb is off by about 25 % of the peak and elutes early: its first-guess `k_eq` is too small.
# Tutorial 03 estimates the isotherm parameters from these data.
