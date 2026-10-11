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
# # Tutorial 03 — Calibrate a model from experiments
#
# Tutorial 02 showed that the first-guess model elutes the mAb too early. Here we calibrate it,
# step by step, the way it is done in the lab:
#
# 1. **Tracer pulses** give the system dead volume and porosities (moments), and the axial
#    dispersion and mixer volume (fit of the dextran pulse).
# 2. **Yamamoto method**: four low-load gradient runs give starting values for the SMA
#    parameters ν and k_eq.
# 3. **Inverse fit** to the preparative run (UV and both fraction assays) refines them.
# 4. **Uncertainty**: confidence intervals and parameter correlations.
#
# All data are synthetic (`tests/data/`, made by `tests/data/generate_fixtures.py`), so we know the
# true mAb parameters (k_eq = 1.5·10⁻², ν = 6). That lets us see what calibration can and cannot
# tell us.

# %%
from pathlib import Path

import numpy as np

from chromsim import api

ROOT = Path.cwd() if (Path.cwd() / "examples").is_dir() else Path.cwd().parent
DATA, EXAMPLES = ROOT / "tests" / "data", ROOT / "examples"
ML = 1e-6  # m³ per mL; ChromSim works in SI internally

project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
calibration = api.load_conductivity_calibration(
    EXAMPLES / "calibrations" / "nacl_conductivity.yaml"
)
flow = project.resolved_method()[0].flow_rate_m3_per_s  # 150 cm/h on this column

# %% [markdown]
# ## 1. Tracer pulses
#
# Each run injects 0.1 mL at the injection mark. We subtract a straight-line UV baseline through
# the signal before the injection and at the end of the run, then analyse the peak moments.

# %%
tracer_runs = []
tracer_files = {
    "system": "tracer_system",
    "non_penetrating": "tracer_dextran",
    "pore_penetrating": "tracer_acetone",
}
for kind, file in tracer_files.items():
    run = api.read_chromatogram(DATA / "akta" / f"{file}.csv")
    injection = run.injections_m3[0]
    uv = api.subtract_baseline(run.signal("uv"), (0.0, injection - 0.05 * ML), (8 * ML, 9 * ML))
    tracer_runs.append(
        api.TracerRun(file, kind, uv, injection_start_m3=injection, injection_volume_m3=0.1 * ML)
    )

tracers = api.analyse_tracers(tracer_runs, project.column, flow)
print(f"system volume (to the UV detector) {tracers.system_volume_m3 / ML:.3f} mL")
eps_c, eps_t, eps_p = tracers.col_porosity, tracers.total_porosity, tracers.par_porosity
print(f"ε_c = {eps_c:.3f}, ε_t = {eps_t:.3f}, ε_p = {eps_p:.3f}")
for note in tracers.notes:
    print("NOTE:", note)

# %% [markdown]
# The **first moments** (system volume and porosities) are robust. The **variances** are not on
# this small column: the mixer broadens the peak much more than the column does, and baseline
# noise in the long tails swamps the difference. The analysis says so (see the note above).
# Instead of trusting moment-based D_ax and mixer volume, we fit a model of the dextran pulse.
# The tracer is pore-excluded, so we use an LRM without pore volume. Its UV response, 400 mAU
# per mM from the injected standard, is entered as an extinction coefficient.

# %%
project = tracers.apply_to(project)  # porosities and system volume from the moments
dextran_project = api.Project.model_validate(
    {
        "name": "dextran pulse",
        "components": [{"name": "dextran", "mw_kDa": 1.0, "extinction_coeff_L_per_g_cm": 2.0}],
        "column": {
            "model": "LRM",
            "length": "20 cm",
            "diameter": "0.5 cm",
            "axial_dispersion": "1e-6 m^2/s",
        },
        "resin": project.resin.model_dump(mode="json") | {"par_porosity": 1e-9},  # no pore volume
        "system": project.system.model_dump(mode="json"),
        "isotherm": {"model": "Linear", "parameters": {}},
        "method": [
            {"name": "equilibrate", "duration": "1 mL", "flow": "150 cm/h", "inlet": {}},
            {
                "name": "inject",
                "duration": "0.1 mL",
                "flow": "150 cm/h",
                "inlet": {"dextran": "1 mM"},
            },
            {"name": "elute", "duration": "7.9 mL", "flow": "150 cm/h", "inlet": {}},
        ],
        "solver": {"n_col": 400},
    }
)
dextran_uv = tracer_runs[1].signal
tracer_fit = api.inverse_fit(
    [
        api.FitExperiment(
            "dextran", dextran_project, [api.SignalTarget(dextran_uv, window=(1 * ML, 6 * ML))]
        )
    ],
    [
        api.FitParameter(
            path="column.axial_dispersion",
            lower=1e-8,
            upper=1e-5,
            scale="log",
            unit="m^2/s",
            name="D_ax",
        ),
        api.FitParameter(
            path="system.mixer_volume", lower=0.01, upper=1.0, unit="mL", name="mixer"
        ),
        api.FitParameter(
            path="system.tubing_volume", lower=0.01, upper=1.0, unit="mL", name="tubing"
        ),
    ],
)
print(
    api.parameter_uncertainty(tracer_fit).table.to_string(
        index=False, float_format=lambda v: f"{v:.4g}"
    )
)

# %% [markdown]
# Mixer and D_ax match the values used to make the data (0.2 mL and 2·10⁻⁷ m²/s, up to the
# numerical dispersion of the discretised model). The tubing comes out 0.27 mL instead of
# 0.25 mL: it absorbs the small error of ε_c from the moments, because both move the peak. We
# write the fitted values into the project.

# %%
project = api.set_values(
    project,
    [
        ("column.axial_dispersion", tracer_fit.values["D_ax"], "m^2/s"),
        ("system.mixer_volume", tracer_fit.values["mixer"], "mL"),
        ("system.tubing_volume", tracer_fit.values["tubing"], "mL"),
    ],
)

# %% [markdown]
# ## 2. Yamamoto method
#
# For each low-load run we need the gradient slope g and the salt concentration I_R at the column
# outlet when the mAb peak elutes. The UV and conductivity detectors sit behind different delay
# volumes, which `elution_salt_at_peak` corrects.

# %%
lge_runs = []
for gradient_cv in (10, 20, 30, 40):
    run = api.read_chromatogram(DATA / "akta" / f"lge_{gradient_cv}cv.csv")
    uv = api.subtract_baseline(run.signal("uv"), (2 * ML, 15 * ML))
    salt = api.conductivity_to_salt(run.signal("conductivity"), calibration)
    elution_salt = api.elution_salt_at_peak(
        uv, salt, project.system.uv_delay_volume_m3, project.system.conductivity_delay_volume_m3
    )
    method = project.model_copy(deep=True)
    method.method[3].duration = f"{gradient_cv} CV"
    slope = api.gradient_slope(method.resolved_method()[3])
    lge_runs.append(api.LgeRun(f"{gradient_cv} CV", slope, elution_salt))
    print(f"{gradient_cv:>3} CV gradient: I_R = {elution_salt:.0f} mM")

yamamoto = api.yamamoto_fit(
    lge_runs,
    project.column.volume_m3,
    project.resin.total_porosity,
    project.resin.ionic_capacity_mol_per_m3,
)
print(f"Yamamoto: ν = {yamamoto.nu:.2f}, k_eq = {yamamoto.k_eq:.2e}, R² = {yamamoto.r_squared:.5f}")
api.plot_yamamoto(yamamoto)

# %% [markdown]
# The Yamamoto relation assumes equilibrium, but our column has film and pore diffusion
# resistance. These delay the peak a little, so ν comes out slightly high and k_eq low. That's
# good enough as a starting point for the inverse fit.

# %% [markdown]
# ## 3. Inverse fit to the preparative run
#
# We fit the mAb's k_eq (on a log scale) and ν to three kinds of data from run 01 at once:
#
# - the baseline-corrected UV trace during elution (normalised residuals),
# - SEC-HPLC monomer and HMW per fraction (weighted by the replicate SD),
# - total protein per fraction (weighted by the SD).
#
# A coarser grid (30 × 6 cells) keeps each simulation at a few seconds. Tutorial 01 shows how to
# check that the grid is fine enough.

# %%
run01 = api.read_chromatogram(DATA / "akta" / "cex_mab_hmw_run01.csv")
uv01 = api.subtract_baseline(run01.signal("uv"), (2 * ML, 15 * ML), (140 * ML, 147 * ML))
sec = api.read_assay(
    DATA / "assays" / "sec_hplc.xlsx", EXAMPLES / "assay_profiles" / "sec_hplc_long.yaml"
)
total = api.read_assay(
    DATA / "assays" / "total_protein_a280.csv",
    EXAMPLES / "assay_profiles" / "total_protein_a280.yaml",
)

start = project.model_copy(deep=True)
start.isotherm.parameters["mAb"]["k_eq"] = yamamoto.k_eq
start.isotherm.parameters["mAb"]["nu"] = yamamoto.nu
start.solver.n_col, start.solver.n_par = 30, 6

experiment = api.FitExperiment(
    "run01",
    start,
    targets=[
        api.SignalTarget(uv01, window=(60 * ML, 140 * ML)),
        api.AssayTarget(sec, fractions=run01.fractions),
        api.AssayTarget(total),
    ],
)
parameters = [
    api.FitParameter(
        path="isotherm.parameters.mAb.k_eq", lower=1e-3, upper=1e-1, scale="log", name="k_eq"
    ),
    api.FitParameter(path="isotherm.parameters.mAb.nu", lower=3.0, upper=10.0, name="nu"),
]


def progress(evaluation, cost):
    if evaluation % 5 == 1:
        print(f"  evaluation {evaluation:>2}: cost {cost:.4g}")


fit = api.inverse_fit([experiment], parameters, callback=progress)
print(f"{fit.message} after {fit.n_evaluations} evaluations")

# %% [markdown]
# ## 4. Uncertainty

# %%
report = api.parameter_uncertainty(fit)
print(report.table.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("highly correlated:", report.high_correlations)
api.plot_correlation(report)

# %% [markdown]
# ## 5. Does the calibrated model reproduce the data?
#
# This is what calibration has to deliver. The first guess was about 25 % off (Tutorial 02).
# The assays have 3–4 % noise.

# %%
fitted_project, fitted_result = fit.projects["run01"], fit.results["run01"]
for assay, fractions in ((sec, run01.fractions), (total, None)):
    table = api.compare_with_assay(fitted_result, assay, fractions=fractions, warn=False).table
    table = table[~table["below_loq"]]
    for analyte, rows in table.groupby("analyte"):
        rms = np.sqrt(
            np.mean(((rows["simulated"] - rows["measured"]) / rows["measured"].max()) ** 2)
        )
        print(f"{analyte:<14} RMS deviation {rms:.1%} of the peak value")
        assert rms < 0.05, f"{analyte} not reproduced within the assay noise"

# %% [markdown]
# ## 6. A lesson about k_eq and ν
#
# Since the data are synthetic, we know the true values: k_eq = 1.5·10⁻² and ν = 6.

# %%
for name, true_value in {"k_eq": 1.5e-2, "nu": 6.0}.items():
    row = report.table.set_index("parameter").loc[name]
    interval = f"{row['ci95_low']:.4g}–{row['ci95_high']:.4g}"
    print(f"{name}: fitted {row['value']:.4g} (95 % CI {interval}), true {true_value:g}")

# %% [markdown]
# The fit reproduces the data at noise level, yet k_eq is about 20 % below the true value and
# outside its "95 %" interval. The two facts are consistent:
#
# - **k_eq and ν are almost perfectly correlated** (|R| > 0.99, see the heat map). Raising ν
#   while lowering k_eq barely changes the elution salt, so the data determine a *combination*
#   of the two very well, and each of them alone poorly.
# - **Confidence intervals only cover measurement noise, not model error.** The data were made
#   with a 100 × 25 grid; we fit on 30 × 6. That shifts the peaks by only ~0.04 %, but along the
#   flat k_eq–ν valley this moves k_eq by tens of percent. Real data have larger model errors.
#
# In practice:
#
# 1. Report correlated parameters together, and use the calibrated model for predictions in the
#    range where it was calibrated.
# 2. Fit on the grid you will predict with, after checking it with `api.check_discretisation`.
# 3. To separate k_eq and ν, add experiments that respond differently to them, e.g. several
#    gradient slopes (the information the Yamamoto method uses) or isocratic runs. Pass them
#    all to `inverse_fit` as separate `FitExperiment`s.

# %% [markdown]
# ## 7. The calibrated model against the data

# %%
comparison = api.compare_with_assay(fitted_result, sec, fractions=run01.fractions)
api.plot_fraction_comparison(comparison, "monomer", result=fitted_result, components=["mAb"])

# %%
uv_fitted = api.predict_uv(fitted_result, fitted_project)
api.plot_overlay(
    api.Experiment("run01", {"UV (measured)": uv01}, fractions=run01.fractions),
    simulated=[uv_fitted],
)
