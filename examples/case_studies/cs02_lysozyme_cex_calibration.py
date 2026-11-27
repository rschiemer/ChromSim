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
# # Case study 02 — Lysozyme on a strong cation exchanger, calibrated from measured data
#
# Real measurements by S. Menzel and J. Buyel (Fraunhofer IME, Aachen), published with
# CADET-MCMC (CC BY-NC-SA 4.0; non-commercial use, attribution). Fetch them first:
#
#     pixi run python examples/case_studies/fetch_data.py fraunhofer_lysozyme
#
# The set is complete for a full model calibration:
#
# | Files | Experiment | Used for |
# |---|---|---|
# | `bypass{Acetone,Dextran}{1,2,3}.dat` | tracer pulses, column bypassed | system volume |
# | `columnDextran{1,2,3}.dat` | pore-excluded dextran pulses | ε_c, axial dispersion |
# | `columnAcetone{1,2,3}.dat` | pore-penetrating acetone pulses | ε_t, ε_p |
# | `columnLysozyme{0,5,10,30,60,120}cv.dat` | lysozyme load–wash–gradient elution | SMA ν, k_eq |
#
# What we know about the runs (from the CADET-MCMC set-up files): column 2.5 cm × 0.7 cm i.d.
# (0.962 mL), particles 45 µm, flow 0.5 mL/min; tracer injections 12 s (0.1 mL); lysozyme runs
# load 0.2 mM lysozyme for 12.6 s at 20 mM salt, wash until 591 s, then a linear gradient to
# 526 mM salt and a 526 mM strip. The ionic capacity was titrated: 10 mM × 19.252 mL per column.
# Every file holds two columns, time (s) and outlet concentration (mM), **no header, no
# conductivity** — so this case study also shows how to work with data that are not ÄKTA exports.
#
# Workflow: (1) tracer moments, (2) fit of the dextran pulse, (3) Yamamoto method, (4) inverse
# fit to three gradients, (5) validation on two gradients and the step elution not used in
# the fit.

# %%
from pathlib import Path

from chromsim import api

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
DATA = HERE / "data" / "fraunhofer_lysozyme"
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)
if not (DATA / "columnLysozyme30cv.dat").is_file():
    raise SystemExit(f"Data not found in {DATA}; run examples/case_studies/fetch_data.py first.")

ML = 1e-6  # m³ per mL; ChromSim works in SI internally
FLOW = 0.5 * ML / 60.0  # 0.5 mL/min in m³/s
INJECTION = 12.0 * FLOW  # 12 s tracer injections
COLUMN = api.Column.model_validate({"model": "GRM", "length": "2.5 cm", "diameter": "0.7 cm"})


def read_run(name: str, signal_name: str = "c"):
    """One measured run: time (s) and concentration (mM), put on the pumped-volume axis."""
    mapping = {
        "header": False,
        "delimiter": "whitespace",
        "x_column": 0,
        "x_unit": "s",
        "signals": [{"column": 1, "name": signal_name, "unit": "mM", "kind": "concentration"}],
    }
    experiment = api.read_chromatogram(DATA / f"{name}.dat", "generic_csv", mapping=mapping)
    return api.time_to_volume(experiment.signals[signal_name], FLOW)


# %% [markdown]
# ## 1. Tracer moments
#
# Two things matter with these data:
#
# - **Each tracer needs its own bypass run.** The dextran bypass peak is centred at 0.32 mL,
#   the acetone bypass peak at 0.24 mL. Large, slowly diffusing dextran is spread much more by
#   laminar flow in the tubing (Taylor dispersion) than acetone, and its long tail shifts the
#   first moment. Subtracting the acetone system volume from the dextran peak would give an ε_c
#   that is 10 % too high.
# - **Baselines decide the moments.** An offset of 1 % of the peak height, integrated over the
#   400 s tail, moves the first moment by tens of µL. We draw a straight baseline through the
#   10 s before the injection and the last minute of each run, and compute moments only up to
#   300 s (where every peak is back at baseline).


# %%
def tracer_run(name: str, kind: str):
    signal = api.subtract_baseline(read_run(name), (0.0, 10 * FLOW), (480 * FLOW, 540 * FLOW))
    return api.TracerRun(name, kind, signal, 0.0, INJECTION, window=(0.0, 300 * FLOW))


tracer_results = {}
for tracer, kind in (("Dextran", "non_penetrating"), ("Acetone", "pore_penetrating")):
    runs = [tracer_run(f"bypass{tracer}{i}", "system") for i in (1, 2, 3)]
    runs += [tracer_run(f"column{tracer}{i}", kind) for i in (1, 2, 3)]
    tracer_results[tracer] = api.analyse_tracers(runs, COLUMN, FLOW)
    for peak in tracer_results[tracer].peaks:
        print(f"{peak.name:<15} first moment {peak.first_moment_m3 / ML:6.3f} mL")

eps_c = tracer_results["Dextran"].col_porosity
eps_t = tracer_results["Acetone"].total_porosity
eps_p = (eps_t - eps_c) / (1.0 - eps_c)
system_volume_salt = tracer_results["Acetone"].system_volume_m3  # small molecules like salt
print(f"\nε_c = {eps_c:.3f}, ε_t = {eps_t:.3f}, ε_p = {eps_p:.3f}")
print(
    f"system volume: acetone {system_volume_salt / ML:.3f} mL, "
    f"dextran {tracer_results['Dextran'].system_volume_m3 / ML:.3f} mL"
)
assert 0.38 < eps_c < 0.45 and 0.90 < eps_t < 0.95, "porosities outside the expected range"

# %% [markdown]
# The replicates agree to about 1 %. The moment analysis also warns that the column contributes
# only about a third of the dextran peak variance, so D_ax from moments is unreliable. As in
# Tutorial 03 we fit a model of the dextran pulse instead: a pore-excluded tracer through
# mixer, tubing and an LRM column. The measured signal is a concentration, so it is compared
# directly with the simulated outlet (`SignalTarget(..., component="dextran")`).

# %%
RESIN = {
    "name": "strong cation exchanger (Fraunhofer IME)",
    "particle_radius": "45 um",
    # titrated capacity per solid volume: 10 mM × 19.252 mL / (V_c (1 − ε_t))
    "ionic_capacity": f"{10.0 * 19.252 / (COLUMN.volume_m3 / ML * (1.0 - eps_t)):.1f} mM",
    "col_porosity": eps_c,
    "par_porosity": eps_p,
}
dextran_project = api.Project.model_validate(
    {
        "name": "dextran pulse",
        "components": [{"name": "dextran"}],
        "column": {
            "model": "LRM",
            "length": "2.5 cm",
            "diameter": "0.7 cm",
            "axial_dispersion": "5e-7 m^2/s",
        },
        "resin": RESIN | {"par_porosity": 1e-9},  # dextran does not enter the pores
        "system": {"mixer_volume": "0.2 mL", "tubing_volume": "0.1 mL"},
        "isotherm": {"model": "Linear", "parameters": {}},
        "method": [
            {
                "name": "inject",
                "duration": "12 s",
                "flow": "0.5 mL/min",
                "inlet": {"dextran": "2.7e-3 mM"},
            },
            {"name": "elute", "duration": "528 s", "flow": "0.5 mL/min", "inlet": {}},
        ],
        "solver": {"n_col": 200},
    }
)
dextran_targets = [
    api.SignalTarget(
        tracer_run(f"columnDextran{i}", "non_penetrating").signal,
        window=(0.0, 300 * FLOW),
        component="dextran",
    )
    for i in (1, 2, 3)
]
dextran_fit = api.inverse_fit(
    [api.FitExperiment("dextran", dextran_project, dextran_targets)],
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
            path="system.mixer_volume", lower=0.01, upper=0.5, unit="mL", name="mixer"
        ),
        api.FitParameter(
            path="system.tubing_volume", lower=0.0, upper=0.5, unit="mL", name="tubing"
        ),
    ],
)
print(
    api.parameter_uncertainty(dextran_fit).table.to_string(
        index=False, float_format=lambda v: f"{v:.4g}"
    )
)

# %% [markdown]
# ## 2. The lysozyme project
#
# One project per gradient: the files are named by nominal gradient length, but the gradients
# actually ran 6.0, 10.4, 30, 60 and 124 CV (from the section times), and we use the real ones.
# Salt and lysozyme use the measured system; the protein binds with the SMA isotherm in
# rapid equilibrium (no `k_kin`). Transport parameters are typical values for a 45 µm resin; the
# fit below refines the lysozyme pore diffusion.

# %%
GRADIENT_END_S = {0: 591.0, 5: 1284.0, 10: 1795.0, 30: 4055.0, 60: 7518.0, 120: 14950.0}
WASH_END_S, LOAD_S = 591.0, 12.6


def lysozyme_project(nominal_cv: int, isotherm=None, pore_diffusion: str = "1e-11 m^2/s"):
    """Project of the run with the given nominal gradient (0 = step elution)."""
    elution = (
        [{"name": "step", "duration": "2000 s", "flow": "0.5 mL/min", "inlet": {"salt": "526 mM"}}]
        if nominal_cv == 0
        else [
            {
                "name": "gradient",
                "duration": f"{GRADIENT_END_S[nominal_cv] - WASH_END_S} s",
                "flow": "0.5 mL/min",
                "inlet": {"salt": "20 mM"},
                "inlet_end": {"salt": "526 mM"},
            },
            {
                "name": "strip",
                "duration": "1500 s",
                "flow": "0.5 mL/min",
                "inlet": {"salt": "526 mM"},
            },
        ]
    )
    return api.Project.model_validate(
        {
            "name": f"lysozyme {nominal_cv} CV (nominal)",
            "components": [{"name": "salt"}, {"name": "lysozyme", "mw_kDa": 14.3}],
            "column": {
                "model": "GRM",
                "length": "2.5 cm",
                "diameter": "0.7 cm",
                "axial_dispersion": f"{dextran_fit.values['D_ax']} m^2/s",
                "film_diffusion": {"salt": "1e-4 m/s", "lysozyme": "2e-5 m/s"},
                "pore_diffusion": {"salt": "1e-9 m^2/s", "lysozyme": pore_diffusion},
            },
            "resin": RESIN,
            "system": {
                "mixer_volume": f"{dextran_fit.values['mixer']} mL",
                "tubing_volume": f"{dextran_fit.values['tubing']} mL",
            },
            "isotherm": {
                "model": "SMA",
                "parameters": {"lysozyme": isotherm or {"k_eq": 1e-2, "nu": 4.0, "sigma": 10.0}},
            },
            "method": [
                {
                    "name": "load",
                    "duration": f"{LOAD_S} s",
                    "flow": "0.5 mL/min",
                    "inlet": {"salt": "20 mM", "lysozyme": "0.2 mM"},
                },
                {
                    "name": "wash",
                    "duration": f"{WASH_END_S - LOAD_S} s",
                    "flow": "0.5 mL/min",
                    "inlet": {"salt": "20 mM"},
                },
                *elution,
            ],
            "solver": {"n_col": 40, "n_par": 8},
        }
    )


measured = {cv: read_run(f"columnLysozyme{cv}cv", "lysozyme") for cv in GRADIENT_END_S}

# %% [markdown]
# ## 3. Yamamoto method
#
# No conductivity was recorded, so the outlet salt is the programmed inlet salt shifted by the
# volume a salt front needs to reach the detector: the (acetone) system volume plus the
# pore-accessible column volume ε_t·V_c (`api.inlet_profile`).

# %%
salt_delay = system_volume_salt + eps_t * COLUMN.volume_m3
lge_runs = []
for cv in (5, 10, 30, 60, 120):
    project = lysozyme_project(cv)
    salt = api.inlet_profile(project, "salt", delay_volume_m3=salt_delay)
    elution_salt = api.elution_salt_at_peak(measured[cv], salt)
    gradient = project.resolved_method()[2]
    actual_cv = gradient.flow_rate_m3_per_s * gradient.duration_s / COLUMN.volume_m3
    lge_runs.append(api.LgeRun(f"{actual_cv:.1f} CV", api.gradient_slope(gradient), elution_salt))
    print(f"{actual_cv:6.1f} CV gradient: I_R = {elution_salt:5.1f} mM")

yamamoto = api.yamamoto_fit(
    lge_runs, COLUMN.volume_m3, eps_t, lysozyme_project(30).resin.ionic_capacity_mol_per_m3
)
print(f"Yamamoto: ν = {yamamoto.nu:.2f}, k_eq = {yamamoto.k_eq:.3g}, R² = {yamamoto.r_squared:.4f}")
api.plot_yamamoto(yamamoto).write_html(OUTPUT / "cs02_yamamoto.html")
assert yamamoto.r_squared > 0.99

# %% [markdown]
# Five gradients over a 20-fold range of slopes lie on one straight line (R² > 0.99), as the
# Yamamoto relation predicts for a protein binding by charge.
#
# ## 4. Inverse fit to three gradients
#
# Starting from the Yamamoto values, we fit k_eq, ν and the lysozyme pore diffusion to the 10,
# 30 and 60 CV runs at once (only the elution peaks; normalised residuals). The 6 and 124 CV
# gradients and the step elution are kept back for validation.

# %%
start_isotherm = {"k_eq": yamamoto.k_eq, "nu": yamamoto.nu, "sigma": 10.0}
fit_runs = (10, 30, 60)
experiments = []
for cv in fit_runs:
    peak_volume = api.peak_maximum(measured[cv].x, measured[cv].values)[0]
    window = (peak_volume - 4 * ML, peak_volume + 6 * ML)
    experiments.append(
        api.FitExperiment(
            f"{cv} CV",
            lysozyme_project(cv, start_isotherm),
            [api.SignalTarget(measured[cv], window=window, component="lysozyme")],
        )
    )
parameters = [
    api.FitParameter(
        path="isotherm.parameters.lysozyme.k_eq", lower=1e-4, upper=10.0, scale="log", name="k_eq"
    ),
    api.FitParameter(path="isotherm.parameters.lysozyme.nu", lower=1.0, upper=10.0, name="nu"),
    # Upper bound above free diffusion (~1e-10 m²/s for lysozyme) on purpose, as in the
    # original CADET-MCMC study: see the discussion of D_p below.
    api.FitParameter(
        path="column.pore_diffusion.lysozyme",
        lower=1e-12,
        upper=1e-8,
        scale="log",
        unit="m^2/s",
        name="D_p",
    ),
]
fit = api.inverse_fit(experiments, parameters, max_evaluations=80)
report = api.parameter_uncertainty(fit)
print(f"{fit.message} after {fit.n_evaluations} evaluations")
print(report.table.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("highly correlated:", report.high_correlations)

# %% [markdown]
# **About D_p.** D_p runs to its upper bound, two orders of magnitude above the free
# diffusivity of lysozyme in water (about 1·10⁻¹⁰ m²/s), which pore diffusion cannot exceed.
# The measured peaks are narrower than the general rate model can make them with physical pore
# diffusion; refining the grid (120 × 16 instead of 40 × 8) changes the peak width by < 2 %,
# so this is not numerical. Once intraparticle transport no longer limits the peak width, the
# data only say "D_p is large", hence the bound. The data ask for an extra, fast intraparticle
# transport path; surface diffusion of bound protein is the usual explanation. Treat D_p here
# as an *effective* parameter for these conditions, not a material property. (The k_eq–ν
# result and the peak positions do not depend on it.)
#
# ## 5. Validation on runs not used in the fit
#
# The calibrated model predicts the 6 and 124 CV gradients and the step elution. We compare
# peak positions, and plot measured (solid) against simulated (dashed) for all six runs.

# %%
fitted_isotherm = {"k_eq": fit.values["k_eq"], "nu": fit.values["nu"], "sigma": 10.0}
fitted_dp = f"{fit.values['D_p']} m^2/s"
print(f"{'run':>10} {'peak measured':>14} {'peak simulated':>15} {'deviation':>10}")
simulated_signals, deviations = [], {}
for cv in GRADIENT_END_S:
    project = lysozyme_project(cv, fitted_isotherm, fitted_dp)
    result = api.simulate(project)
    simulated = api.Signal(
        f"{cv} CV (simulated)",
        "concentration",
        "volume",
        result.volume_m3 + project.system.uv_delay_volume_m3,
        result.outlet("lysozyme"),
        "mol/m**3",
        "mM",
    )
    simulated_signals.append(simulated)
    peak_measured = api.peak_maximum(measured[cv].x, measured[cv].values)[0]
    peak_simulated = api.peak_maximum(simulated.x, simulated.values)[0]
    deviations[cv] = (peak_simulated - peak_measured) / (peak_measured - WASH_END_S * FLOW)
    role = "fit" if cv in fit_runs else "prediction"
    print(
        f"{cv:>7} CV {peak_measured / ML:11.2f} mL {peak_simulated / ML:12.2f} mL "
        f"{deviations[cv]:9.1%}  ({role})"
    )

figure = api.plot_overlay(
    api.Experiment("lysozyme", {f"{cv} CV": measured[cv] for cv in GRADIENT_END_S}),
    simulated=simulated_signals,
    kinds=("concentration",),
)
figure.write_html(OUTPUT / "cs02_lysozyme_overlay.html")
gradient_runs = [cv for cv in GRADIENT_END_S if cv != 0]
assert max(abs(deviations[cv]) for cv in gradient_runs) < 0.05, "gradient peaks not reproduced"

# %% [markdown]
# **Gradients.** The peak positions of all five gradient runs, including the 6 and 124 CV
# predictions, agree with the measurements within about 2 % of the elution volume. As in
# Tutorial 03, k_eq and ν are strongly correlated; the data pin down their combination (the
# elution salt), which is what predictions need.
#
# **Step elution: an open discrepancy.** The model elutes the step run about 1.5 mL early. The
# measured peak is sharp (0.8 mL wide) but starts 1.8 mL after the salt front should leave the
# column, so it is shifted rather than broadened, which points to the step having started
# later than the 591 s in the set-up files rather than to a model error. Check such timing
# assumptions against the run log when you have it; these published files do not include one.
#
# Figures are in `examples/case_studies/output/`.
