"""Calibration workflows recover known parameters from synthetic data (SPEC §9 workflows tier).

Data are simulated with CADET from known "true" parameters (plus noise where stated) and fed to
the workflows as if measured. Tolerances are justified next to each assertion.
"""

import numpy as np
import pytest
from project_examples import make_physics_project, make_project
from synthetic import through_system

from chromsim import api
from chromsim.analysis.peaks import peak_maximum
from chromsim.data.experiment import Signal
from chromsim.estimation.objectives import FitExperiment, SignalTarget
from chromsim.estimation.parameters import FitParameter
from chromsim.estimation.uncertainty import parameter_uncertainty
from chromsim.estimation.workflows.inverse_fit import inverse_fit
from chromsim.estimation.workflows.tracer import TracerRun, analyse_tracers
from chromsim.estimation.workflows.yamamoto import LgeRun, gradient_slope, yamamoto_fit

pytestmark = pytest.mark.slow

# --- tracer ----------------------------------------------------------------------------------

TRUE_MIXER, TRUE_TUBING = 0.30e-6, 0.20e-6  # m³
EPS_C, EPS_P, D_AX = 0.38, 0.60, 1e-7


def tracer_signal(project, kind):
    """Column pulse from CADET, passed through the exact mixer+tubing response (synthetic.py)."""
    result = api.simulate(project)
    outlet = through_system(result.volume_m3, result.outlet("tracer"), TRUE_MIXER, TRUE_TUBING)
    return Signal(kind, "uv", "volume", result.volume_m3, outlet, "AU", "mAU")


def tracer_projects():
    resin = {"name": "r", "particle_radius": "30 um", "col_porosity": EPS_C}
    non_penetrating = make_physics_project("LRM", n_col=500)
    # A pore-excluded tracer behaves like a resin without pore volume (ε_p → 0) in the LRM.
    non_penetrating.resin = resin | {"par_porosity": 1e-9}
    non_penetrating.column.axial_dispersion = f"{D_AX} m^2/s"
    penetrating = make_physics_project("GRM", n_col=100, n_par=10)
    penetrating.resin = resin | {"par_porosity": EPS_P}
    return api.Project.model_validate(non_penetrating.model_dump()), api.Project.model_validate(
        penetrating.model_dump()
    )


def test_tracer_analysis_recovers_porosities_dispersion_and_system():
    non_penetrating, penetrating = tracer_projects()
    volume = np.linspace(0.0, 40e-6, 40001)
    pulse = np.where((volume >= 0.5e-6) & (volume <= 0.6e-6), 1.0, 0.0)  # 0.1 mL from 0.5 mL
    system = Signal(
        "system",
        "uv",
        "volume",
        volume,
        through_system(volume, pulse, TRUE_MIXER, TRUE_TUBING),
        "AU",
        "mAU",
    )
    runs = [
        TracerRun("system", "system", system, 0.5e-6, 0.1e-6),
        TracerRun(
            "dextran", "non_penetrating", tracer_signal(non_penetrating, "dextran"), 0.0, 0.1e-6
        ),
        TracerRun(
            "acetone", "pore_penetrating", tracer_signal(penetrating, "acetone"), 0.0, 0.1e-6
        ),
    ]
    analysis = analyse_tracers(runs, non_penetrating.column, flow_rate_m3_per_s=1e-6 / 60)

    # Moments of exact/converged signals: porosities and volumes within 0.5 %.
    assert analysis.system_volume_m3 == pytest.approx(TRUE_MIXER + TRUE_TUBING, rel=5e-3)
    assert analysis.mixer_volume_m3 == pytest.approx(TRUE_MIXER, rel=1e-2)
    assert analysis.col_porosity == pytest.approx(EPS_C, rel=5e-3)
    assert analysis.par_porosity == pytest.approx(EPS_P, rel=1e-2)
    # D_ax from the variance also contains the solver's numerical dispersion (500 cells, WENO);
    # measured 0.2 % high here. 3 % allows for platform differences without hiding real errors.
    assert analysis.axial_dispersion_m2_per_s == pytest.approx(D_AX, rel=0.03)


# --- Yamamoto ------------------------------------------------------------------------------------


def lge_run(gradient_cv):
    project = make_project()
    project.column.model = "LRM"
    project.column.axial_dispersion = "2e-7 m^2/s"
    project.isotherm.parameters["mAb"] = {
        "k_eq": 1e-2,
        "nu": 6.0,
        "sigma": 30.0,
    }  # rapid equilibrium
    project.method[1].duration = "0.1 mL"
    project.method[1].inlet["mAb"] = "0.1 g/L"  # low load
    project.method[3].duration = f"{gradient_cv} CV"
    project.solver.n_col = 200
    result = api.simulate(project)
    peak_volume, _ = peak_maximum(result.volume_m3, result.outlet("mAb"))
    elution_salt = float(np.interp(peak_volume, result.volume_m3, result.outlet("salt")))
    return LgeRun(
        f"{gradient_cv} CV", gradient_slope(project.resolved_method()[3]), elution_salt
    ), project


def test_yamamoto_recovers_sma_parameters():
    runs, projects = zip(*(lge_run(cv) for cv in (10, 20, 40, 60)), strict=True)
    project = projects[0]
    result = yamamoto_fit(
        list(runs), project.column.volume_m3, project.resin.total_porosity, 1200.0
    )
    # Equilibrium column, low load: the Yamamoto relation is exact up to dispersion and the
    # neglected starting salt; measured: ν 0.06 % and k_eq 0.9 % off.
    assert result.nu == pytest.approx(6.0, rel=5e-3)
    assert result.k_eq == pytest.approx(1e-2, rel=5e-2)
    assert result.r_squared > 0.9999


# --- inverse fit ---------------------------------------------------------------------------------

TRUE = {"k_eq": 1.0, "q_max": 10.0}
TRUE_FEEDS = {"low": 1.05, "high": 1.90}  # mM; nominal 1.0 and 2.0 (load concentration error)


def breakthrough_project(feed_mM):
    project = make_physics_project(
        "LRM",
        components=["protein"],
        isotherm={"model": "Langmuir", "parameters": {"protein": dict(TRUE, k_kin=0.5)}},
        method=[
            {
                "name": "load",
                "duration": "40 min",
                "flow": "1 mL/min",
                "inlet": {"protein": f"{feed_mM} mM"},
            },
            {"name": "wash", "duration": "20 min", "flow": "1 mL/min", "inlet": {}},
        ],
        n_col=40,
    )
    data = project.model_dump(mode="json")
    data["components"] = [{"name": "protein", "mw_kDa": 50.0, "extinction_coeff_L_per_g_cm": 1.0}]
    return api.Project.model_validate(data)


def measured_uv(project, seed):
    """UV of the true model with 0.5 %-of-peak Gaussian noise."""
    result = api.simulate(project)
    uv = api.predict_uv(result, project)
    rng = np.random.default_rng(seed)
    noisy = uv.values + rng.normal(0.0, 0.005 * uv.values.max(), uv.values.size)
    return Signal("UV", "uv", "volume", uv.x[::5], noisy[::5], "AU", "mAU")


def test_inverse_fit_recovers_shared_and_per_experiment_parameters():
    experiments = []
    for seed, (name, feed) in enumerate(TRUE_FEEDS.items()):
        nominal = {"low": 1.0, "high": 2.0}[name]
        experiments.append(
            FitExperiment(
                name,
                breakthrough_project(nominal),
                [SignalTarget(measured_uv(breakthrough_project(feed), seed))],
            )
        )
    parameters = [
        FitParameter(
            path="isotherm.parameters.protein.k_eq",
            lower=0.05,
            upper=20.0,
            scale="log",
            initial=0.3,
            name="k_eq",
        ),
        FitParameter(
            path="isotherm.parameters.protein.q_max",
            lower=1.0,
            upper=50.0,
            initial=6.0,
            name="q_max",
        ),
        FitParameter(
            path="method.0.inlet.protein",
            lower=0.5,
            upper=3.0,
            unit="mM",
            per_experiment=True,
            name="feed",
        ),
    ]
    fit = inverse_fit(experiments, parameters, max_evaluations=100)
    report = parameter_uncertainty(fit)
    table = report.table.set_index("parameter")

    truth = {"k_eq": 1.0, "q_max": 10.0, "feed[low]": 1.05, "feed[high]": 1.90}
    for label, value in truth.items():
        # With 0.5 % noise on ~500 points per run the parameters are determined to well under 1 %;
        # require 2 % and agreement with the reported uncertainty (|error| < 4 standard errors).
        assert fit.values[label] == pytest.approx(value, rel=0.02), label
        assert abs(fit.values[label] - value) < 4 * table.loc[label, "std_error"], label
    assert fit.success
    assert not table["at_bound"].any()
