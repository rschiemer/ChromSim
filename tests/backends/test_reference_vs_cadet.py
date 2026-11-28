"""Reference backend against CADET on shared benchmark cases (SPEC §5b acceptance).

Requirements: outlet profiles within 1 % NRMSE of CADET, mass balance closure ≤ 0.1 %. Both
backends solve the same equations on the same grid (100 cells); the remaining differences are the
spatial schemes (van Leer vs WENO) and time stepping. Measured NRMSE: ≤ 0.4 % (tracers),
≤ 0.13 % (with binding).
"""

import numpy as np
import pytest
from project_examples import make_binding_project, make_physics_project

from chromsim import api
from chromsim.backends.checks import nrmse
from chromsim.errors import SimulationSetupError

pytestmark = pytest.mark.slow

NRMSE_LIMIT = 0.01  # SPEC §5b
MASS_BALANCE_LIMIT = 1e-3  # SPEC §5b

SALT = [{"name": "salt"}]
LANGMUIR_2 = {
    "model": "Langmuir",
    "parameters": {
        "A": {"k_eq": 1.0, "q_max": 10.0, "k_kin": 1.0},
        "B": {"k_eq": 0.3, "q_max": 8.0, "k_kin": 2.0},
    },
}
SMA = {"model": "SMA", "parameters": {"A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0}}}
MPM_PARAMETERS = {"A": {"k_eq": 1e3, "q_max": 10.0, "beta": 2.0, "gamma": 1e-3, "k_kin": 1.0}}


def with_model(project, column_model, backend="cadet", system=None):
    data = project.model_dump(mode="json")
    data["column"]["model"] = column_model
    data["solver"] = data["solver"] | {"n_col": 100, "backend": backend}
    if system:
        data["system"] = system
    return api.Project.model_validate(data)


def benchmark_cases():
    tracer = make_physics_project("LRM")
    linear = make_binding_project(
        [{"name": "A"}],
        {"model": "Linear", "parameters": {"A": {"k_eq": 2.0, "k_kin": 1.0}}},
        {"A": "1 mM"},
    )
    langmuir = make_binding_project(
        [{"name": "A"}, {"name": "B"}], LANGMUIR_2, {"A": "1 mM", "B": "1 mM"}
    )
    sma = make_binding_project(
        SALT + [{"name": "A"}], SMA, {"salt": "50 mM", "A": "0.5 mM"}, elute={"salt": "1000 mM"}
    )
    cases = []
    for model in ("LRM", "LRMP"):
        cases += [
            pytest.param(tracer, model, None, ["tracer"], id=f"tracer-{model}"),
            pytest.param(linear, model, None, ["A"], id=f"linear-{model}"),
            pytest.param(langmuir, model, None, ["A", "B"], id=f"langmuir-{model}"),
            pytest.param(sma, model, None, ["salt", "A"], id=f"sma-{model}"),
        ]
    system = {"mixer_volume": "0.3 mL", "tubing_volume": "0.2 mL"}
    cases.append(pytest.param(sma, "LRM", system, ["salt", "A"], id="sma-LRM-with-system"))
    return cases


@pytest.mark.parametrize(("project", "column_model", "system", "components"), benchmark_cases())
def test_reference_matches_cadet(project, column_model, system, components):
    cadet = api.simulate(with_model(project, column_model, "cadet", system))
    reference = api.simulate(with_model(project, column_model, "reference", system))
    for component in components:
        on_cadet_grid = np.interp(cadet.time_s, reference.time_s, reference.outlet(component))
        assert nrmse(on_cadet_grid, cadet.outlet(component)) < NRMSE_LIMIT, component


@pytest.mark.parametrize("column_model", ["LRM", "LRMP"])
def test_python_isotherm_matches_its_cadet_twin(column_model):
    # MixedModeLangmuir (Python, reference backend) has the same equations as CADET's MPM.
    common = {
        "components": SALT + [{"name": "A"}],
        "load": {"salt": "50 mM", "A": "0.1 mM"},
        "elute": {"salt": "1000 mM"},
        "column_model": column_model,
        "n_col": 100,
    }
    python = make_binding_project(
        isotherm={"model": "MixedModeLangmuir", "parameters": MPM_PARAMETERS}, **common
    )
    cadet = make_binding_project(isotherm={"model": "MPM", "parameters": MPM_PARAMETERS}, **common)
    python_result = api.simulate(with_model(python, column_model, "reference"))
    cadet_result = api.simulate(with_model(cadet, column_model, "cadet"))
    on_cadet_grid = np.interp(cadet_result.time_s, python_result.time_s, python_result.outlet("A"))
    assert nrmse(on_cadet_grid, cadet_result.outlet("A")) < NRMSE_LIMIT


COLLOIDAL_2 = {
    "model": "Colloidal",
    "parameters": {
        "A": {
            "radius": 4.5e-9,
            "ln_ke_salt_factor": 40.0,
            "ln_ke_salt_exponent": 0.4,
            "ln_ke_exp_factor": 1.0,
            "ln_ke_exp_salt_coefficient": -1e-3,
            "b_pp_salt_factor": 5.0,
            "b_pp_salt_exponent": 0.2,
            "b_pp_exp_factor": 2.0,
            "b_pp_exp_salt_coefficient": -2e-3,
            "k_kin": 1.0,
        },
        "B": {
            "radius": 3e-9,
            "ln_ke_salt_factor": 30.0,
            "ln_ke_salt_exponent": 0.4,
            "b_pp_salt_factor": 3.0,
            "b_pp_salt_exponent": 0.0,
            "k_kin": 2.0,
        },
    },
    "global_parameters": {
        "phase_ratio": 5e7,
        "kappa_factor": 1e-10,
        "kappa_exponent": 0.3,
        "kappa_constant": 5e-10,
    },
}


@pytest.mark.parametrize("column_model", ["LRM", "LRMP"])
def test_colloidal_on_reference_backend_matches_cadet(column_model):
    # Two proteins of different size, every salt term switched on, and a load where the lattice
    # term matters (dropping it changes the outlet by 16-29 % NRMSE): this pins the mapping onto
    # CADET's COL_* names, which CADET's own docs describe inconsistently (ADR 0011). Measured
    # NRMSE <= 0.3 %. Much higher loads (q beyond close packing) are too stiff for the
    # reference backend's finite-difference Jacobian.
    project = make_binding_project(
        SALT + [{"name": "A"}, {"name": "B"}],
        COLLOIDAL_2,
        {"salt": "50 mM", "A": "0.1 mM", "B": "0.1 mM"},
        elute={"salt": "1000 mM"},
        column_model=column_model,
        n_col=100,
    )
    cadet = api.simulate(with_model(project, column_model, "cadet"))
    reference = api.simulate(with_model(project, column_model, "reference"))
    for component in ("A", "B"):
        on_cadet_grid = np.interp(cadet.time_s, reference.time_s, reference.outlet(component))
        assert nrmse(on_cadet_grid, cadet.outlet(component)) < NRMSE_LIMIT, component


@pytest.mark.parametrize("column_model", ["LRM", "LRMP"])
def test_sas_on_reference_backend_matches_cadet(column_model):
    # Strong self-association (k_sa c ≈ 5 at the load concentration): pins k_a2 = k_sa·k_a.
    sas = {
        "model": "SAS",
        "parameters": {
            "A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "self_association": 10.0, "k_kin": 10.0}
        },
    }
    project = make_binding_project(
        SALT + [{"name": "A"}],
        sas,
        {"salt": "50 mM", "A": "0.5 mM"},
        elute={"salt": "1000 mM"},
        column_model=column_model,
        n_col=100,
    )
    cadet = api.simulate(with_model(project, column_model, "cadet"))
    reference = api.simulate(with_model(project, column_model, "reference"))
    on_cadet_grid = np.interp(cadet.time_s, reference.time_s, reference.outlet("A"))
    assert nrmse(on_cadet_grid, cadet.outlet("A")) < NRMSE_LIMIT


@pytest.mark.parametrize("column_model", ["LRM", "LRMP"])
def test_reference_mass_balance(column_model):
    tracer = api.simulate(with_model(make_physics_project("LRM"), column_model, "reference"))
    assert np.trapezoid(tracer.outlet("tracer"), tracer.time_s) == pytest.approx(
        6.0, rel=MASS_BALANCE_LIMIT
    )
    project = with_model(
        make_binding_project(
            SALT + [{"name": "A"}], SMA, {"salt": "50 mM", "A": "0.5 mM"}, elute={"salt": "1000 mM"}
        ),
        column_model,
        "reference",
    )
    result = api.simulate(project)
    load = project.resolved_method()[1]
    injected = load.inlet_start_mol_per_m3[1] * load.duration_s
    assert np.trapezoid(result.outlet("A"), result.time_s) == pytest.approx(
        injected, rel=MASS_BALANCE_LIMIT
    )


def test_reference_backend_refuses_what_it_cannot_do():
    sma = make_binding_project(SALT + [{"name": "A"}], SMA, {"salt": "50 mM", "A": "0.5 mM"})
    with pytest.raises(SimulationSetupError, match="supports LRM and LRMP"):
        api.simulate(with_model(sma, "GRM", "reference"))
    equilibrium = make_binding_project(
        [{"name": "A"}],
        {"model": "Langmuir", "parameters": {"A": {"k_eq": 1.0, "q_max": 10.0}}},
        {"A": "1 mM"},
    )
    with pytest.raises(SimulationSetupError, match="give k_kin"):
        api.simulate(with_model(equilibrium, "LRM", "reference"))
    giex = make_binding_project(
        SALT + [{"name": "pH"}, {"name": "A"}],
        {
            "model": "GIEX",
            "parameters": {"A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0}},
        },
        {"salt": "50 mM", "pH": "5 pH", "A": "0.1 mM"},
    )
    with pytest.raises(SimulationSetupError, match="no Python flux"):
        api.simulate(with_model(giex, "LRM", "reference"))
