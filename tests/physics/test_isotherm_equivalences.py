"""Each additional isotherm reduces to a simpler one in a special case (exact identities).

Same column, grid and solver for both sides, so outlet profiles must agree to solver precision.
This checks the mapping of ChromSim parameters onto CADET-Process names physically.
"""

import numpy as np
import pytest
from project_examples import make_binding_project

from chromsim import api
from chromsim.backends.checks import nrmse

# Both simulations solve the same equations on the same grid; only round-off and adaptive time
# stepping differ (CADET reltol 1e-6). 1e-4 is far below any mapping error (factors of order 1).
IDENTITY_NRMSE = 1e-4

SALT = [{"name": "salt"}]
AB = [{"name": "A"}, {"name": "B"}]


def outlet(project, component):
    return api.simulate(project).outlet(component)


def test_mpm_without_salt_dependence_is_langmuir():
    # β = 0, γ = 0: e^{γ c_s} = 1 and c_s^β = 1.
    load = {"salt": "50 mM", "A": "1 mM"}
    mpm = make_binding_project(
        SALT + [{"name": "A"}],
        {
            "model": "MPM",
            "parameters": {
                "A": {"k_eq": 2.0, "q_max": 10.0, "beta": 0.0, "gamma": 0.0, "k_kin": 1.0}
            },
        },
        load,
    )
    langmuir = make_binding_project(
        SALT + [{"name": "A"}],
        {"model": "Langmuir", "parameters": {"A": {"k_eq": 2.0, "q_max": 10.0, "k_kin": 1.0}}},
        load,
    )
    assert nrmse(outlet(mpm, "A"), outlet(langmuir, "A")) < IDENTITY_NRMSE


@pytest.mark.parametrize(("nu_lin", "k_eq_lin"), [(0.0, 0.0), (0.4, 0.0), (0.0, 0.3)])
def test_giex_at_constant_ph_is_sma(nu_lin, k_eq_lin):
    # At constant pH p: ν_eff = ν + ν_lin p and k_eq,eff = k_eq exp(k_eq_lin p).
    ph = 5.0
    components = SALT + [{"name": "pH"}, {"name": "A"}]
    giex_parameters = {
        "k_eq": 1e-2,
        "nu": 4.0,
        "nu_lin": nu_lin,
        "sigma": 20.0,
        "k_eq_lin": k_eq_lin,
        "k_kin": 10.0,
    }
    giex = make_binding_project(
        components,
        {"model": "GIEX", "parameters": {"A": giex_parameters}},
        {"salt": "50 mM", "pH": f"{ph} pH", "A": "0.1 mM"},
        elute={"salt": "1000 mM", "pH": f"{ph} pH"},
    )
    sma_parameters = {
        "k_eq": 1e-2 * np.exp(k_eq_lin * ph),
        "nu": 4.0 + nu_lin * ph,
        "sigma": 20.0,
        "k_kin": 10.0,
    }
    sma = make_binding_project(
        components,
        {"model": "SMA", "parameters": {"A": sma_parameters}},
        {"salt": "50 mM", "pH": f"{ph} pH", "A": "0.1 mM"},
        elute={"salt": "1000 mM", "pH": f"{ph} pH"},
    )
    assert nrmse(outlet(giex, "A"), outlet(sma, "A")) < IDENTITY_NRMSE


def test_bi_langmuir_with_identical_sites_is_langmuir_with_double_capacity():
    # Two identical, independent site types each hold half:
    # dq/dt = k_a c (2 q_max)(1 − q/(2 q_max)) − k_d q.
    # Two components with different parameters also check the state-major ordering of the arrays.
    load = {"A": "1 mM", "B": "0.5 mM"}
    bi = make_binding_project(
        AB,
        {
            "model": "BiLangmuir",
            "parameters": {
                "A": {
                    "k_eq_1": 2.0,
                    "q_max_1": 5.0,
                    "k_eq_2": 2.0,
                    "q_max_2": 5.0,
                    "k_kin_1": 1.0,
                    "k_kin_2": 1.0,
                },
                "B": {
                    "k_eq_1": 0.5,
                    "q_max_1": 3.0,
                    "k_eq_2": 0.5,
                    "q_max_2": 3.0,
                    "k_kin_1": 2.0,
                    "k_kin_2": 2.0,
                },
            },
        },
        load,
    )
    langmuir = make_binding_project(
        AB,
        {
            "model": "Langmuir",
            "parameters": {
                "A": {"k_eq": 2.0, "q_max": 10.0, "k_kin": 1.0},
                "B": {"k_eq": 0.5, "q_max": 6.0, "k_kin": 2.0},
            },
        },
        load,
    )
    bi_result, langmuir_result = api.simulate(bi), api.simulate(langmuir)
    for component in ("A", "B"):
        assert (
            nrmse(bi_result.outlet(component), langmuir_result.outlet(component)) < IDENTITY_NRMSE
        )


@pytest.mark.parametrize("model", ["HICConstantWaterActivity", "HICWaterOnHydrophobicSurfaces"])
def test_hic_with_one_ligand_and_no_water_term_is_langmuir(model):
    # ν = 1, β₀ = 0 (β = 0): dq/dt = k_a c (1 − q/q_max) − k_d q, i.e. Langmuir with
    # k_eq,Langmuir = k_eq,HIC / q_max (Langmuir multiplies k_a by q_max).
    load = {"salt": "1000 mM", "A": "1 mM"}
    hic = make_binding_project(
        SALT + [{"name": "A"}],
        {
            "model": model,
            "parameters": {"A": {"k_eq": 20.0, "q_max": 10.0, "nu": 1.0, "k_kin": 1.0}},
            "global_parameters": {"beta_0": 0.0, "beta_1": 0.0},
        },
        load,
    )
    langmuir = make_binding_project(
        SALT + [{"name": "A"}],
        {"model": "Langmuir", "parameters": {"A": {"k_eq": 2.0, "q_max": 10.0, "k_kin": 1.0}}},
        load,
    )
    assert nrmse(outlet(hic, "A"), outlet(langmuir, "A")) < IDENTITY_NRMSE


def test_colloidal_without_lattice_and_salt_dependence_is_linear():
    # b_pp = 0 removes the lattice term and k_b = 0 makes ln K_e = k_a constant:
    # dq/dt = k_kin (c − q/K_e) is Linear with k_eq = K_e and k_d = k_kin/K_e.
    load = {"salt": "50 mM", "A": "1 mM"}
    k_e, k_kin = 2.0, 1.0
    colloidal = make_binding_project(
        SALT + [{"name": "A"}],
        {
            "model": "Colloidal",
            "parameters": {
                "A": {
                    "radius": 4e-9,
                    "ln_ke_salt_factor": float(np.log(k_e)),
                    "ln_ke_salt_exponent": 0.0,
                    "b_pp_salt_factor": 0.0,
                    "b_pp_salt_exponent": 0.0,
                    "k_kin": k_kin,
                }
            },
            "global_parameters": {
                "phase_ratio": 5e7,
                "kappa_factor": 0.0,
                "kappa_exponent": 0.0,
                "kappa_constant": 1e-9,
            },
        },
        load,
    )
    linear = make_binding_project(
        SALT + [{"name": "A"}],
        {"model": "Linear", "parameters": {"A": {"k_eq": k_e, "k_kin": k_kin / k_e}}},
        load,
    )
    assert nrmse(outlet(colloidal, "A"), outlet(linear, "A")) < IDENTITY_NRMSE


def test_sas_without_self_association_is_sma():
    # k_sa = 0 removes the dimerisation term: (k_a1 + k_a2 c) c → k_a1 c.
    load = {"salt": "50 mM", "A": "0.5 mM"}
    sma = {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0}
    projects = [
        make_binding_project(
            SALT + [{"name": "A"}],
            {"model": model, "parameters": {"A": values}},
            load,
            elute={"salt": "1000 mM"},
        )
        for model, values in (("SMA", sma), ("SAS", sma | {"self_association": 0.0}))
    ]
    assert nrmse(outlet(projects[1], "A"), outlet(projects[0], "A")) < IDENTITY_NRMSE
