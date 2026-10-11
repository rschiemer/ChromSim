"""Isotherm descriptors and registry (SPEC §5): parameter checks and the mapping to CADET values."""

import pytest
from project_examples import make_project
from pydantic import ValidationError

from chromsim.errors import ChromSimError
from chromsim.isotherms.base import Parameter, validate_parameters
from chromsim.isotherms.cadet_native import CadetNativeIsotherm, Langmuir, Linear, StericMassAction
from chromsim.isotherms.registry import get_isotherm, list_isotherms, register_isotherm

# --- registry -------------------------------------------------------------------------------------


def test_built_in_isotherms_are_registered():
    assert {"Linear", "Langmuir", "SMA"} <= set(list_isotherms())
    assert isinstance(get_isotherm("SMA"), StericMassAction)


def test_unknown_isotherm_lists_available_ones():
    with pytest.raises(
        ChromSimError, match="unknown isotherm 'SMAA'. Available isotherms: .*Langmuir.*SMA"
    ):
        get_isotherm("SMAA")


def test_duplicate_registration_is_rejected():
    class AnotherSMA(CadetNativeIsotherm):
        name = "SMA"
        cadet_class = "StericMassAction"
        parameters: dict[str, Parameter] = {}

        def cadet_parameters(self, component_names, parameters, ionic_capacity_mol_per_m3):
            return {}

    with pytest.raises(ChromSimError, match="already registered"):
        register_isotherm(AnotherSMA)


# --- parameter validation -------------------------------------------------------------------------


def test_unknown_parameter_is_rejected():
    with pytest.raises(ValueError, match="'mAb'.*unknown parameter 'kq'.*k_eq, k_kin"):
        validate_parameters(Linear(), ["salt", "mAb"], {"mAb": {"kq": 1.0}})


def test_missing_required_parameter_is_rejected():
    with pytest.raises(ValueError, match="'mAb'.*missing.*q_max"):
        validate_parameters(Langmuir(), ["mAb"], {"mAb": {"k_eq": 1.0}})


def test_parameter_bounds_are_checked():
    params = {"mAb": {"k_eq": 0.01, "nu": 60.0, "sigma": 30.0}}
    with pytest.raises(ValueError, match="nu = 60.*between 0 and 50"):
        validate_parameters(StericMassAction(), ["salt", "mAb"], params)


def test_k_kin_must_be_given_for_all_or_none():
    params = {"a": {"k_eq": 1.0, "k_kin": 10.0}, "b": {"k_eq": 2.0}}
    with pytest.raises(ValueError, match="k_kin.*all binding components or none"):
        validate_parameters(Linear(), ["a", "b"], params)


def test_sma_salt_must_be_component_zero():
    # Exact wording from SPEC §10.
    params = {"mAb": {"k_eq": 0.01, "nu": 6.0, "sigma": 30.0}}
    with pytest.raises(
        ValueError, match="Salt must be component 0 for SMA; found 'mAb' at index 0"
    ):
        validate_parameters(StericMassAction(), ["mAb", "salt"], params)


def test_sma_salt_cannot_have_binding_parameters():
    params = {"salt": {"k_eq": 0.01, "nu": 6.0, "sigma": 30.0}}
    with pytest.raises(ValueError, match="component 0 \\('salt'\\) is the salt"):
        validate_parameters(StericMassAction(), ["salt", "mAb"], params)


# --- project-level checks -------------------------------------------------------------------------


def test_project_rejects_unknown_isotherm():
    with pytest.raises(ValidationError, match="unknown isotherm 'Freundlich'"):
        make_project(isotherm={"model": "Freundlich", "parameters": {}})


def test_project_sma_needs_ionic_capacity():
    resin = {"name": "r", "particle_radius": "45 um", "col_porosity": 0.36, "par_porosity": 0.85}
    with pytest.raises(ValidationError, match="SMA needs resin.ionic_capacity"):
        make_project(resin=resin)


def test_project_sma_salt_position():
    components = [{"name": "mAb", "mw_kDa": 148}, {"name": "salt"}]
    with pytest.raises(ValidationError, match="Salt must be component 0 for SMA"):
        make_project(components=components)


# --- mapping to CADET-Process values --------------------------------------------------------------


def test_linear_rapid_equilibrium_mapping():
    isotherm = Linear()
    params = {"A": {"k_eq": 2.0}}
    assert not isotherm.is_kinetic(params)
    assert isotherm.bound_states(["salt", "A"], params) == [0, 1]
    assert isotherm.cadet_parameters(["salt", "A"], params, None) == {
        "adsorption_rate": [0.0, 2.0],  # k_a = k_eq when rapid equilibrium (k_d = 1)
        "desorption_rate": [0.0, 1.0],
    }


def test_linear_kinetic_mapping():
    isotherm = Linear()
    params = {"A": {"k_eq": 2.0, "k_kin": 5.0}}
    assert isotherm.is_kinetic(params)
    values = isotherm.cadet_parameters(["A"], params, None)
    assert values == {"adsorption_rate": [10.0], "desorption_rate": [5.0]}  # k_a = k_eq k_kin


def test_langmuir_mapping():
    params = {"A": {"k_eq": 0.5, "q_max": 100.0}}
    values = Langmuir().cadet_parameters(["salt", "A"], params, None)
    assert values == {
        "adsorption_rate": [0.0, 0.5],
        "desorption_rate": [0.0, 1.0],
        "capacity": [0.0, 100.0],
    }


def test_sma_mapping_uses_ionic_capacity_as_reference_concentrations():
    params = {"mAb": {"k_eq": 1e-2, "nu": 6.0, "sigma": 30.0, "k_kin": 1e2}}
    isotherm = StericMassAction()
    values = isotherm.cadet_parameters(["salt", "mAb"], params, 1200.0)
    assert values == {
        "adsorption_rate": [0.0, 1.0],  # k_a = k_eq k_kin = 1e-2 × 1e2
        "desorption_rate": [0.0, 100.0],
        "characteristic_charge": [0.0, 6.0],
        "steric_factor": [0.0, 30.0],
        "capacity": 1200.0,
        "reference_liquid_phase_conc": 1200.0,
        "reference_solid_phase_conc": 1200.0,
    }
    assert isotherm.bound_states(["salt", "mAb", "hcp"], params) == [1, 1, 0]


def test_sma_initial_bound_state_is_all_salt():
    # Column starts equilibrated: every site carries a salt counter-ion, q_salt = Λ.
    q = StericMassAction().initial_bound_state(["salt", "mAb"], {"mAb": {}}, 1200.0)
    assert q == [1200.0, 0.0]


def test_modifier_components():
    assert StericMassAction().modifier_components(["salt", "mAb"]) == ["salt"]
    assert Langmuir().modifier_components(["a", "b"]) == []
