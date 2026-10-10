"""Pin every CADET-Process name used by ``chromsim/backends/cadet_backend.py`` and
``chromsim/isotherms/cadet_native.py`` (CLAUDE.md rule 4).

If a CADET-Process release renames something, these tests say exactly what changed.
"""

import tempfile
from pathlib import Path

import CADETProcess.processModel as process_model
import pytest
from CADETProcess.simulator import Cadet

from chromsim.isotherms.registry import get_isotherm, list_isotherms

COMPONENTS = process_model.ComponentSystem(["salt", "protein"])


@pytest.mark.parametrize("isotherm_name", ["Linear", "Langmuir", "SMA"])
def test_isotherm_cadet_class_and_attributes_exist(isotherm_name):
    isotherm = get_isotherm(isotherm_name)
    binding = getattr(process_model, isotherm.cadet_class)(COMPONENTS, name="binding")
    example_values = {"k_eq": 1.0, "nu": 5.0, "sigma": 10.0, "q_max": 10.0}
    parameters = {"protein": {k: v for k, v in example_values.items() if k in isotherm.parameters}}
    values = isotherm.cadet_parameters(["salt", "protein"], parameters, 1000.0)
    for attribute in values:
        assert attribute in binding.parameters, f"{isotherm.cadet_class} has no `{attribute}`"
    for attribute in ("is_kinetic", "bound_states"):
        assert hasattr(binding, attribute), f"{isotherm.cadet_class} has no `{attribute}`"


def test_all_registered_cadet_isotherms_are_checked():
    # When a new CADET-native isotherm is registered, add it to the parametrisation above.
    assert set(list_isotherms()) == {"Linear", "Langmuir", "SMA"}


def test_no_binding_class_exists():
    assert process_model.NoBinding(COMPONENTS, name="binding") is not None


@pytest.mark.parametrize(
    ("class_name", "attributes"),
    [
        (
            "LumpedRateModelWithoutPores",
            ["length", "diameter", "total_porosity", "axial_dispersion", "c", "q"],
        ),
        (
            "LumpedRateModelWithPores",
            [
                "length",
                "diameter",
                "bed_porosity",
                "particle_porosity",
                "particle_radius",
                "axial_dispersion",
                "film_diffusion",
                "c",
                "cp",
                "q",
            ],
        ),
        (
            "GeneralRateModel",
            [
                "length",
                "diameter",
                "bed_porosity",
                "particle_porosity",
                "particle_radius",
                "axial_dispersion",
                "film_diffusion",
                "pore_diffusion",
                "c",
                "cp",
                "q",
            ],
        ),
        ("TubularReactor", ["length", "diameter", "axial_dispersion", "c"]),
        ("Cstr", ["init_liquid_volume", "c"]),
        ("Inlet", ["c", "flow_rate"]),
    ],
)
def test_unit_operation_parameters(class_name, attributes):
    unit = getattr(process_model, class_name)(COMPONENTS, name="unit")
    for attribute in attributes:
        assert attribute in unit.parameters, f"{class_name} has no parameter `{attribute}`"


@pytest.mark.parametrize(
    ("class_name", "attributes"),
    [
        ("LumpedRateModelWithoutPores", ["ncol"]),
        ("LumpedRateModelWithPores", ["ncol"]),
        ("GeneralRateModel", ["ncol", "npar"]),
        ("TubularReactor", ["ncol"]),
    ],
)
def test_discretization_parameters(class_name, attributes):
    unit = getattr(process_model, class_name)(COMPONENTS, name="unit")
    for attribute in attributes:
        assert attribute in unit.discretization.parameters


def test_unit_binding_model_attribute():
    unit = process_model.GeneralRateModel(COMPONENTS, name="unit")
    assert hasattr(unit, "binding_model")


def test_process_events_take_polynomial_inlet_concentrations():
    flow_sheet = process_model.FlowSheet(COMPONENTS)
    inlet = process_model.Inlet(COMPONENTS, name="inlet")
    outlet = process_model.Outlet(COMPONENTS, name="outlet")
    flow_sheet.add_unit(inlet)
    flow_sheet.add_unit(outlet)
    flow_sheet.add_connection(inlet, outlet)
    process = process_model.Process(flow_sheet, "pinning")
    process.cycle_time = 10.0
    process.add_event("c", "flow_sheet.inlet.c", [[1.0, 0.5], [0.0, 0.0]], 0.0)
    process.add_event("flow", "flow_sheet.inlet.flow_rate", 1e-8, 0.0)
    assert [event.name for event in process.events] == ["c", "flow"]
    assert [unit.name for unit in flow_sheet.units] == ["inlet", "outlet"]


def test_simulator_time_integrator_abstol():
    with tempfile.TemporaryDirectory() as temp_dir:
        simulator = Cadet(temp_dir=Path(temp_dir))
        assert "abstol" in simulator.time_integrator_parameters.parameters
