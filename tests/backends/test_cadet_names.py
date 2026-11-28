"""Pin every CADET-Process name used by ``chromsim/backends/cadet_backend.py`` and
``chromsim/isotherms/cadet_native.py`` (CLAUDE.md rule 4).

If a CADET-Process release renames something, these tests say exactly what changed.
"""

import tempfile
from pathlib import Path

import CADETProcess.processModel as process_model
import pytest
from CADETProcess.simulator import Cadet

from chromsim.isotherms.cadet_native import CadetNativeIsotherm
from chromsim.isotherms.registry import get_isotherm, list_isotherms

COMPONENTS = process_model.ComponentSystem(["salt", "protein"])


CADET_NATIVE = [n for n in list_isotherms() if isinstance(get_isotherm(n), CadetNativeIsotherm)]


@pytest.mark.parametrize("isotherm_name", CADET_NATIVE)
def test_isotherm_cadet_class_and_attributes_exist(isotherm_name):
    # Every registered CADET-native isotherm, using its own example configuration (so new
    # isotherms are pinned automatically). tests/plugins/ runs them in simulations too.
    isotherm = get_isotherm(isotherm_name)
    example = isotherm.example
    names = [component["name"] for component in example["components"]]
    binding = getattr(process_model, isotherm.cadet_class)(
        process_model.ComponentSystem(names), name="binding"
    )
    values = isotherm.cadet_parameters(
        names, example["parameters"], 1000.0, example.get("global_parameters")
    )
    for attribute in values:
        assert attribute in binding.parameters, f"{isotherm.cadet_class} has no `{attribute}`"
    for attribute in ("is_kinetic", "bound_states"):
        assert hasattr(binding, attribute), f"{isotherm.cadet_class} has no `{attribute}`"


def test_built_in_cadet_isotherms_are_all_pinned():
    expected = {"Linear", "Langmuir", "SMA", "GIEX", "MPM", "BiLangmuir"}
    expected |= {"HICConstantWaterActivity", "HICWaterOnHydrophobicSurfaces", "Colloidal", "SAS"}
    assert expected <= set(CADET_NATIVE)


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
