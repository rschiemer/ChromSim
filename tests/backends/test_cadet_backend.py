"""Behaviour of the CADET backend: flow-sheet assembly, errors, result axes and provenance."""

import tempfile

import numpy as np
import pytest
from project_examples import PHYSICS_COLUMN_VOLUME_ML, make_physics_project, make_project

from chromsim.backends import SimulationResult, get_backend
from chromsim.backends.cadet_backend import build_process, check_simulatable
from chromsim.backends.checks import check_discretisation, nrmse
from chromsim.errors import SimulationSetupError
from chromsim.io.project_file import project_from_yaml


def example_with_transport():
    project = make_project()
    project.column.axial_dispersion = "2e-7 m^2/s"
    project.column.film_diffusion = {"salt": "1e-4 m/s", "mAb": "2e-6 m/s"}
    project.column.pore_diffusion = {"salt": "1e-10 m^2/s", "mAb": "5e-12 m^2/s"}
    project.solver.n_col, project.solver.n_par = 20, 5
    return project


# --- setup checks ---------------------------------------------------------------------------------


def test_missing_transport_parameters_are_all_listed():
    with pytest.raises(SimulationSetupError) as excinfo:
        check_simulatable(make_project())
    message = str(excinfo.value)
    assert "column.axial_dispersion is not set" in message
    assert "column.film_diffusion is missing for salt, mAb" in message
    assert "column.pore_diffusion is missing for salt, mAb" in message


def test_lrm_needs_only_axial_dispersion():
    project = make_physics_project("LRM")
    project.column.film_diffusion = None
    project.column.pore_diffusion = None
    check_simulatable(project)  # no error


def test_unknown_backend():
    assert get_backend("reference").name == "reference"
    with pytest.raises(SimulationSetupError, match="Unknown backend 'cadett'.*cadet, reference"):
        get_backend("cadett")


# --- flow sheet -----------------------------------------------------------------------------------


def test_flow_sheet_without_system_volumes():
    process = build_process(make_physics_project("LRM"))
    assert [unit.name for unit in process.flow_sheet.units] == ["inlet", "column", "outlet"]


def test_flow_sheet_with_system_volumes():
    project = example_with_transport()
    process = build_process(project)
    assert [unit.name for unit in process.flow_sheet.units] == [
        "inlet",
        "mixer",
        "tubing_in",
        "column",
        "tubing_out",
        "outlet",
    ]
    tubing = process.flow_sheet.tubing_in
    # 0.125 mL in 0.5 mm tubing: L = 1.25e-7 m³ / (π × 0.25e-3² m²) = 0.6366 m
    assert tubing.length == pytest.approx(0.63662, rel=1e-4)
    assert process.flow_sheet.mixer.init_liquid_volume == pytest.approx(0.2e-6)


def test_events_encode_gradient_slope():
    process = build_process(example_with_transport())
    elute = next(event for event in process.events if event.name == "step3_elute_c")
    # 50 → 500 mM over 20 CV = 9600 s at 150 cm/h: slope 450/9600 mol/m³/s for salt
    assert elute.state[0][0] == pytest.approx(50.0)
    assert elute.state[0][1] == pytest.approx(450.0 / 9600.0)
    assert elute.state[1] == pytest.approx([0.0, 0.0])
    assert len(process.events) == 2 * 4  # concentration + flow rate per step
    assert process.cycle_time == pytest.approx(project_end_s(example_with_transport()))


def project_end_s(project):
    return project.resolved_method()[-1].end_s


def test_initial_state_is_clean_column_in_starting_salt():
    column = build_process(example_with_transport()).flow_sheet.column
    assert column.c == pytest.approx([50.0, 0.0])
    assert column.cp == pytest.approx([50.0, 0.0])
    assert column.q == pytest.approx([1200.0, 0.0])  # SMA: all sites carry salt


def test_sma_maps_through_to_cadet():
    binding = build_process(example_with_transport()).flow_sheet.column.binding_model
    assert type(binding).__name__ == "StericMassAction"
    assert binding.is_kinetic
    assert binding.capacity == pytest.approx(1200.0)
    assert binding.characteristic_charge == pytest.approx([0.0, 6.0])


# --- results --------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tracer_result() -> SimulationResult:
    return get_backend("cadet").simulate(make_physics_project("LRM"))


def test_result_axes(tracer_result):
    result = tracer_result
    assert result.component_names == ["tracer"]
    assert result.time_s[0] == 0.0
    assert result.time_s[-1] == pytest.approx(6.0 + 1800.0)
    assert result.outlet_mol_per_m3.shape == (len(result.time_s), 1)
    # 1 mL/min for 1806 s = 30.1 mL = 3.8325 CV of 7.854 mL
    assert result.volume_m3[-1] == pytest.approx(30.1e-6)
    assert result.column_volumes[-1] == pytest.approx(30.1 / PHYSICS_COLUMN_VOLUME_ML)
    assert np.all(np.diff(result.volume_m3) >= 0)


def test_result_outlet_units():
    project = example_with_transport()
    result = get_backend("cadet").simulate(project)
    mab_mM = result.outlet("mAb")
    # g/L = mol/m³ × kg/mol; mw 148 kDa → 148 kg/mol
    assert result.outlet("mAb", "g/L") == pytest.approx(mab_mM * 148.0)
    assert result.outlet("salt", "M") == pytest.approx(result.outlet("salt") / 1000.0)


def test_result_provenance(tracer_result):
    provenance = tracer_result.provenance
    assert provenance["backend"] == "cadet"
    assert provenance["versions"]["CADET-Core"].startswith("5.")
    assert {"ChromSim", "CADET-Process", "Python"} <= set(provenance["versions"])
    assert provenance["created_by"]
    assert provenance["solver"]["n_col"] == 100
    # The stored project file reproduces the simulated project.
    assert project_from_yaml(provenance["project_yaml"]) == make_physics_project("LRM")


def test_simulation_leaves_no_files_in_working_directory(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    tempdir_before = tempfile.gettempdir()
    get_backend("cadet").simulate(make_physics_project("LRM", n_col=20))
    assert list(tmp_path.iterdir()) == []
    assert tempfile.gettempdir() == tempdir_before


# --- discretisation check -------------------------------------------------------------------------


def test_nrmse():
    reference = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
    assert nrmse(reference, reference) == 0.0
    assert nrmse(reference + 0.2, reference) == pytest.approx(0.1)  # 0.2 / range 2


def test_check_discretisation_doubles_grid_and_converges():
    check = check_discretisation(make_physics_project("GRM", n_col=50, n_par=5))
    assert (check.refined_n_col, check.refined_n_par) == (100, 10)
    # A smooth tracer peak on 50 cells is already within 1 % of the doubled grid.
    assert 0.0 < check.max_nrmse < 0.01
