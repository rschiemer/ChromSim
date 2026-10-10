"""Real CADET checks behind ``chromsim doctor``.

The first group pins every CADET-Process name used in ``chromsim/diagnostics.py`` (CLAUDE.md
rule 4): if a CADET-Process release renames one, these tests fail with a clear message instead of
``doctor`` breaking in an obscure way.
"""

import shutil
import subprocess

from CADETProcess.processModel import (
    ComponentSystem,
    FlowSheet,
    Inlet,
    LumpedRateModelWithoutPores,
    Outlet,
    Process,
)
from CADETProcess.simulator import Cadet

from chromsim import api, diagnostics

# Wall-time limit for the smoke simulation. It takes ~0.1 s on a laptop; 5 s leaves room for slow
# CI runners and first-call overhead while still catching a solver that is badly misconfigured.
SMOKE_TEST_MAX_WALL_TIME_S = 5.0


def test_cadet_process_inlet_parameter_names():
    inlet = Inlet(ComponentSystem(1), name="inlet")
    assert {"flow_rate", "c"} <= set(inlet.parameters)


def test_cadet_process_lrm_parameter_names():
    column = LumpedRateModelWithoutPores(ComponentSystem(1), name="column")
    assert {"length", "diameter", "total_porosity", "axial_dispersion"} <= set(
        column.required_parameters
    )
    assert "ncol" in column.discretization.parameters


def test_cadet_process_flow_sheet_and_process_api():
    component_system = ComponentSystem(1)
    flow_sheet = FlowSheet(component_system)
    inlet = Inlet(component_system, name="inlet")
    outlet = Outlet(component_system, name="outlet")
    flow_sheet.add_unit(inlet)
    flow_sheet.add_unit(outlet)
    flow_sheet.add_connection(inlet, outlet)
    process = Process(flow_sheet, "pinning")
    process.cycle_time = 10.0
    process.add_event("start", "flow_sheet.inlet.c", [1.0], 0.0)
    assert process.cycle_time == 10.0


def test_cadet_simulator_api():
    # `SimulationResults.exit_flag`, `.exit_message` and `.solution.<unit>.outlet` are exercised by
    # the real simulation in `test_smoke_simulation_conserves_mass_quickly`.
    for attribute in ("check_cadet", "simulate", "version"):
        assert hasattr(Cadet, attribute), f"CADET-Process `Cadet` has no `{attribute}`"


def test_cadet_core_is_major_version_5():
    # SPEC §3: CADET-Core 5.x; moving to 6.x requires an ADR.
    assert diagnostics.find_cadet_core_version().split(".")[0] == "5"


def test_smoke_simulation_conserves_mass_quickly():
    result = diagnostics.run_cadet_smoke_test()
    assert abs(result.mass_recovery - 1.0) <= diagnostics.MASS_RECOVERY_TOLERANCE
    assert result.wall_time_s < SMOKE_TEST_MAX_WALL_TIME_S


def test_api_doctor_is_healthy():
    report = api.doctor()
    assert report.ok, report.problems
    for name in ("Python", "ChromSim", "CADET-Process", "CADET-Python", "CADET-Core"):
        assert report.versions[name] not in ("not installed", "not found")


def test_chromsim_doctor_command_end_to_end():
    # Runs the installed console entry point exactly as a user would.
    executable = shutil.which("chromsim")
    assert executable is not None, "`chromsim` command not on PATH; run `pixi install`"
    completed = subprocess.run([executable, "doctor"], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "CADET-Core" in completed.stdout
    assert "CADET test simulation: OK" in completed.stdout
