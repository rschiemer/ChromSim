"""Environment diagnostics behind ``chromsim doctor``.

``run_doctor`` collects the versions of Python, ChromSim, CADET-Process, CADET-Python and
CADET-Core, then runs a tiny CADET simulation to prove the solver works end to end.

The smoke test simulates a rectangular tracer pulse through a lumped rate model column without
binding and checks that all the injected mass leaves the column (mass balance). It takes well
under a second of wall time.

All CADET-Process names used here (``Inlet.flow_rate``, ``LumpedRateModelWithoutPores.length``,
``Cadet.version`` …) were read from the installed CADET-Process 0.12 and are pinned by
``tests/backends/test_doctor_cadet.py`` (CLAUDE.md rule 4).
"""

import contextlib
import importlib.metadata
import io
import logging
import platform
import sys
import time
from dataclasses import dataclass, field

import numpy as np

import chromsim
from chromsim.errors import EnvironmentCheckError

logger = logging.getLogger(__name__)

# Smoke-test tolerance on mass recovery. CADET recovers the pulse to ~1e-8 relative error on
# this case; 1 % is loose enough for any solver/platform yet catches a broken installation.
MASS_RECOVERY_TOLERANCE = 0.01


@dataclass
class SmokeTestResult:
    """Outcome of the CADET smoke-test simulation.

    Attributes
    ----------
    wall_time_s : float
        Wall-clock time of the CADET simulation call, in s.
    mass_recovery : float
        Eluted mass divided by injected mass, dimensionless (1.0 = perfect mass balance).
    """

    wall_time_s: float
    mass_recovery: float


@dataclass
class DoctorReport:
    """Result of ``chromsim doctor``.

    Attributes
    ----------
    versions : dict[str, str]
        Package name → version string. Missing packages are reported as ``"not installed"``.
    platform : str
        Operating system and machine architecture.
    smoke_test : SmokeTestResult or None
        Result of the CADET test simulation, ``None`` if it could not be run.
    problems : list[str]
        Actionable descriptions of everything that went wrong. Empty if the environment is healthy.
    """

    versions: dict[str, str]
    platform: str
    smoke_test: SmokeTestResult | None = None
    problems: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """``True`` if all checks passed."""
        return self.smoke_test is not None and not self.problems


def _package_version(distribution_name: str) -> str:
    """Return the installed version of a distribution, or ``"not installed"``."""
    try:
        return importlib.metadata.version(distribution_name)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def collect_versions() -> dict[str, str]:
    """Collect the versions of Python and the Python packages ChromSim depends on.

    CADET-Core is not a Python package; its version is added by :func:`run_doctor` after
    the CADET binary has been located.

    Returns
    -------
    dict[str, str]
        Name → version string.
    """
    return {
        "Python": platform.python_version(),
        "ChromSim": chromsim.__version__,
        "CADET-Process": _package_version("CADET-Process"),
        "CADET-Python": _package_version("CADET-Python"),
    }


def find_cadet_core_version() -> str:
    """Locate the CADET-Core binary through CADET-Process and return its version.

    Returns
    -------
    str
        CADET-Core version, e.g. ``"5.1.1"``.

    Raises
    ------
    EnvironmentCheckError
        If CADET-Process is missing or cannot find a working CADET-Core installation.
    """
    try:
        from CADETProcess.simulator import Cadet
    except ImportError as error:
        raise EnvironmentCheckError(
            f"CADET-Process cannot be imported ({error}). Run `pixi install` in the ChromSim "
            "folder to set up the environment."
        ) from error

    try:
        simulator = Cadet()
        # check_cadet() prints its own success message; keep it out of the doctor report.
        captured_output = io.StringIO()
        with contextlib.redirect_stdout(captured_output):
            cadet_works = simulator.check_cadet()
        logger.debug("CADET self-check output: %s", captured_output.getvalue().strip())
    except Exception as error:  # CADET-Process raises plain exceptions for a missing binary
        raise EnvironmentCheckError(
            f"CADET-Core could not be found or started ({error}). Install it with `pixi install` "
            "(it comes from the conda-forge package `cadet`)."
        ) from error
    if not cadet_works:
        raise EnvironmentCheckError(
            "CADET-Core was found but failed its own self-check. Reinstall the environment with "
            "`pixi clean && pixi install`."
        )
    return str(simulator.version)


def run_cadet_smoke_test() -> SmokeTestResult:
    """Simulate a tracer pulse through a non-binding LRM column with CADET.

    A rectangular pulse of feed concentration c_feed = 1 mol/m³ and duration t_pulse = 60 s is
    injected into a 10 cm × 1 cm column (total porosity 0.6, axial dispersion 1e-7 m²/s) at
    Q = 1e-8 m³/s. The column residence time is about 470 s; the simulation runs for 1200 s, so
    the whole pulse has left the column. Mass recovery is

    .. math:: R = \\frac{\\int_0^{t_{end}} c_{out}(t)\\,dt}{c_{feed}\\, t_{pulse}}

    (the flow rate cancels because it is constant).

    Returns
    -------
    SmokeTestResult
        Wall time of the simulation and the mass recovery.

    Raises
    ------
    EnvironmentCheckError
        If CADET-Process is missing, the simulation fails, or mass recovery is outside
        ``1 ± MASS_RECOVERY_TOLERANCE``.
    """
    try:
        from CADETProcess.processModel import (
            ComponentSystem,
            FlowSheet,
            Inlet,
            LumpedRateModelWithoutPores,
            Outlet,
            Process,
        )
        from CADETProcess.simulator import Cadet
    except ImportError as error:
        raise EnvironmentCheckError(
            f"CADET-Process cannot be imported ({error}). Run `pixi install` in the ChromSim "
            "folder to set up the environment."
        ) from error

    feed_concentration = 1.0  # mol/m³
    pulse_duration = 60.0  # s

    component_system = ComponentSystem(["tracer"])

    inlet = Inlet(component_system, name="inlet")
    inlet.flow_rate = 1.0e-8  # m³/s

    column = LumpedRateModelWithoutPores(component_system, name="column")
    column.length = 0.10  # m
    column.diameter = 0.01  # m
    column.total_porosity = 0.6  # -
    column.axial_dispersion = 1.0e-7  # m²/s
    column.discretization.ncol = 50

    outlet = Outlet(component_system, name="outlet")

    flow_sheet = FlowSheet(component_system)
    flow_sheet.add_unit(inlet)
    flow_sheet.add_unit(column)
    flow_sheet.add_unit(outlet)
    flow_sheet.add_connection(inlet, column)
    flow_sheet.add_connection(column, outlet)

    process = Process(flow_sheet, "doctor_tracer_pulse")
    process.cycle_time = 1200.0  # s
    process.add_event("pulse_start", "flow_sheet.inlet.c", [feed_concentration], 0.0)
    process.add_event("pulse_end", "flow_sheet.inlet.c", [0.0], pulse_duration)

    start = time.perf_counter()
    try:
        results = Cadet().simulate(process)
    except Exception as error:  # CADET-Process raises plain exceptions on solver failure
        raise EnvironmentCheckError(
            f"The CADET test simulation failed ({error}). Check the CADET-Core installation with "
            "`pixi clean && pixi install`."
        ) from error
    wall_time_s = time.perf_counter() - start

    if results.exit_flag != 0:
        raise EnvironmentCheckError(
            f"The CADET test simulation ended with exit flag {results.exit_flag}: "
            f"{results.exit_message}. Reinstall the environment with `pixi clean && pixi install`."
        )

    outlet_solution = results.solution.outlet.outlet
    time_s = outlet_solution.time
    c_out = outlet_solution.solution[:, 0]
    mass_recovery = float(np.trapezoid(c_out, time_s) / (feed_concentration * pulse_duration))

    if abs(mass_recovery - 1.0) > MASS_RECOVERY_TOLERANCE:
        raise EnvironmentCheckError(
            f"The CADET test simulation recovered {mass_recovery:.2%} of the injected mass "
            f"(expected 100 ± {MASS_RECOVERY_TOLERANCE:.0%}). The solver is not working correctly; "
            "reinstall the environment with `pixi clean && pixi install`."
        )

    logger.info(
        "CADET smoke test: %.3f s wall time, mass recovery %.6f", wall_time_s, mass_recovery
    )
    return SmokeTestResult(wall_time_s=wall_time_s, mass_recovery=mass_recovery)


def run_doctor() -> DoctorReport:
    """Check the ChromSim environment.

    Never raises for environment problems; they are collected in ``DoctorReport.problems`` so the
    versions can always be shown.

    Returns
    -------
    DoctorReport
        Versions, platform, smoke-test result and any problems found.
    """
    report = DoctorReport(
        versions=collect_versions(),
        platform=f"{platform.system()} {platform.machine()} ({sys.platform})",
    )
    try:
        report.versions["CADET-Core"] = find_cadet_core_version()
    except EnvironmentCheckError as error:
        report.versions["CADET-Core"] = "not found"
        report.problems.append(str(error))
        return report

    try:
        report.smoke_test = run_cadet_smoke_test()
    except EnvironmentCheckError as error:
        report.problems.append(str(error))
    return report
