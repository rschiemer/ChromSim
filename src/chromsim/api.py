"""Public API of ChromSim.

This is the only module the CLI and the web UI import from (SPEC §4, §8). Functions here take and
return serialisable objects (dataclasses, pydantic models, file paths) so they can later be exposed
as a REST API without redesign (SPEC §8.1).
"""

from chromsim.analysis.detectors import predict_conductivity, predict_uv
from chromsim.analysis.fractionation import AssayComparison, compare_with_assay, fraction_averages
from chromsim.backends import SimulationResult, get_backend
from chromsim.backends.checks import DiscretisationCheck
from chromsim.backends.checks import check_discretisation as _check_discretisation
from chromsim.core.column import Column, ColumnModel
from chromsim.core.components import Component
from chromsim.core.method import ResolvedStep, Step
from chromsim.core.project import IsothermConfig, Project, SolverSettings
from chromsim.core.resin import Resin
from chromsim.core.system import System
from chromsim.data.assays import (
    AssayMappingProfile,
    FractionAssay,
    load_mapping_profile,
    read_assay,
)
from chromsim.data.experiment import Experiment, Fraction, Signal
from chromsim.data.importers import list_importers, read_chromatogram
from chromsim.data.preprocessing import (
    ConductivityCalibration,
    conductivity_to_salt,
    load_conductivity_calibration,
    resample,
    subtract_baseline,
    uv_to_concentration,
    volume_to_time,
)
from chromsim.diagnostics import DoctorReport, SmokeTestResult, run_doctor
from chromsim.errors import (
    ChromSimError,
    EnvironmentCheckError,
    IsothermConfigError,
    ProjectFileError,
    SimulationError,
    SimulationSetupError,
    UnitError,
)
from chromsim.io.project_file import load_project, project_from_yaml, project_to_yaml, save_project
from chromsim.isotherms.registry import list_isotherms
from chromsim.plotting import plot_chromatogram, plot_fraction_comparison, plot_overlay
from chromsim.units import Quantity

__all__ = [
    # errors
    "ChromSimError",
    "EnvironmentCheckError",
    "IsothermConfigError",
    "ProjectFileError",
    "SimulationError",
    "SimulationSetupError",
    "UnitError",
    # environment
    "DoctorReport",
    "SmokeTestResult",
    "doctor",
    # domain model
    "Column",
    "ColumnModel",
    "Component",
    "IsothermConfig",
    "Project",
    "Quantity",
    "ResolvedStep",
    "Resin",
    "SolverSettings",
    "Step",
    "System",
    # project files
    "load_project",
    "project_from_yaml",
    "project_to_yaml",
    "save_project",
    # isotherms
    "list_isotherms",
    # simulation
    "DiscretisationCheck",
    "SimulationResult",
    "check_discretisation",
    "plot_chromatogram",
    "simulate",
    # experimental data
    "Experiment",
    "Fraction",
    "Signal",
    "list_importers",
    "read_chromatogram",
    "AssayMappingProfile",
    "FractionAssay",
    "load_mapping_profile",
    "read_assay",
    "ConductivityCalibration",
    "load_conductivity_calibration",
    "conductivity_to_salt",
    "resample",
    "subtract_baseline",
    "uv_to_concentration",
    "volume_to_time",
    # comparing simulation and experiment
    "AssayComparison",
    "compare_with_assay",
    "fraction_averages",
    "plot_fraction_comparison",
    "plot_overlay",
    "predict_conductivity",
    "predict_uv",
]


def doctor() -> DoctorReport:
    """Check that ChromSim and CADET are installed and working.

    Collects the versions of Python, ChromSim, CADET-Process, CADET-Python and CADET-Core, and runs
    a sub-second CADET test simulation (tracer pulse through a non-binding column, checked for mass
    balance).

    Returns
    -------
    DoctorReport
        Versions, platform, smoke-test result and a list of problems; ``report.ok`` is ``True``
        if everything works.

    Examples
    --------
    >>> from chromsim import api
    >>> report = api.doctor()
    >>> report.ok
    True
    """
    return run_doctor()


def simulate(project: Project) -> SimulationResult:
    """Simulate a project with the backend named in ``project.solver.backend``.

    Parameters
    ----------
    project : Project
        A complete project (all transport parameters set for its column model).

    Returns
    -------
    SimulationResult
        Outlet concentrations over time and volume, with provenance.

    Raises
    ------
    SimulationSetupError
        If parameters are missing; the message lists all of them.
    SimulationError
        If the solver fails.
    """
    return get_backend(project.solver.backend).simulate(project)


def check_discretisation(project: Project) -> DiscretisationCheck:
    """Simulate on the project grid and on a doubled grid and report the change (SPEC §7.1).

    Parameters
    ----------
    project : Project
        The project to check (not modified).

    Returns
    -------
    DiscretisationCheck
        NRMSE per component between the two grids; well below 0.01 means the grid is fine enough.
    """
    return _check_discretisation(project)
