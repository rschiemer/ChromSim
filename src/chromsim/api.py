"""Public API of ChromSim.

This is the only module the CLI and the web UI import from (SPEC §4, §8). Functions here take and
return serialisable objects (dataclasses, pydantic models, file paths) so they can later be exposed
as a REST API without redesign (SPEC §8.1).
"""

from chromsim.core.column import Column, ColumnModel
from chromsim.core.components import Component
from chromsim.core.method import ResolvedStep, Step
from chromsim.core.project import IsothermConfig, Project, SolverSettings
from chromsim.core.resin import Resin
from chromsim.core.system import System
from chromsim.diagnostics import DoctorReport, SmokeTestResult, run_doctor
from chromsim.errors import ChromSimError, EnvironmentCheckError, ProjectFileError, UnitError
from chromsim.io.project_file import load_project, project_from_yaml, project_to_yaml, save_project
from chromsim.units import Quantity

__all__ = [
    # errors
    "ChromSimError",
    "EnvironmentCheckError",
    "ProjectFileError",
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
