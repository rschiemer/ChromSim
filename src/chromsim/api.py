"""Public API of ChromSim.

This is the only module the CLI and the web UI import from (SPEC §4, §8). Functions here take and
return serialisable objects (dataclasses, pydantic models, file paths) so they can later be exposed
as a REST API without redesign (SPEC §8.1).
"""

from chromsim.diagnostics import DoctorReport, SmokeTestResult, run_doctor
from chromsim.errors import ChromSimError, EnvironmentCheckError

__all__ = [
    "ChromSimError",
    "DoctorReport",
    "EnvironmentCheckError",
    "SmokeTestResult",
    "doctor",
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
