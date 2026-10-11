"""Simulator backends behind the ``SimulatorBackend`` protocol: CADET and reference.

See SPEC §4. Use :func:`get_backend` to pick one by the name in ``project.solver.backend``.
"""

from chromsim.backends.base import SimulationResult, SimulatorBackend
from chromsim.backends.cadet_backend import CadetBackend
from chromsim.backends.reference import ReferenceBackend
from chromsim.errors import SimulationSetupError

__all__ = [
    "CadetBackend",
    "ReferenceBackend",
    "SimulationResult",
    "SimulatorBackend",
    "get_backend",
]


def get_backend(name: str) -> SimulatorBackend:
    """Return the simulation backend called ``name``.

    Raises
    ------
    SimulationSetupError
        For unknown names.
    """
    if name == "cadet":
        return CadetBackend()
    if name == "reference":
        return ReferenceBackend()
    raise SimulationSetupError(f"Unknown backend '{name}'; available: cadet, reference.")
