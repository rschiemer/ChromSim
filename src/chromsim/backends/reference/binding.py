"""Binding: dq/dt from the isotherm's flux function, with shape and sanity checks."""

import numpy as np

from chromsim.errors import SimulationError, SimulationSetupError
from chromsim.isotherms.python_isotherm import BindingParameters


def check_reference_capable(isotherm, kinetic: bool) -> None:
    """The reference backend needs a ``flux`` method and kinetic binding.

    Raises
    ------
    SimulationSetupError
        With what to change.
    """
    if not hasattr(isotherm, "flux"):
        raise SimulationSetupError(
            f"Isotherm '{isotherm.name}' has no Python flux function, so the reference backend "
            "cannot simulate it; use `solver: {backend: cadet}`."
        )
    if not kinetic:
        raise SimulationSetupError(
            f"The reference backend integrates binding kinetics: give k_kin for every binding "
            f"component of '{isotherm.name}' (a large k_kin approximates rapid equilibrium)."
        )


def binding_rate(
    isotherm, c: np.ndarray, q: np.ndarray, parameters: BindingParameters
) -> np.ndarray:
    """dq/dt from the isotherm, shape (n_components, n_cells)."""
    rate = isotherm.flux(c, q, parameters)
    if rate.shape != q.shape:
        raise SimulationError(
            f"Isotherm '{isotherm.name}' returned dq/dt of shape {rate.shape}, expected {q.shape}."
        )
    return rate
