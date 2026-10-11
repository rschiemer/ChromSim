"""Parameters to estimate, addressed by their path in a project (SPEC §7.3).

A path is the project-file location joined by dots, with list positions as numbers::

    isotherm.parameters.mAb.k_eq      # isotherm parameter (plain number)
    resin.col_porosity                # plain number
    column.axial_dispersion           # quantity: give `unit`, e.g. "m^2/s"
    column.film_diffusion.mAb         # quantity in a per-component table
    method.1.inlet.mAb                # load concentration of step 1 (quantity, e.g. "g/L")

Fitting works in a scaled space x ∈ [0, 1] (SPEC §6):

- lin: θ = θ_lo + x (θ_hi − θ_lo)
- log: log10 θ = log10 θ_lo + x (log10 θ_hi − log10 θ_lo)
"""

from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from chromsim import units
from chromsim.core.project import Project
from chromsim.errors import ChromSimError
from chromsim.units import Quantity


class FitParameter(BaseModel):
    """One parameter to estimate.

    Attributes
    ----------
    path : str
        Location in the project, e.g. ``"isotherm.parameters.mAb.k_eq"``.
    lower, upper : float
        Bounds, in ``unit`` for quantities.
    scale : {"lin", "log"}
        Fit the value (lin) or its logarithm (log; for rate and equilibrium constants).
    unit : str, optional
        Unit of bounds and values for quantity fields (``column.axial_dispersion`` …).
    initial : float, optional
        Start value; default: the value in the project.
    per_experiment : bool
        One value per experiment (e.g. load concentration) instead of one shared value.
    name : str, optional
        Short name for reports; default the path.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    path: str
    lower: float
    upper: float
    scale: Literal["lin", "log"] = "lin"
    unit: str | None = None
    initial: float | None = None
    per_experiment: bool = False
    name: str | None = None

    @model_validator(mode="after")
    def _bounds(self) -> "FitParameter":
        if not self.lower < self.upper:
            raise ValueError(
                f"{self.path}: lower bound {self.lower:g} must be below upper {self.upper:g}."
            )
        if self.scale == "log" and self.lower <= 0:
            raise ValueError(f"{self.path}: log scale needs a positive lower bound.")
        if self.initial is not None and not self.lower <= self.initial <= self.upper:
            raise ValueError(f"{self.path}: initial value {self.initial:g} is outside the bounds.")
        return self

    @property
    def label(self) -> str:
        """Name used in reports."""
        return self.name or self.path

    def to_scaled(self, value: float) -> float:
        """Parameter value → x ∈ [0, 1]."""
        if self.scale == "log":
            return (np.log10(value) - np.log10(self.lower)) / (
                np.log10(self.upper) - np.log10(self.lower)
            )
        return (value - self.lower) / (self.upper - self.lower)

    def from_scaled(self, x: float) -> float:
        """x ∈ [0, 1] → parameter value."""
        if self.scale == "log":
            return 10.0 ** (
                np.log10(self.lower) + x * (np.log10(self.upper) - np.log10(self.lower))
            )
        return self.lower + x * (self.upper - self.lower)

    def derivative(self, x: float) -> float:
        """dθ/dx, for transforming uncertainties from the scaled space."""
        if self.scale == "log":
            return (
                self.from_scaled(x) * np.log(10.0) * (np.log10(self.upper) - np.log10(self.lower))
            )
        return self.upper - self.lower


def _walk(data, parts: list[str], path: str):
    """Container and key of the last path element inside nested dicts/lists."""
    node = data
    for depth, part in enumerate(parts[:-1]):
        node = _child(node, part, parts[: depth + 1], path)
    last: str | int = parts[-1]
    if isinstance(node, list):
        last = int(parts[-1])
    elif last not in node:
        available = ", ".join(map(str, node))
        raise ChromSimError(f"Parameter path '{path}': '{last}' not found. Available: {available}.")
    return node, last


def _child(node, part, so_far, path):
    try:
        return node[int(part)] if isinstance(node, list) else node[part]
    except (KeyError, IndexError, ValueError):
        options = range(len(node)) if isinstance(node, list) else list(node)
        raise ChromSimError(
            f"Parameter path '{path}': '{'.'.join(so_far)}' not found. "
            f"Available at that level: {', '.join(map(str, options))}."
        ) from None


def get_value(project: Project, path: str, unit: str | None = None) -> float:
    """Current value at ``path`` (in ``unit`` for quantities).

    Raises
    ------
    ChromSimError
        If the path does not exist, or a quantity is requested without ``unit``.
    """
    container, key = _walk(project.model_dump(mode="json"), path.split("."), path)
    value = container[key]
    if isinstance(value, str):
        if unit is None:
            example = Quantity.parse(value).unit
            raise ChromSimError(
                f"'{path}' is a quantity ({value}); give its unit, e.g. unit=\"{example}\"."
            )
        quantity = Quantity.parse(value)
        return float(units.ureg.Quantity(quantity.magnitude, quantity.unit).to(unit).magnitude)
    return float(value)


def set_values(project: Project, assignments: list[tuple[str, float, str | None]]) -> Project:
    """A validated copy of ``project`` with several values replaced.

    Parameters
    ----------
    project : Project
        Not modified.
    assignments : list of (path, value, unit)
        ``unit`` is required for quantity fields and ignored otherwise.
    """
    data = project.model_dump(mode="json")
    for path, value, unit in assignments:
        container, key = _walk(data, path.split("."), path)
        if isinstance(container[key], str) or unit is not None:
            if unit is None:
                raise ChromSimError(f"'{path}' is a quantity; give the FitParameter a unit.")
            container[key] = str(Quantity(magnitude=float(value), unit=unit))
        else:
            container[key] = float(value)
    return Project.model_validate(data)
