"""Isotherms written in Python, solved by the reference backend (SPEC §5b).

Write the binding kinetics as vectorised numpy and register the class::

    @register_isotherm
    class MyIsotherm(PythonIsotherm):
        \"\"\"<equation in LaTeX and a literature reference>\"\"\"
        name = "my_isotherm"
        parameters = {"k_eq": Parameter(...), ...}
        example = {...}                       # used by the plugin contract tests

        def flux(self, c, q, p):
            return ...                        # dq/dt, same shape as q

Shapes: ``c`` and ``q`` are (n_components, n_points). ``c`` is the liquid concentration the
binding sees (pore liquid; bulk liquid for the LRM), ``q`` the bound concentration per solid
volume; rows of non-binding components are zero and their flux must be zero. ``p`` is a
:class:`BindingParameters` with every per-component parameter as an array over components (0 for
non-binding components).

Python isotherms are always kinetic: the reference backend integrates dq/dt. For a
rapid-equilibrium model, use a large rate constant.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np

from chromsim.isotherms.base import Parameter


@dataclass
class BindingParameters:
    """Isotherm parameters prepared for vectorised flux functions.

    Attributes
    ----------
    per_component : dict[str, numpy.ndarray]
        Parameter name → value per component (shape (n_components,), 0 where not given).
    global_ : dict[str, float]
        Isotherm-wide parameters.
    ionic_capacity : float or None
        Λ, mol/m³.
    binding : numpy.ndarray
        Boolean per component: has binding parameters.
    """

    per_component: dict[str, np.ndarray]
    global_: dict[str, float] = field(default_factory=dict)
    ionic_capacity: float | None = None
    binding: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))

    def __getitem__(self, name: str) -> np.ndarray:
        """Per-component parameter as a column vector (n_components, 1), ready to broadcast."""
        return self.per_component[name][:, np.newaxis]


def binding_parameters(
    component_names: list[str],
    parameters: dict[str, dict[str, float]],
    global_parameters: dict[str, float] | None,
    ionic_capacity: float | None,
) -> BindingParameters:
    """Collect project isotherm parameters into arrays over components."""
    names = sorted({name for values in parameters.values() for name in values})
    per_component = {
        name: np.array([parameters.get(c, {}).get(name, 0.0) for c in component_names], dtype=float)
        for name in names
    }
    return BindingParameters(
        per_component=per_component,
        global_=dict(global_parameters or {}),
        ionic_capacity=ionic_capacity,
        binding=np.array([c in parameters for c in component_names]),
    )


class PythonIsotherm(ABC):
    """Base class for isotherms implemented in Python (reference backend only)."""

    name: ClassVar[str]
    parameters: ClassVar[dict[str, Parameter]]
    requires_salt: ClassVar[bool] = False
    requires_ionic_capacity: ClassVar[bool] = False
    n_modifiers: ClassVar[int] = 0
    example: ClassVar[dict]

    def bound_states(self, component_names, parameters) -> list[int]:
        """1 for components with parameters, else 0."""
        return [1 if name in parameters else 0 for name in component_names]

    def modifier_components(self, component_names) -> list[str]:
        """Leading modifier components (salt, pH) that start at their first-step value."""
        return list(component_names[: self.n_modifiers])

    def initial_bound_state(
        self, component_names, parameters, ionic_capacity_mol_per_m3
    ) -> list[float]:
        """Bound concentrations of an equilibrated, protein-free column: nothing bound."""
        return [0.0] * sum(self.bound_states(component_names, parameters))

    def is_kinetic(self, parameters) -> bool:
        """Python isotherms are integrated as kinetics."""
        return True

    @abstractmethod
    def flux(self, c: np.ndarray, q: np.ndarray, p: BindingParameters) -> np.ndarray:
        """Binding rate dq/dt, shape (n_components, n_points); zero rows for non-binding ones."""
