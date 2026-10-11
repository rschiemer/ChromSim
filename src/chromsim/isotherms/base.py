"""Building blocks shared by all isotherms: parameter descriptions and parameter checks.

An isotherm is described by a class with

- ``name``: the name used in project files (``isotherm: {model: SMA}``),
- ``parameters``: ``{parameter name: Parameter}``, one set per binding component,
- ``requires_salt``: whether component 0 must be the salt,
- ``requires_ionic_capacity``: whether the resin's ionic capacity Λ is used.

Components that have a parameter set in the project bind; all others are non-binding.
"""

from dataclasses import dataclass
from typing import ClassVar, Protocol

from chromsim.errors import IsothermConfigError


@dataclass(frozen=True)
class Parameter:
    """Description of one isotherm parameter.

    Attributes
    ----------
    unit : str
        Unit of the value in project files. ChromSim isotherms use SI-compatible units
        (``mol/m^3`` = mM, ``m^3/mol`` = 1/mM, ``1/s``), so values are passed on unchanged.
    description : str
        Short meaning, shown in docs and UI forms.
    bounds : tuple of (float or None, float or None)
        Allowed range (inclusive); ``None`` means unbounded.
    log_scale : bool
        Whether fitting should work on log10 of the value.
    required : bool
        Whether it must be given (by every binding component, if per component).
    per_component : bool
        One value per binding component (``isotherm.parameters``), or one value for the whole
        isotherm (``isotherm.global_parameters``).
    """

    unit: str
    description: str
    bounds: tuple[float | None, float | None] = (0.0, None)
    log_scale: bool = False
    required: bool = True
    per_component: bool = True


class IsothermBase(Protocol):
    """What the registry and the project checks need from any isotherm."""

    name: ClassVar[str]
    parameters: ClassVar[dict[str, Parameter]]
    requires_salt: ClassVar[bool]
    requires_ionic_capacity: ClassVar[bool]
    # Leading components that modulate binding and never bind themselves: salt (1) or
    # salt and the pH pseudo-component (2).
    n_modifiers: ClassVar[int]

    def bound_states(self, component_names: list[str], parameters: dict) -> list[int]:
        """Number of bound states per component (per site for multi-site isotherms)."""
        ...

    def modifier_components(self, component_names: list[str]) -> list[str]:
        """Leading modifier components that start at their first-step inlet value."""
        ...

    def initial_bound_state(
        self, component_names: list[str], parameters: dict, ionic_capacity_mol_per_m3: float | None
    ) -> list[float]:
        """Bound concentrations of an equilibrated, protein-free column, one per bound state."""
        ...

    def is_kinetic(self, parameters: dict) -> bool:
        """Whether binding is kinetic (else rapid equilibrium)."""
        ...


def _bounds_text(bounds: tuple[float | None, float | None]) -> str:
    lower, upper = bounds
    if lower is not None and upper is not None:
        return f"between {lower:g} and {upper:g}"
    if lower is not None:
        return f"at least {lower:g}"
    return f"at most {upper:g}"


def _check_value(isotherm, owner: str, name: str, value: float) -> None:
    lower, upper = isotherm.parameters[name].bounds
    if (lower is not None and value < lower) or (upper is not None and value > upper):
        raise IsothermConfigError(
            f"{isotherm.name}, {owner}: {name} = {value:g} is outside its bounds; it must be "
            f"{_bounds_text(isotherm.parameters[name].bounds)}."
        )


def validate_parameters(
    isotherm: IsothermBase,
    component_names: list[str],
    parameters: dict[str, dict[str, float]],
    global_parameters: dict[str, float] | None = None,
) -> None:
    """Check isotherm parameters against the isotherm's description.

    Parameters
    ----------
    isotherm : IsothermBase
        The isotherm.
    component_names : list of str
        Project components in order.
    parameters : dict
        Component name → {parameter name → value}, as in the project file.
    global_parameters : dict, optional
        Parameter name → value for parameters that are not per component.

    Raises
    ------
    IsothermConfigError
        On misplaced or binding modifier components (salt, pH), unknown or missing parameters,
        values out of bounds, or ``k_kin`` given for only some binding components.
    """
    global_parameters = global_parameters or {}
    if isotherm.requires_salt:
        for index, name in enumerate(component_names):
            if name.lower() == "salt" and index != 0:
                raise IsothermConfigError(
                    f"Salt must be component 0 for {isotherm.name}; found '{component_names[0]}' "
                    "at index 0. Move salt to the top of the components list."
                )
    n_modifiers = getattr(isotherm, "n_modifiers", 1 if isotherm.requires_salt else 0)
    roles = ["the salt", "the pH pseudo-component"]
    if len(component_names) < n_modifiers + 1:
        raise IsothermConfigError(
            f"{isotherm.name} needs {' and '.join(roles[:n_modifiers])} as the first component(s), "
            "followed by at least one binding component."
        )
    for index in range(n_modifiers):
        if component_names[index] in parameters:
            raise IsothermConfigError(
                f"component {index} ('{component_names[index]}') is {roles[index]} for "
                f"{isotherm.name} and cannot have binding parameters; put "
                f"{' and '.join(roles[:n_modifiers])} first and remove their parameters."
            )

    per_component = [n for n, p in isotherm.parameters.items() if p.per_component]
    global_names = [n for n, p in isotherm.parameters.items() if not p.per_component]
    for component_name, values in parameters.items():
        for parameter_name, value in values.items():
            if parameter_name not in per_component:
                where = " (it is a global parameter)" if parameter_name in global_names else ""
                raise IsothermConfigError(
                    f"{isotherm.name}, component '{component_name}': unknown parameter "
                    f"'{parameter_name}'{where}. {isotherm.name} parameters per component: "
                    f"{', '.join(per_component)}."
                )
            _check_value(isotherm, f"component '{component_name}'", parameter_name, value)
        missing = [
            name
            for name in per_component
            if isotherm.parameters[name].required and name not in values
        ]
        if missing:
            raise IsothermConfigError(
                f"{isotherm.name}, component '{component_name}': missing parameter(s) "
                f"{', '.join(missing)}."
            )

    for parameter_name, value in global_parameters.items():
        if parameter_name not in global_names:
            raise IsothermConfigError(
                f"{isotherm.name}: unknown global parameter '{parameter_name}'. Global parameters: "
                f"{', '.join(global_names) or 'none'}."
            )
        _check_value(isotherm, "global", parameter_name, value)
    missing_global = [
        name
        for name in global_names
        if isotherm.parameters[name].required and name not in global_parameters
    ]
    if missing_global:
        raise IsothermConfigError(
            f"{isotherm.name}: missing global parameter(s) {', '.join(missing_global)}; add them "
            "under isotherm.global_parameters."
        )

    # Kinetic parameters (k_kin, or k_kin_1, k_kin_2 … per binding site): all or none.
    kinetic_names = [name for name in isotherm.parameters if name.startswith("k_kin")]
    if kinetic_names and parameters:
        given = {
            name: [p for p in kinetic_names if p in values] for name, values in parameters.items()
        }
        complete = [name for name, present in given.items() if len(present) == len(kinetic_names)]
        partial = [name for name, present in given.items() if present and name not in complete]
        if partial or (complete and len(complete) != len(parameters)):
            raise IsothermConfigError(
                f"{isotherm.name}: give {', '.join(kinetic_names)} for all binding components or "
                f"none (none means rapid equilibrium); given for {', '.join(complete + partial)}."
            )
