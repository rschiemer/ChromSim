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
        Whether every binding component must give it.
    """

    unit: str
    description: str
    bounds: tuple[float | None, float | None] = (0.0, None)
    log_scale: bool = False
    required: bool = True


class IsothermBase(Protocol):
    """What the registry and the project checks need from any isotherm."""

    name: ClassVar[str]
    parameters: ClassVar[dict[str, Parameter]]
    requires_salt: ClassVar[bool]
    requires_ionic_capacity: ClassVar[bool]


def _bounds_text(bounds: tuple[float | None, float | None]) -> str:
    lower, upper = bounds
    if lower is not None and upper is not None:
        return f"between {lower:g} and {upper:g}"
    if lower is not None:
        return f"at least {lower:g}"
    return f"at most {upper:g}"


def validate_parameters(
    isotherm: IsothermBase, component_names: list[str], parameters: dict[str, dict[str, float]]
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

    Raises
    ------
    IsothermConfigError
        On a misplaced salt component, unknown or missing parameters, values out of bounds, or
        ``k_kin`` given for only some binding components.
    """
    if isotherm.requires_salt:
        for index, name in enumerate(component_names):
            if name.lower() == "salt" and index != 0:
                raise IsothermConfigError(
                    f"Salt must be component 0 for {isotherm.name}; found '{component_names[0]}' "
                    "at index 0. Move salt to the top of the components list."
                )
        if component_names[0] in parameters:
            raise IsothermConfigError(
                f"component 0 ('{component_names[0]}') is the salt for {isotherm.name} and cannot "
                "have binding parameters; put the salt first and remove its parameters."
            )

    allowed = ", ".join(isotherm.parameters)
    for component_name, values in parameters.items():
        for parameter_name, value in values.items():
            if parameter_name not in isotherm.parameters:
                raise IsothermConfigError(
                    f"{isotherm.name}, component '{component_name}': unknown parameter "
                    f"'{parameter_name}'. {isotherm.name} parameters: {allowed}."
                )
            lower, upper = isotherm.parameters[parameter_name].bounds
            if (lower is not None and value < lower) or (upper is not None and value > upper):
                raise IsothermConfigError(
                    f"{isotherm.name}, component '{component_name}': {parameter_name} = {value:g} "
                    f"is outside its bounds; it must be "
                    f"{_bounds_text(isotherm.parameters[parameter_name].bounds)}."
                )
        missing = [
            name
            for name, parameter in isotherm.parameters.items()
            if parameter.required and name not in values
        ]
        if missing:
            raise IsothermConfigError(
                f"{isotherm.name}, component '{component_name}': missing parameter(s) "
                f"{', '.join(missing)}."
            )

    if "k_kin" in isotherm.parameters and parameters:
        with_k_kin = [name for name, values in parameters.items() if "k_kin" in values]
        if with_k_kin and len(with_k_kin) != len(parameters):
            raise IsothermConfigError(
                f"{isotherm.name}: give k_kin for all binding components or none (none means "
                f"rapid equilibrium); it is given only for {', '.join(with_k_kin)}."
            )
