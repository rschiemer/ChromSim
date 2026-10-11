"""Units at the boundary, SI inside (SPEC §2.4).

This is the only module in ChromSim that converts units. Users write quantities as text with a
unit, for example ``"20 cm"``, ``"150 cm/h"``, ``"5 CV"`` or ``"1 g/L"``. Everything after this
module works in SI units as CADET expects:

==========================  ==================  ====================================
quantity                    SI unit             common user units
==========================  ==================  ====================================
length                      m                   cm, mm, um
volume                      m³                  mL, L
time                        s                   s, min, h
linear (superficial) vel.   m/s                 cm/h
volumetric flow rate        m³/s                mL/min
molar concentration         mol/m³              mM (identical to mol/m³), M
mass concentration          kg/m³               g/L, mg/mL (identical to g/L)
diffusivity                 m²/s                m^2/s, cm^2/s
film transfer coefficient   m/s                 m/s, cm/s
rate constant               1/s                 1/s, 1/min
specific absorptivity       m²/kg               L/(g*cm)
conductivity                S/m                 mS/cm
absorbance (detector)       AU                  mAU
==========================  ==================  ====================================

Column volumes (``CV``) are not a pint unit, because converting them needs the column
volume. :func:`duration_to_si` handles them explicitly.

Equations
---------
Column cross-section and packed-bed volume (CLAUDE.md, domain conventions):

.. math:: A = \\pi (d/2)^2, \\qquad V_c = A \\, L

Volumetric flow rate from linear (superficial) velocity:

.. math:: Q = u \\, A
"""

import math
import re
from typing import Annotated, Any, overload

import numpy as np
import pint
from pydantic import (
    AfterValidator,
    BaseModel,
    ConfigDict,
    ValidationError,
    field_validator,
    model_serializer,
    model_validator,
)

from chromsim.errors import UnitError

# The single pint registry used by ChromSim. It is created once and never modified afterwards.
ureg = pint.UnitRegistry()

COLUMN_VOLUME_UNIT = "CV"
# pH as a pseudo-component "concentration" (CLAUDE.md: placed directly after salt). Its value is
# the pH number itself and is passed to CADET unchanged.
PH_UNIT = "pH"

# Dimension name -> SI unit. The SI unit also defines the dimension (e.g. m**3/s = volume/time).
SI_UNITS = {
    "length": "m",
    "volume": "m**3",
    "time": "s",
    "velocity": "m/s",
    "volumetric_flow": "m**3/s",
    "molar_concentration": "mol/m**3",
    "mass_concentration": "kg/m**3",
    "diffusivity": "m**2/s",
    "rate": "1/s",
    "specific_absorptivity": "m**2/kg",
    "conductivity": "S/m",
    "pressure": "Pa",
    "dimensionless": "dimensionless",
}

# Dimension name -> units suggested in error messages.
EXAMPLE_UNITS = {
    "length": "cm or mm",
    "volume": "mL",
    "time": "min or s",
    "velocity": "cm/h",
    "volumetric_flow": "mL/min",
    "molar_concentration": "mM",
    "mass_concentration": "g/L",
    "diffusivity": "m^2/s",
    "rate": "1/s",
    "specific_absorptivity": "L/(g*cm)",
    "conductivity": "mS/cm",
    "pressure": "MPa or bar",
    "dimensionless": "dimensionless",
    "column_volumes": "CV",
    "ph": "pH",
}

_NUMBER_THEN_UNIT = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$")


def format_magnitude(value: float) -> str:
    """Format a number without losing precision: ``20.0 -> "20"``, ``0.25 -> "0.25"``."""
    if value.is_integer() and abs(value) < 1e15:
        return str(int(value))
    return repr(value)


class Quantity(BaseModel):
    """A number with a unit, exactly as the user wrote it.

    The unit text is kept unchanged so project files round-trip without loss; convert to SI with
    the functions in this module.

    Attributes
    ----------
    magnitude : float
        The number, in ``unit``.
    unit : str
        Any unit pint understands (``"cm"``, ``"mL/min"``, ``"mM"`` …) or ``"CV"``.

    Examples
    --------
    >>> Quantity.parse("150 cm/h")
    Quantity(magnitude=150.0, unit='cm/h')
    >>> str(Quantity(magnitude=0.25, unit="mL"))
    '0.25 mL'
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    magnitude: float
    unit: str

    @model_validator(mode="before")
    @classmethod
    def _split_text(cls, value: Any) -> Any:
        if isinstance(value, str):
            match = _NUMBER_THEN_UNIT.match(value)
            if match is None:
                raise UnitError(
                    f"'{value}' does not start with a number; "
                    "write a value and a unit, e.g. '20 cm'."
                )
            number, unit = match.groups()
            if not unit:
                raise UnitError(
                    f"'{value}' has no unit; write e.g. '{number} cm' or '{number} mL'."
                )
            return {"magnitude": float(number), "unit": unit}
        if isinstance(value, int | float):
            raise UnitError(f"'{value}' has no unit; write e.g. '{value} cm' or '{value} mL'.")
        return value

    @field_validator("unit")
    @classmethod
    def _known_unit(cls, unit: str) -> str:
        unit = unit.strip()
        if unit in (COLUMN_VOLUME_UNIT, PH_UNIT):
            return unit
        try:
            ureg.Quantity(1.0, unit)
        except (pint.errors.PintError, ValueError, SyntaxError, TypeError) as error:
            raise UnitError(
                f"unknown unit '{unit}'; use e.g. cm, mL, min, cm/h, mL/min, mM or g/L."
            ) from error
        return unit

    @model_serializer
    def _to_text(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"{format_magnitude(self.magnitude)} {self.unit}"

    @classmethod
    def parse(cls, value: "str | dict | Quantity") -> "Quantity":
        """Create a Quantity from text such as ``"20 cm"``.

        Raises
        ------
        UnitError
            If the text has no number, no unit, or an unknown unit.
        """
        try:
            return cls.model_validate(value)
        except ValidationError as error:
            first = error.errors()[0]
            original = first.get("ctx", {}).get("error")
            if isinstance(original, UnitError):
                raise original from None
            raise UnitError(f"cannot read quantity {value!r}: {first['msg']}") from None


# --- dimension checks -----------------------------------------------------------------------------


def _describe(dimension: str) -> str:
    return dimension.replace("_", " ")


def check_dimension(quantity: Quantity, *allowed: str, what: str | None = None) -> str:
    """Check that a quantity has one of the allowed dimensions and return which one.

    Parameters
    ----------
    quantity : Quantity
        The quantity to check.
    *allowed : str
        Dimension names: keys of ``SI_UNITS`` or ``"column_volumes"``.
    what : str, optional
        What the quantity is, for the error message (e.g. ``"duration"``). Defaults to the
        allowed dimension names.

    Returns
    -------
    str
        The matching dimension name.

    Raises
    ------
    UnitError
        If the quantity has none of the allowed dimensions.
    """
    # Values assigned into dicts of a model (e.g. ``step.inlet["mAb"] = "1 g/L"``) are not
    # validated by pydantic, so text may arrive here; accept it like any other input.
    quantity = Quantity.parse(quantity)
    if quantity.unit == COLUMN_VOLUME_UNIT:
        if "column_volumes" in allowed:
            return "column_volumes"
    elif quantity.unit == PH_UNIT:
        if "ph" in allowed:
            return "ph"
    else:
        dimensionality = ureg.Quantity(1.0, quantity.unit).dimensionality
        for dimension in allowed:
            if dimension in ("column_volumes", "ph"):
                continue
            if dimensionality == ureg.Quantity(1.0, SI_UNITS[dimension]).dimensionality:
                return dimension

    description = what or " or ".join(_describe(dimension) for dimension in allowed)
    examples = ", ".join(EXAMPLE_UNITS[dimension] for dimension in allowed)
    raise UnitError(f"'{quantity}' is not a {description}; use e.g. {examples}.")


# --- conversion to and from SI --------------------------------------------------------------------


def to_si(quantity: Quantity, *allowed: str) -> float:
    """Convert a quantity to the SI unit of its dimension.

    Parameters
    ----------
    quantity : Quantity
        Quantity in any unit.
    *allowed : str
        Allowed dimension names (see ``SI_UNITS``).

    Returns
    -------
    float
        Magnitude in SI units (e.g. m for a length, mol/m³ for a molar concentration).

    Raises
    ------
    UnitError
        If the dimension is not allowed, or the quantity is in column volumes (use
        :func:`duration_to_si`, which knows the column volume and flow rate).
    """
    quantity = Quantity.parse(quantity)
    dimension = check_dimension(quantity, *allowed)
    if dimension == "column_volumes":
        raise UnitError(
            f"'{quantity}' is in column volumes; converting it needs the column volume and flow "
            "rate, so use duration_to_si()."
        )
    return float(ureg.Quantity(quantity.magnitude, quantity.unit).to(SI_UNITS[dimension]).magnitude)


@overload
def from_si(value_si: float, unit: str) -> float: ...


@overload
def from_si(value_si: np.ndarray, unit: str) -> np.ndarray: ...


def from_si(value_si, unit: str):
    """Convert an SI value into ``unit``, e.g. ``from_si(8.18e-9, "mL/min") -> 0.491``.

    Parameters
    ----------
    value_si : float or numpy.ndarray
        Value(s) in the SI unit of the same dimension as ``unit``.
    unit : str
        Target unit understood by pint.

    Returns
    -------
    float or numpy.ndarray
        Value(s) expressed in ``unit``.
    """
    si_unit = ureg.Quantity(1.0, unit).to_base_units().units
    converted = ureg.Quantity(value_si, si_unit).to(unit).magnitude
    return float(converted) if np.ndim(converted) == 0 else np.asarray(converted, dtype=float)


def molar_mass_from_kda(mw_kDa: float) -> float:
    """Convert a molecular weight in kDa to a molar mass in kg/mol (1 kDa ≙ 1 kg/mol).

    Parameters
    ----------
    mw_kDa : float
        Molecular weight, kDa.

    Returns
    -------
    float
        Molar mass, kg/mol.
    """
    return mw_kDa * 1.0  # 1 kDa = 1000 g/mol = 1 kg/mol


def concentration_to_si(quantity: Quantity, molar_mass_kg_per_mol: float | None = None) -> float:
    """Convert a molar or mass concentration to mol/m³.

    Parameters
    ----------
    quantity : Quantity
        Concentration, e.g. ``"50 mM"`` or ``"1 g/L"``.
    molar_mass_kg_per_mol : float, optional
        Molar mass of the component, kg/mol. Required for mass concentrations.

    Returns
    -------
    float
        Concentration, mol/m³ (identical to mM).

    Raises
    ------
    UnitError
        If the quantity is not a concentration, or is a mass concentration and no molar mass
        is given.
    """
    quantity = Quantity.parse(quantity)
    dimension = check_dimension(
        quantity, "molar_concentration", "mass_concentration", "ph", what="concentration"
    )
    if dimension == "ph":
        return quantity.magnitude
    if dimension == "molar_concentration":
        return to_si(quantity, "molar_concentration")
    if molar_mass_kg_per_mol is None:
        raise UnitError(
            f"'{quantity}' is a mass concentration; converting it to mol/m³ needs the component's "
            "molecular weight. Set `mw_kDa` for this component or give the concentration in mM."
        )
    return to_si(quantity, "mass_concentration") / molar_mass_kg_per_mol


def concentration_from_si(value_mol_per_m3, unit: str, molar_mass_kg_per_mol: float | None = None):
    """Convert a concentration in mol/m³ to ``unit`` (molar, e.g. mM, or mass, e.g. g/L).

    Parameters
    ----------
    value_mol_per_m3 : float or numpy.ndarray
        Concentration, mol/m³.
    unit : str
        Target unit.
    molar_mass_kg_per_mol : float, optional
        Molar mass, kg/mol. Required for mass-concentration units.

    Returns
    -------
    float or numpy.ndarray
        Concentration in ``unit``.
    """
    dimension = check_dimension(
        Quantity(magnitude=1.0, unit=unit),
        "molar_concentration",
        "mass_concentration",
        "ph",
        what="concentration",
    )
    if dimension == "ph":
        return value_mol_per_m3 * 1.0
    if dimension == "molar_concentration":
        return value_mol_per_m3 * from_si(1.0, unit)
    if molar_mass_kg_per_mol is None:
        raise UnitError(
            f"Converting to {unit} needs the component's molecular weight; set `mw_kDa` or use mM."
        )
    return value_mol_per_m3 * molar_mass_kg_per_mol * from_si(1.0, unit)


def column_cross_section(diameter_m: float) -> float:
    """Column cross-section :math:`A = \\pi (d/2)^2`, m².

    Parameters
    ----------
    diameter_m : float
        Inner column diameter d, m.
    """
    d = diameter_m
    return math.pi * (d / 2.0) ** 2


def column_volume(diameter_m: float, length_m: float) -> float:
    """Packed-bed column volume :math:`V_c = \\pi (d/2)^2 L`, m³ (1 CV).

    Parameters
    ----------
    diameter_m : float
        Inner column diameter d, m.
    length_m : float
        Packed-bed length L, m.
    """
    return column_cross_section(diameter_m) * length_m


def linear_velocity_to_flow_rate(velocity_m_per_s: float, cross_section_m2: float) -> float:
    """Volumetric flow rate :math:`Q = u A`, m³/s, from superficial velocity u (m/s) and A (m²)."""
    return velocity_m_per_s * cross_section_m2


def flow_rate_to_linear_velocity(flow_rate_m3_per_s: float, cross_section_m2: float) -> float:
    """Superficial velocity :math:`u = Q / A`, m/s, from flow rate Q (m³/s) and A (m²)."""
    return flow_rate_m3_per_s / cross_section_m2


def flow_rate_to_si(flow: Quantity, cross_section_m2: float) -> float:
    """Convert a flow given as linear velocity (cm/h) or volumetric flow (mL/min) to m³/s.

    Parameters
    ----------
    flow : Quantity
        Linear (superficial) velocity or volumetric flow rate.
    cross_section_m2 : float
        Column cross-section, m². Only used for linear velocities.

    Returns
    -------
    float
        Volumetric flow rate, m³/s.
    """
    dimension = check_dimension(flow, "velocity", "volumetric_flow", what="flow rate")
    if dimension == "velocity":
        return linear_velocity_to_flow_rate(to_si(flow, "velocity"), cross_section_m2)
    return to_si(flow, "volumetric_flow")


def convert_flow(flow: Quantity, to_unit: str, cross_section_m2: float) -> Quantity:
    """Convert between linear velocity and volumetric flow, e.g. cm/h ↔ mL/min.

    Parameters
    ----------
    flow : Quantity
        Linear velocity or volumetric flow rate.
    to_unit : str
        Target unit, a velocity (e.g. ``"cm/h"``) or a volumetric flow (e.g. ``"mL/min"``).
    cross_section_m2 : float
        Column cross-section, m².

    Returns
    -------
    Quantity
        The flow expressed in ``to_unit``.
    """
    flow_rate_m3_per_s = flow_rate_to_si(flow, cross_section_m2)
    target = Quantity(magnitude=1.0, unit=to_unit)
    if check_dimension(target, "velocity", "volumetric_flow", what="flow rate") == "velocity":
        velocity = flow_rate_to_linear_velocity(flow_rate_m3_per_s, cross_section_m2)
        return Quantity(magnitude=from_si(velocity, to_unit), unit=to_unit)
    return Quantity(magnitude=from_si(flow_rate_m3_per_s, to_unit), unit=to_unit)


def duration_to_si(duration: Quantity, flow_rate_m3_per_s: float, column_volume_m3: float) -> float:
    """Convert a step duration given in time, volume or column volumes to seconds.

    Parameters
    ----------
    duration : Quantity
        Duration as time (``"10 min"``), volume (``"2 mL"``) or column volumes (``"5 CV"``).
    flow_rate_m3_per_s : float
        Volumetric flow rate during the step, m³/s. Not used for time units.
    column_volume_m3 : float
        Packed-bed column volume, m³. Only used for CV.

    Returns
    -------
    float
        Duration, s.

    Raises
    ------
    UnitError
        If the quantity is not a duration, or it is volume-based and the flow rate is not positive.
    """
    duration = Quantity.parse(duration)
    dimension = check_dimension(duration, "time", "volume", "column_volumes", what="duration")
    if dimension == "time":
        return to_si(duration, "time")

    if dimension == "column_volumes":
        volume_m3 = duration.magnitude * column_volume_m3
    else:
        volume_m3 = to_si(duration, "volume")
    if flow_rate_m3_per_s <= 0.0:
        raise UnitError(
            f"'{duration}' is given as a volume, but the flow rate is zero or negative; give the "
            "duration in time units (e.g. min) instead."
        )
    return volume_m3 / flow_rate_m3_per_s


# --- detector signals -----------------------------------------------------------------------------


def signal_to_si(values, unit: str):
    """Convert detector values (UV, conductivity, %B, pH, pressure …) to ChromSim's internal units.

    Parameters
    ----------
    values : numpy.ndarray
        Values as exported, in ``unit``.
    unit : str
        Unit text from the export, e.g. ``"mAU"``, ``"mS/cm"``, ``"%"``, ``"pH"``, ``"MPa"``.

    Returns
    -------
    tuple of (numpy.ndarray, str)
        Converted values and their unit: absorbance in ``"AU"``, conductivity in ``"S/m"``,
        percentages as a fraction (``"1"``), pH as ``"pH"``, other pint-compatible units in SI.
        Units ChromSim does not know (e.g. offset temperatures) are returned unchanged.
    """
    text = unit.strip()
    if text in ("", "pH"):
        return values * 1.0, "pH"
    if text in ("%", "percent"):
        return values / 100.0, "1"
    if text.endswith("AU"):  # mAU, AU (absorbance units, not astronomical units)
        return values * ureg.Quantity(1.0, text).to("absorbance_unit").magnitude, "AU"
    try:
        quantity = ureg.Quantity(1.0, text)
    except (pint.errors.PintError, ValueError, SyntaxError, TypeError):
        return values * 1.0, text
    for dimension in ("conductivity", "pressure", "volume", "time", "volumetric_flow"):
        si_unit = SI_UNITS[dimension]
        if quantity.dimensionality == ureg.Quantity(1.0, si_unit).dimensionality:
            return values * quantity.to(si_unit).magnitude, si_unit
    return values * 1.0, text


def signal_from_si(values_si, si_unit: str, unit: str):
    """Inverse of :func:`signal_to_si` for display, e.g. S/m → mS/cm or AU → mAU."""
    if si_unit == "pH" or si_unit == unit:
        return values_si * 1.0
    if si_unit == "1" and unit in ("%", "percent"):
        return values_si * 100.0
    target = "absorbance_unit" if si_unit == "AU" else si_unit
    return values_si * ureg.Quantity(1.0, target).to(unit).magnitude


# --- annotated types for pydantic models ---------------------------------------------------------


def _require(*allowed: str, what: str | None = None, allow_zero: bool = True):
    """Build a pydantic validator that checks dimension and sign of a Quantity field."""

    def validate(quantity: Quantity) -> Quantity:
        check_dimension(quantity, *allowed, what=what)
        if quantity.magnitude < 0 or (quantity.magnitude == 0 and not allow_zero):
            sign = "zero or negative" if not allow_zero else "negative"
            raise UnitError(f"'{quantity}' must not be {sign}.")
        return quantity

    return validate


Length = Annotated[Quantity, AfterValidator(_require("length", allow_zero=False))]
Volume = Annotated[Quantity, AfterValidator(_require("volume"))]
Duration = Annotated[
    Quantity,
    AfterValidator(_require("time", "volume", "column_volumes", what="duration", allow_zero=False)),
]
FlowRate = Annotated[
    Quantity, AfterValidator(_require("velocity", "volumetric_flow", what="flow rate"))
]
Concentration = Annotated[
    Quantity,
    AfterValidator(
        _require("molar_concentration", "mass_concentration", "ph", what="concentration")
    ),
]
MolarConcentration = Annotated[Quantity, AfterValidator(_require("molar_concentration"))]
Diffusivity = Annotated[Quantity, AfterValidator(_require("diffusivity"))]
# Film mass-transfer coefficient k_f has the dimension of a velocity (m/s).
FilmTransferCoefficient = Annotated[
    Quantity, AfterValidator(_require("velocity", what="film transfer coefficient"))
]
