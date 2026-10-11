"""Unit handling at the boundary (SPEC §2.4, §11 Phase 1).

Expected values are worked out by hand in the comments, independently of pint, so a units bug
cannot hide behind the same library computing both sides.
"""

import math

import pytest

from chromsim import units
from chromsim.errors import ChromSimError, UnitError
from chromsim.units import Quantity

# Column used throughout: 20 cm long, 0.5 cm inner diameter (SPEC §4.3).
#   cross-section A = π (0.25 cm)² = 0.196350 cm²  = 1.96350e-5 m²
#   column volume V = A · 20 cm    = 3.92699 cm³    = 3.92699e-6 m³
DIAMETER_M = 0.005
LENGTH_M = 0.20
AREA_M2 = math.pi * 0.0025**2
VOLUME_M3 = AREA_M2 * LENGTH_M
# 150 cm/h × 0.196350 cm² = 29.4524 cm³/h = 0.490874 mL/min = 8.18123e-9 m³/s
FLOW_150_CM_PER_H_M3_PER_S = 1.5 / 3600 * AREA_M2


# --- Quantity parsing and formatting ------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "magnitude", "unit"),
    [
        ("20 cm", 20.0, "cm"),
        ("  45 um ", 45.0, "um"),
        ("1.5e-3 m", 1.5e-3, "m"),
        ("150 cm/h", 150.0, "cm/h"),
        ("5 CV", 5.0, "CV"),
        ("-2 mM", -2.0, "mM"),
        (".5 mL", 0.5, "mL"),
    ],
)
def test_quantity_parses_value_and_unit(text, magnitude, unit):
    quantity = Quantity.parse(text)
    assert quantity.magnitude == magnitude
    assert quantity.unit == unit


@pytest.mark.parametrize("text", ["20 cm", "0.25 mL", "5 CV", "150 cm/h", "1e-08 m^2/s", "1 g/L"])
def test_quantity_string_round_trip(text):
    assert str(Quantity.parse(text)) == text


def test_quantity_accepts_mapping_and_instance():
    assert Quantity.model_validate({"magnitude": 2, "unit": "mL"}) == Quantity.parse("2 mL")
    quantity = Quantity.parse("2 mL")
    assert Quantity.model_validate(quantity) == quantity


def test_quantity_serialises_to_string():
    assert Quantity.parse("0.25 mL").model_dump() == "0.25 mL"


@pytest.mark.parametrize(
    ("text", "hint"),
    [
        ("20", "has no unit"),
        (20, "has no unit"),
        ("cm", "does not start with a number"),
        ("20 furlongz", "unknown unit"),
    ],
)
def test_quantity_rejects_bad_input(text, hint):
    with pytest.raises(UnitError, match=hint):
        Quantity.parse(text)


def test_unit_error_is_a_chromsim_error_and_a_value_error():
    # ChromSimError for callers; ValueError so pydantic reports it as a validation error.
    assert issubclass(UnitError, ChromSimError)
    assert issubclass(UnitError, ValueError)


# --- dimension checks and SI conversion ----------------------------------------------------------


def test_check_dimension_returns_matching_name():
    assert units.check_dimension(Quantity.parse("2 mL"), "time", "volume") == "volume"
    assert units.check_dimension(Quantity.parse("5 CV"), "volume", "column_volumes") == (
        "column_volumes"
    )


def test_check_dimension_rejects_wrong_dimension_with_actionable_message():
    with pytest.raises(UnitError, match=r"'5 mL' is not a length.*e\.g\. cm"):
        units.check_dimension(Quantity.parse("5 mL"), "length")


@pytest.mark.parametrize(
    ("text", "dimension", "expected_si"),
    [
        ("20 cm", "length", 0.20),  # m
        ("45 um", "length", 45e-6),  # m
        ("0.25 mL", "volume", 0.25e-6),  # m³
        ("5 min", "time", 300.0),  # s
        ("1200 mM", "molar_concentration", 1200.0),  # mM == mol/m³ (CLAUDE.md)
        ("1 M", "molar_concentration", 1000.0),  # mol/m³
        ("150 cm/h", "velocity", 1.5 / 3600),  # m/s
        ("1 mL/min", "volumetric_flow", 1e-6 / 60),  # m³/s
        ("1e-7 m^2/s", "diffusivity", 1e-7),  # m²/s
        ("100 1/s", "rate", 100.0),  # 1/s
    ],
)
def test_to_si(text, dimension, expected_si):
    assert units.to_si(Quantity.parse(text), dimension) == pytest.approx(expected_si, rel=1e-12)


def test_to_si_refuses_column_volumes_without_context():
    with pytest.raises(UnitError, match="column volume"):
        units.to_si(Quantity.parse("5 CV"), "volume", "column_volumes")


def test_from_si():
    # 8.18123e-9 m³/s → 0.490874 mL/min
    assert units.from_si(FLOW_150_CM_PER_H_M3_PER_S, "mL/min") == pytest.approx(0.490874, rel=1e-5)
    assert units.from_si(1200.0, "mM") == pytest.approx(1200.0, rel=1e-12)


# --- concentrations -------------------------------------------------------------------------------


def test_molar_mass_from_kda():
    assert units.molar_mass_from_kda(148.0) == pytest.approx(148.0)  # 1 kDa ≙ 1 kg/mol


@pytest.mark.parametrize("text", ["1 g/L", "1 mg/mL", "1000 mg/L"])
def test_mass_concentration_converted_with_molar_mass(text):
    # 1 kg/m³ ÷ 148 kg/mol = 6.75676e-3 mol/m³
    result = units.concentration_to_si(Quantity.parse(text), molar_mass_kg_per_mol=148.0)
    assert result == pytest.approx(1.0 / 148.0, rel=1e-12)


def test_molar_concentration_ignores_molar_mass():
    assert units.concentration_to_si(Quantity.parse("50 mM")) == pytest.approx(50.0)


def test_mass_concentration_without_molar_mass_is_an_error():
    with pytest.raises(UnitError, match="molecular weight"):
        units.concentration_to_si(Quantity.parse("1 g/L"))


# --- column geometry ------------------------------------------------------------------------------


def test_column_geometry():
    assert units.column_cross_section(DIAMETER_M) == pytest.approx(1.963495e-5, rel=1e-6)
    assert units.column_volume(DIAMETER_M, LENGTH_M) == pytest.approx(3.926991e-6, rel=1e-6)


# --- flow: cm/h <-> mL/min (Phase 1 acceptance) --------------------------------------------------


def test_linear_velocity_to_volumetric_flow():
    flow_m3_per_s = units.flow_rate_to_si(Quantity.parse("150 cm/h"), AREA_M2)
    assert flow_m3_per_s == pytest.approx(8.181231e-9, rel=1e-6)
    assert units.from_si(flow_m3_per_s, "mL/min") == pytest.approx(0.490874, rel=1e-5)


def test_volumetric_flow_to_linear_velocity():
    # 0.490874 mL/min on a 0.5 cm column → 150 cm/h
    flow_m3_per_s = units.flow_rate_to_si(Quantity.parse("0.490874 mL/min"), AREA_M2)
    velocity_m_per_s = units.flow_rate_to_linear_velocity(flow_m3_per_s, AREA_M2)
    assert units.from_si(velocity_m_per_s, "cm/h") == pytest.approx(150.0, rel=1e-5)


def test_volumetric_flow_does_not_depend_on_column():
    flow = Quantity.parse("1 mL/min")
    assert units.flow_rate_to_si(flow, AREA_M2) == units.flow_rate_to_si(flow, 1.0)


@pytest.mark.parametrize("unit", ["cm/h", "mL/min"])
def test_convert_flow_round_trip(unit):
    other = "mL/min" if unit == "cm/h" else "cm/h"
    original = Quantity(magnitude=150.0, unit=unit)
    there = units.convert_flow(original, other, AREA_M2)
    back = units.convert_flow(there, unit, AREA_M2)
    assert back.magnitude == pytest.approx(150.0, rel=1e-12)


def test_flow_rejects_non_flow_units():
    with pytest.raises(UnitError, match="not a flow"):
        units.flow_rate_to_si(Quantity.parse("150 cm"), AREA_M2)


# --- durations: CV / mL / min (Phase 1 acceptance) -----------------------------------------------


def test_duration_in_column_volumes():
    # 1 CV takes L/u = 20 cm ÷ 150 cm/h = 0.1333 h = 480 s, so 5 CV = 2400 s.
    seconds = units.duration_to_si(
        Quantity.parse("5 CV"), FLOW_150_CM_PER_H_M3_PER_S, column_volume_m3=VOLUME_M3
    )
    assert seconds == pytest.approx(2400.0, rel=1e-12)


def test_duration_in_volume():
    # 2 mL ÷ 29.4524 mL/h = 0.067906 h = 244.462 s
    seconds = units.duration_to_si(
        Quantity.parse("2 mL"), FLOW_150_CM_PER_H_M3_PER_S, column_volume_m3=VOLUME_M3
    )
    assert seconds == pytest.approx(244.462, rel=1e-5)


def test_duration_in_time_ignores_flow():
    seconds = units.duration_to_si(Quantity.parse("10 min"), 123.0, column_volume_m3=VOLUME_M3)
    assert seconds == pytest.approx(600.0)


def test_column_volumes_with_volumetric_flow():
    # 1 CV = 3.92699 mL at 0.5 mL/min → 7.85398 min = 471.239 s
    flow_m3_per_s = units.flow_rate_to_si(Quantity.parse("0.5 mL/min"), AREA_M2)
    seconds = units.duration_to_si(Quantity.parse("1 CV"), flow_m3_per_s, VOLUME_M3)
    assert seconds == pytest.approx(471.239, rel=1e-5)


def test_volume_based_duration_with_zero_flow_is_an_error():
    with pytest.raises(UnitError, match="flow rate is zero"):
        units.duration_to_si(Quantity.parse("5 CV"), 0.0, VOLUME_M3)


def test_duration_rejects_wrong_dimension():
    with pytest.raises(UnitError, match="not a duration"):
        units.duration_to_si(Quantity.parse("5 mM"), 1e-8, VOLUME_M3)


def test_ph_pseudo_concentration():
    # pH is carried as a pseudo-component; its "concentration" is the pH number itself.
    ph = Quantity.parse("5.5 pH")
    assert units.concentration_to_si(ph) == 5.5
    assert units.concentration_from_si(5.5, "pH") == 5.5
    with pytest.raises(UnitError, match="not a length"):
        units.check_dimension(ph, "length")
