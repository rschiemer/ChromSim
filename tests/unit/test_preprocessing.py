"""Signal preprocessing: baselines, Lambert–Beer, conductivity calibration, axes (SPEC §7.2)."""

import numpy as np
import pytest
from project_examples import make_project
from pydantic import ValidationError

from chromsim.backends.base import pumped_volume
from chromsim.data.experiment import Signal
from chromsim.data.preprocessing import (
    ConductivityCalibration,
    conductivity_to_salt,
    load_conductivity_calibration,
    resample,
    subtract_baseline,
    uv_to_concentration,
    volume_to_time,
)
from chromsim.errors import ChromSimError

CALIBRATION = ConductivityCalibration(
    points=[("0 mM", "1 mS/cm"), ("100 mM", "11 mS/cm"), ("500 mM", "41 mS/cm")]
)


def signal(values, unit="AU", x=None):
    values = np.asarray(values, dtype=float)
    x = np.arange(values.size, dtype=float) if x is None else np.asarray(x, dtype=float)
    return Signal("s", "uv", "volume", x, values, unit, unit)


def test_calibration_interpolates_and_extrapolates():
    # 1 mS/cm = 0.1 S/m. Segment 0–100 mM: 0.1 → 1.1 S/m; segment 100–500 mM: 1.1 → 4.1 S/m.
    assert CALIBRATION.salt_from_conductivity([0.6, 1.1, 2.6]) == pytest.approx(
        [50.0, 100.0, 300.0]
    )
    assert CALIBRATION.salt_from_conductivity(4.85) == pytest.approx(600.0)  # extrapolated
    salt = np.array([0.0, 50.0, 250.0, 700.0])
    assert CALIBRATION.salt_from_conductivity(
        CALIBRATION.conductivity_from_salt(salt)
    ) == pytest.approx(salt)


@pytest.mark.parametrize(
    ("points", "message"),
    [
        ([("0 mM", "1 mS/cm")], "at least two"),
        ([("100 mM", "11 mS/cm"), ("0 mM", "1 mS/cm")], "must increase"),
        ([("0 mM", "1 mS/cm"), ("100 mM", "11 cm")], "not a conductivity"),
    ],
)
def test_calibration_validation(points, message):
    with pytest.raises(ValidationError, match=message):
        ConductivityCalibration(points=points)


def test_conductivity_to_salt():
    salt = conductivity_to_salt(signal([0.6, 2.6], unit="S/m"), CALIBRATION)
    assert salt.kind == "salt" and salt.unit == "mol/m**3"
    assert salt.values == pytest.approx([50.0, 300.0])
    with pytest.raises(ChromSimError, match="expected conductivity"):
        conductivity_to_salt(signal([1.0]), CALIBRATION)


def test_constant_and_linear_baseline():
    x = np.linspace(0.0, 10.0, 101)
    drift = 0.5 + 0.1 * x
    peak = np.exp(-((x - 5.0) ** 2))
    raw = signal(drift + peak, x=x)
    constant = subtract_baseline(signal(0.5 + peak, x=x), (0.0, 1.0))
    assert constant.values[-1] == pytest.approx(0.0, abs=1e-6)
    linear = subtract_baseline(raw, (0.0, 1.0), (9.0, 10.0))
    assert linear.values == pytest.approx(peak, abs=1e-6)
    assert raw.values[0] == pytest.approx(0.5)  # input unchanged


def test_baseline_window_without_points():
    with pytest.raises(ChromSimError, match="contains no points"):
        subtract_baseline(signal([1.0, 2.0]), (5.0, 6.0))


def test_uv_to_concentration_lambert_beer():
    # 0.28 AU / (140 m²/kg × 0.002 m) = 1.0 kg/m³ (= 1 g/L); 1.4 L/(g·cm) = 140 m²/kg
    concentration = uv_to_concentration(signal([0.28]), 140.0, "0.2 cm")
    assert concentration.values == pytest.approx([1.0])
    assert concentration.unit == "kg/m**3" and concentration.kind == "concentration"


def test_resample():
    resampled = resample(signal([0.0, 2.0, 4.0]), [0.5, 1.5])
    assert resampled.values == pytest.approx([1.0, 3.0])


def test_volume_to_time_inverts_pumped_volume():
    steps = make_project().resolved_method()
    times = np.array([0.0, 1000.0, 2500.0, 9000.0])
    assert volume_to_time(pumped_volume(times, steps), steps) == pytest.approx(times)


def test_load_conductivity_calibration(tmp_path):
    from pathlib import Path

    example = Path(__file__).parents[2] / "examples/calibrations/nacl_conductivity.yaml"
    calibration = load_conductivity_calibration(example)
    assert calibration.salt_from_conductivity(4.55) == pytest.approx(500.0)  # 45.5 mS/cm
    with pytest.raises(ChromSimError, match="not found"):
        load_conductivity_calibration(tmp_path / "none.yaml")
    bad = tmp_path / "bad.yaml"
    bad.write_text("points: [[0 mM, 1 mS/cm]]\n")
    with pytest.raises(ChromSimError, match="Invalid conductivity calibration"):
        load_conductivity_calibration(bad)
