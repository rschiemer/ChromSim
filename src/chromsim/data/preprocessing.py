"""Preprocessing of experimental signals (SPEC §7.2).

All functions return new :class:`~chromsim.data.experiment.Signal` objects and never change
their input.

Equations
---------
UV → concentration (Lambert–Beer), with absorbance A (AU), mass extinction coefficient
ε (m²/kg) and path length l (m):

.. math:: c = \\frac{A}{\\varepsilon \\, l} \\quad [\\mathrm{kg/m^3}]

Reference: e.g. Harris, *Quantitative Chemical Analysis*, 9th ed., Freeman 2015, ch. 18.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, ValidationError, field_validator

from chromsim import units
from chromsim.core.method import ResolvedStep
from chromsim.data.experiment import Signal
from chromsim.errors import ChromSimError
from chromsim.units import Concentration, Quantity


class ConductivityCalibration(BaseModel):
    """Salt concentration ↔ conductivity, measured for one buffer system.

    Conductivity of salt solutions is not linear in concentration, so the calibration is a table
    of points, interpolated linearly between points and extrapolated from the end segments.

    Attributes
    ----------
    points : list of (Quantity, Quantity)
        (salt concentration, conductivity) pairs, e.g. ``[["0 mM", "0.1 mS/cm"],
        ["500 mM", "45 mS/cm"]]``; at least two, both increasing.

    Examples
    --------
    YAML::

        points:
          - [0 mM, 0.1 mS/cm]
          - [100 mM, 10.5 mS/cm]
          - [500 mM, 45 mS/cm]
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    points: list[tuple[Concentration, Quantity]]

    @field_validator("points")
    @classmethod
    def _increasing(cls, points):
        if len(points) < 2:
            raise ValueError("a conductivity calibration needs at least two points.")
        salt = [units.concentration_to_si(c) for c, _ in points]
        conductivity = [units.to_si(k, "conductivity") for _, k in points]
        if np.any(np.diff(salt) <= 0) or np.any(np.diff(conductivity) <= 0):
            raise ValueError("calibration points must increase in both salt and conductivity.")
        return points

    @property
    def salt_mol_per_m3(self) -> np.ndarray:
        """Calibration salt concentrations, mol/m³."""
        return np.array([units.concentration_to_si(c) for c, _ in self.points])

    @property
    def conductivity_s_per_m(self) -> np.ndarray:
        """Calibration conductivities, S/m."""
        return np.array([units.to_si(k, "conductivity") for _, k in self.points])

    def salt_from_conductivity(self, conductivity_s_per_m):
        """Salt concentration (mol/m³) for conductivities (S/m)."""
        return _interpolate(conductivity_s_per_m, self.conductivity_s_per_m, self.salt_mol_per_m3)

    def conductivity_from_salt(self, salt_mol_per_m3):
        """Conductivity (S/m) for salt concentrations (mol/m³)."""
        return _interpolate(salt_mol_per_m3, self.salt_mol_per_m3, self.conductivity_s_per_m)


def load_conductivity_calibration(path: str | Path) -> ConductivityCalibration:
    """Read a conductivity calibration (salt ↔ conductivity points) from a YAML file.

    Raises
    ------
    ChromSimError
        If the file is missing or the calibration is invalid.
    """
    path = Path(path)
    try:
        return ConductivityCalibration.model_validate(
            yaml.safe_load(path.read_text(encoding="utf-8"))
        )
    except FileNotFoundError:
        raise ChromSimError(f"Calibration file not found: {path}") from None
    except (yaml.YAMLError, ValidationError) as error:
        raise ChromSimError(f"Invalid conductivity calibration {path}: {error}") from error


def _interpolate(x, x_points, y_points):
    """Piecewise-linear interpolation with linear extrapolation from the end segments."""
    x = np.asarray(x, dtype=float)
    y = np.interp(x, x_points, y_points)
    low_slope = (y_points[1] - y_points[0]) / (x_points[1] - x_points[0])
    high_slope = (y_points[-1] - y_points[-2]) / (x_points[-1] - x_points[-2])
    y = np.where(x < x_points[0], y_points[0] + low_slope * (x - x_points[0]), y)
    return np.where(x > x_points[-1], y_points[-1] + high_slope * (x - x_points[-1]), y)


def subtract_baseline(signal: Signal, *windows: tuple[float, float]) -> Signal:
    """Remove a baseline estimated in signal-free windows.

    Parameters
    ----------
    signal : Signal
        Signal to correct.
    *windows : (float, float)
        One window (x_start, x_end) in the signal's x units (m³ or s): subtract its mean
        (constant baseline). Two windows: subtract the straight line through the two window
        means (linear drift).

    Returns
    -------
    Signal
        Corrected signal.
    """
    if len(windows) not in (1, 2):
        raise ChromSimError("subtract_baseline needs one window (constant) or two (linear drift).")
    centres, means = [], []
    for start, end in windows:
        inside = (signal.x >= start) & (signal.x <= end)
        if not np.any(inside):
            raise ChromSimError(
                f"Baseline window {start:g}–{end:g} contains no points of signal '{signal.name}' "
                f"(its axis runs from {signal.x.min():g} to {signal.x.max():g} in SI units)."
            )
        centres.append(0.5 * (start + end))
        means.append(float(np.mean(signal.values[inside])))
    if len(windows) == 1:
        baseline = np.full_like(signal.values, means[0])
    else:
        slope = (means[1] - means[0]) / (centres[1] - centres[0])
        baseline = means[0] + slope * (signal.x - centres[0])
    return replace(signal, values=signal.values - baseline)


def uv_to_concentration(
    signal: Signal, extinction_coeff_m2_per_kg: float, path_length: Quantity | str = "0.2 cm"
) -> Signal:
    """Convert a UV absorbance signal to mass concentration with Lambert–Beer.

    Parameters
    ----------
    signal : Signal
        Absorbance in AU (baseline-corrected).
    extinction_coeff_m2_per_kg : float
        Mass extinction coefficient ε, m²/kg (``Component.extinction_coeff_m2_per_kg``).
    path_length : Quantity or str
        Flow-cell path length l, default 0.2 cm (2 mm, ÄKTA UV monitors).

    Returns
    -------
    Signal
        Concentration in kg/m³ (= g/L), kind ``"concentration"``.
    """
    if signal.unit != "AU":
        raise ChromSimError(
            f"Signal '{signal.name}' is in {signal.unit}, expected absorbance (AU)."
        )
    length_m = units.to_si(Quantity.parse(path_length), "length")
    epsilon, length = extinction_coeff_m2_per_kg, length_m
    return replace(
        signal,
        kind="concentration",
        values=signal.values / (epsilon * length),
        unit="kg/m**3",
        original_unit="g/L",
    )


def conductivity_to_salt(signal: Signal, calibration: ConductivityCalibration) -> Signal:
    """Convert a conductivity signal (S/m) to salt concentration (mol/m³) via a calibration."""
    if signal.unit != "S/m":
        raise ChromSimError(
            f"Signal '{signal.name}' is in {signal.unit}, expected conductivity (S/m)."
        )
    return replace(
        signal,
        kind="salt",
        values=calibration.salt_from_conductivity(signal.values),
        unit="mol/m**3",
        original_unit="mM",
    )


def resample(signal: Signal, x_new) -> Signal:
    """Linearly interpolate a signal onto new x values (same axis type and units)."""
    x_new = np.asarray(x_new, dtype=float)
    return replace(signal, x=x_new, values=np.interp(x_new, signal.x, signal.values))


def volume_to_time(volume_m3, steps: list[ResolvedStep]):
    """Time (s) at which a pumped volume (m³) is reached, for the flow rates of a method.

    Inverse of :func:`chromsim.backends.base.pumped_volume`.
    """
    knot_times = [steps[0].start_s] + [step.end_s for step in steps]
    knot_volumes = [0.0]
    for step in steps:
        knot_volumes.append(knot_volumes[-1] + step.flow_rate_m3_per_s * step.duration_s)
    return np.interp(volume_m3, knot_volumes, knot_times)
