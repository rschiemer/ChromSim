"""Peak analysis: moments, plate number, HETP and asymmetry (SPEC §7.3).

All functions take x (volume in m³ or time in s) and a baseline-free signal y.

Equations
---------
Moments of a peak y(x):

.. math:: \\mu_0 = \\int y\\,dx, \\qquad \\mu_1 = \\frac{1}{\\mu_0}\\int x\\,y\\,dx, \\qquad
          \\sigma^2 = \\frac{1}{\\mu_0}\\int (x-\\mu_1)^2 y\\,dx

A rectangular injection of width Δ adds Δ/2 to μ_1 and Δ²/12 to σ². Plate number and plate
height: N = μ_1²/σ², HETP = L/N.

Asymmetry factor at 10 % of the peak height: A_s = b/a, with a the width from the front to the
apex and b from the apex to the tail.

References: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
Chromatography*, 2nd ed., Academic Press 2006, ch. 6; Schmidt-Traub, Schulte &
Seidel-Morgenstern (eds.), *Preparative Chromatography*, 2nd ed., Wiley-VCH 2012, ch. 6.
"""

from dataclasses import dataclass

import numpy as np

from chromsim.errors import ChromSimError


@dataclass(frozen=True)
class PeakMoments:
    """Moments of one peak.

    Attributes
    ----------
    area : float
        Zeroth moment μ_0, in units of y × x.
    mean : float
        First moment μ_1 (centre of mass), in x units.
    variance : float
        Second central moment σ², in x units squared.
    """

    area: float
    mean: float
    variance: float

    def corrected_for_injection(self, injection_width: float) -> "PeakMoments":
        """Remove a rectangular injection of width Δ: μ_1 − Δ/2, σ² − Δ²/12."""
        return PeakMoments(
            area=self.area,
            mean=self.mean - injection_width / 2.0,
            variance=self.variance - injection_width**2 / 12.0,
        )


def _window(x, y, window):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if window is None:
        return x, y
    inside = (x >= window[0]) & (x <= window[1])
    if inside.sum() < 3:
        raise ChromSimError(f"Window {window} contains fewer than 3 points.")
    return x[inside], y[inside]


def peak_moments(x, y, window: tuple[float, float] | None = None) -> PeakMoments:
    """Area, first moment and variance of a baseline-free peak.

    Parameters
    ----------
    x : array_like
        Axis (m³ or s).
    y : array_like
        Signal with the baseline removed.
    window : (float, float), optional
        Only use x within this range; noise far from the peak biases the variance.

    Raises
    ------
    ChromSimError
        If the peak area is not positive.
    """
    x, y = _window(x, y, window)
    area = float(np.trapezoid(y, x))
    if area <= 0.0:
        raise ChromSimError("Peak area is not positive; check the baseline and the window.")
    mean = float(np.trapezoid(x * y, x) / area)
    variance = float(np.trapezoid((x - mean) ** 2 * y, x) / area)
    return PeakMoments(area=area, mean=mean, variance=variance)


def plate_number(moments: PeakMoments) -> float:
    """N = μ_1² / σ² (dimensionless)."""
    return moments.mean**2 / moments.variance


def hetp(moments: PeakMoments, column_length_m: float) -> float:
    """Height equivalent to a theoretical plate, HETP = L σ²/μ_1², m."""
    return column_length_m / plate_number(moments)


def peak_maximum(x, y) -> tuple[float, float]:
    """Position and height of the maximum, refined with a parabola through the top three points.

    The parabola y = a x² + b x + c through the highest sample and its neighbours has its vertex
    at x_top = −b / (2a); this removes most of the error of the sampling grid.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    i = int(np.argmax(y))
    if i == 0 or i == len(y) - 1:
        return float(x[i]), float(y[i])
    a, b, c = np.polyfit(x[i - 1 : i + 2], y[i - 1 : i + 2], 2)
    if a >= 0:  # flat or not a maximum: keep the sample
        return float(x[i]), float(y[i])
    x_top = -b / (2 * a)
    if not x[i - 1] <= x_top <= x[i + 1]:
        return float(x[i]), float(y[i])
    return float(x_top), float(a * x_top**2 + b * x_top + c)


def asymmetry_factor(x, y, height_fraction: float = 0.1) -> float:
    """Asymmetry factor A_s = b/a at ``height_fraction`` of the peak height (default 10 %).

    Raises
    ------
    ChromSimError
        If the signal does not fall below that height on both sides of the peak.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    apex = int(np.argmax(y))
    level = height_fraction * y[apex]
    front = np.nonzero(y[:apex] < level)[0]
    tail = np.nonzero(y[apex:] < level)[0]
    if front.size == 0 or tail.size == 0:
        raise ChromSimError(
            f"The peak does not drop below {height_fraction:.0%} of its height on both sides; "
            "widen the data range."
        )
    i_front = front[-1]  # last point below the level before the apex
    i_tail = apex + tail[0]  # first point below the level after the apex
    x_front = np.interp(level, [y[i_front], y[i_front + 1]], [x[i_front], x[i_front + 1]])
    x_tail = np.interp(level, [y[i_tail], y[i_tail - 1]], [x[i_tail], x[i_tail - 1]])
    return float((x_tail - x[apex]) / (x[apex] - x_front))
