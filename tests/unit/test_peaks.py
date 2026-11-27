"""Peak moments, plate number, HETP, asymmetry and peak maximum on analytic peaks."""

import numpy as np
import pytest

from chromsim.analysis.peaks import (
    PeakMoments,
    asymmetry_factor,
    hetp,
    peak_maximum,
    peak_moments,
    plate_number,
)
from chromsim.errors import ChromSimError

X = np.linspace(0.0, 100.0, 20001)


def gaussian(mean, sigma, area=1.0):
    return area * np.exp(-0.5 * ((X - mean) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def test_moments_of_a_gaussian():
    moments = peak_moments(X, gaussian(40.0, 3.0, area=2.0))
    assert moments.area == pytest.approx(2.0, rel=1e-6)
    assert moments.mean == pytest.approx(40.0, rel=1e-9)
    assert moments.variance == pytest.approx(9.0, rel=1e-6)


def test_injection_correction_of_a_rectangular_pulse():
    # A rectangular pulse from 10 to 16 has mean 13 and variance 6²/12 = 3: removing a width-6
    # injection leaves a zero-width peak at 10.
    y = ((X >= 10.0) & (X <= 16.0)).astype(float)
    corrected = peak_moments(X, y).corrected_for_injection(6.0)
    assert corrected.mean == pytest.approx(10.0, abs=1e-3)
    assert corrected.variance == pytest.approx(0.0, abs=1e-2)


def test_plate_number_and_hetp():
    moments = PeakMoments(area=1.0, mean=100.0, variance=25.0)
    assert plate_number(moments) == pytest.approx(400.0)
    assert hetp(moments, column_length_m=0.2) == pytest.approx(0.2 / 400.0)


def test_asymmetry_of_symmetric_and_tailing_peaks():
    assert asymmetry_factor(X, gaussian(50.0, 2.0)) == pytest.approx(1.0, rel=1e-3)
    # exponentially modified Gaussian tails: As > 1
    tailing = np.convolve(gaussian(30.0, 1.0), np.exp(-X / 3.0), mode="full")[: X.size]
    assert asymmetry_factor(X, tailing) > 1.5


def test_asymmetry_needs_both_sides():
    with pytest.raises(ChromSimError, match="both sides"):
        asymmetry_factor(X, gaussian(1.0, 5.0))


def test_peak_maximum_is_refined_between_samples():
    x = np.arange(0.0, 10.0, 1.0)
    y = -((x - 4.3) ** 2)
    position, height = peak_maximum(x, y)
    assert position == pytest.approx(4.3)
    assert height == pytest.approx(0.0, abs=1e-12)


def test_window_and_non_positive_area():
    assert peak_moments(X, gaussian(40.0, 3.0), window=(20.0, 60.0)).mean == pytest.approx(40.0)
    with pytest.raises(ChromSimError, match="not positive"):
        peak_moments(X, -gaussian(40.0, 3.0))


def test_asymmetry_window_ignores_a_larger_peak_outside_it():
    tailing = np.convolve(gaussian(30.0, 1.0), np.exp(-X / 3.0), mode="full")[: X.size]
    y = gaussian(85.0, 2.0) + tailing  # the tailing peak is the global maximum
    assert asymmetry_factor(X, y) > 1.5
    assert asymmetry_factor(X, y, window=(70.0, 100.0)) == pytest.approx(1.0, rel=1e-3)
