"""Yamamoto regression on exact synthetic data and its helpers."""

import numpy as np
import pytest
from project_examples import make_project

from chromsim.data.experiment import Signal
from chromsim.errors import ChromSimError
from chromsim.estimation.workflows.yamamoto import (
    LgeRun,
    elution_salt_at_peak,
    gradient_slope,
    yamamoto_fit,
)

V_C, EPS_T, LAMBDA = 3.927e-6, 0.904, 1200.0


def exact_runs(nu, k_eq, slopes):
    """I_R from GH = I_R^(ν+1) / ((ν+1) k_eq Λ^ν) for each gradient slope."""
    runs = []
    for g in slopes:
        gh = g * (1 - EPS_T) * V_C
        elution_salt = (gh * (nu + 1) * k_eq * LAMBDA**nu) ** (1 / (nu + 1))
        runs.append(LgeRun(f"g={g:g}", g, elution_salt))
    return runs


def test_exact_data_are_recovered():
    slopes = [1e7, 2e7, 4e7, 8e7]  # mol/m³ per m³ (10–80 mM per mL)
    result = yamamoto_fit(exact_runs(6.0, 1e-2, slopes), V_C, EPS_T, LAMBDA)
    assert result.nu == pytest.approx(6.0, rel=1e-9)
    assert result.k_eq == pytest.approx(1e-2, rel=1e-6)
    assert result.r_squared == pytest.approx(1.0)
    assert np.allclose(result.table["residual"], 0.0, atol=1e-10)


def test_needs_two_different_slopes_and_positive_nu():
    with pytest.raises(ChromSimError, match="two different gradient slopes"):
        yamamoto_fit([LgeRun("a", 1e7, 300.0), LgeRun("b", 1e7, 310.0)], V_C, EPS_T, LAMBDA)
    with pytest.raises(ChromSimError, match="ν = .* ≤ 0"):
        yamamoto_fit([LgeRun("a", 1e7, 300.0), LgeRun("b", 2e7, 200.0)], V_C, EPS_T, LAMBDA)


def test_gradient_slope_of_example_project():
    elute = make_project().resolved_method()[3]
    # 450 mM over 20 CV × 3.92699 mL = 78.54 mL → 5.7296 mM/mL = 5.7296e6 mol/m³ per m³
    assert gradient_slope(elute) == pytest.approx(450.0 / (20 * 3.92699e-6), rel=1e-5)
    with pytest.raises(ChromSimError, match="not an increasing salt gradient"):
        gradient_slope(make_project().resolved_method()[2])


def test_elution_salt_corrects_detector_delays():
    x = np.linspace(0.0, 100e-6, 10001)
    uv = Signal("uv", "uv", "volume", x, np.exp(-0.5 * ((x - 50.5e-6) / 2e-6) ** 2), "AU", "mAU")
    salt = Signal("salt", "salt", "volume", x, 3e6 * x, "mol/m**3", "mM")  # 3 mM per mL
    # UV sees the outlet 0.5 mL late, conductivity 1.0 mL late: the outlet peak is at 50 mL, where
    # the conductivity detector reads 51 mL × 3 mM/mL = 153 mM.
    assert elution_salt_at_peak(
        uv, salt, protein_delay_m3=0.5e-6, salt_delay_m3=1.0e-6
    ) == pytest.approx(153.0, rel=1e-4)
    with pytest.raises(ChromSimError, match="conductivity_to_salt"):
        elution_salt_at_peak(uv, uv)
