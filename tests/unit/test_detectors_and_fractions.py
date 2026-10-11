"""Detector predictions and fraction averaging on hand-made results (exact expectations)."""

import numpy as np
import pandas as pd
import pytest
from project_examples import make_project

from chromsim.analysis.detectors import predict_conductivity, predict_uv
from chromsim.analysis.fractionation import compare_with_assay, fraction_averages, interval_average
from chromsim.backends.base import SimulationResult
from chromsim.data.assays import FractionAssay
from chromsim.data.experiment import Fraction, fractions_from_marks
from chromsim.data.preprocessing import ConductivityCalibration
from chromsim.errors import ChromSimError

ML = 1e-6


def linear_result() -> SimulationResult:
    """salt = 50 mM constant; mAb rises linearly with volume: c = 1e-3 mol/m³ per mL."""
    volume = np.linspace(0.0, 10.0, 101) * ML
    mab = 1e-3 * volume / ML
    return SimulationResult(
        component_names=["salt", "mAb"],
        time_s=volume / 1e-8,
        outlet_mol_per_m3=np.column_stack([np.full(volume.size, 50.0), mab]),
        volume_m3=volume,
        column_volume_m3=4 * ML,
        steps=[],
        molar_masses_kg_per_mol=[None, 148.0],
    )


def test_interval_average_is_exact_for_linear_profiles():
    x = np.array([0.0, 1.0, 2.0])
    assert interval_average(x, 3.0 * x, 0.25, 1.75) == pytest.approx(
        3.0
    )  # mean of 3x over [0.25,1.75]
    with pytest.raises(ChromSimError, match="must be after"):
        interval_average(x, x, 1.0, 1.0)


def test_fraction_averages_and_delay():
    result = linear_result()
    fractions = [Fraction("F1", 2 * ML, 4 * ML)]
    table = fraction_averages(result, fractions).set_index("component")
    assert table.loc["mAb", "concentration_mol_per_m3"] == pytest.approx(3e-3)  # mid-point 3 mL
    assert table.loc["salt", "concentration_mol_per_m3"] == pytest.approx(50.0)
    delayed = fraction_averages(result, fractions, delay_volume_m3=1 * ML).set_index("component")
    assert delayed.loc["mAb", "concentration_mol_per_m3"] == pytest.approx(2e-3)  # sees 1–3 mL
    with pytest.raises(ChromSimError, match="after the simulation"):
        fraction_averages(result, [Fraction("late", 9 * ML, 12 * ML)])


def test_fractions_from_marks():
    marks = [(5 * ML, "A2"), (0.0, "Waste"), (3 * ML, "A1"), (7 * ML, "Waste")]
    assert fractions_from_marks(marks, end_m3=10 * ML) == [
        Fraction("A1", 3 * ML, 5 * ML),
        Fraction("A2", 5 * ML, 7 * ML),
    ]


def assay(volumes=True):
    results = pd.DataFrame(
        {
            "fraction": ["F1", "F2", "X9"],
            "analyte": ["protein", "protein", "protein"],
            "value": [0.45, np.nan, 1.0],
            "unit": ["g/L", "g/L", "g/L"],
            "sd": [0.01, 0.01, 0.01],
            "below_loq": [False, True, False],
            "loq": [None, None, None],
        }
    )
    fraction_volumes = {"F1": (2 * ML, 4 * ML), "F2": (4 * ML, 6 * ML)} if volumes else {}
    return FractionAssay("a", results, fraction_volumes, readouts={"protein": ["mAb"]})


def test_compare_with_assay_converts_units_and_reports_mismatches():
    comparison = compare_with_assay(linear_result(), assay())
    table = comparison.table.set_index("fraction")
    # F1: 3e-3 mol/m³ × 148 kg/mol = 0.444 g/L
    assert table.loc["F1", "simulated"] == pytest.approx(0.444)
    assert table.loc["F1", "measured"] == pytest.approx(0.45)
    assert table.loc["F2", "below_loq"]
    assert comparison.assay_fractions_not_collected == ["X9"]
    assert comparison.collected_fractions_without_assay == []


def test_compare_uses_collected_fractions_when_assay_has_no_volumes():
    collected = [Fraction("F1", 2 * ML, 4 * ML), Fraction("F3", 6 * ML, 8 * ML)]
    comparison = compare_with_assay(linear_result(), assay(volumes=False), fractions=collected)
    assert list(comparison.table["fraction"]) == ["F1"]
    assert comparison.assay_fractions_not_collected == ["F2", "X9"]
    assert comparison.collected_fractions_without_assay == ["F3"]
    with pytest.raises(ChromSimError, match="pass fractions="):
        compare_with_assay(linear_result(), assay(volumes=False))


def test_readout_with_unknown_component():
    bad = assay()
    bad.readouts = {"protein": ["IgG"]}
    with pytest.raises(ChromSimError, match="unknown component.*IgG"):
        compare_with_assay(linear_result(), bad)


def test_predict_uv_lambert_beer_and_delay():
    project = make_project()  # mAb: 1.4 L/(g·cm) = 140 m²/kg, 148 kDa; salt no extinction
    project.system.uv_delay_volume = "0.5 mL"
    uv = predict_uv(linear_result(), project, path_length="0.2 cm")
    # A = l ε M c = 0.002 m × 140 m²/kg × 148 kg/mol × c
    assert uv.values == pytest.approx(0.002 * 140.0 * 148.0 * linear_result().outlet("mAb"))
    assert uv.x == pytest.approx(linear_result().volume_m3 + 0.5 * ML)
    assert uv.unit == "AU"


def test_predict_conductivity():
    project = make_project()
    calibration = ConductivityCalibration(points=[("0 mM", "0 mS/cm"), ("100 mM", "10 mS/cm")])
    cond = predict_conductivity(linear_result(), project, calibration)
    assert cond.values == pytest.approx(0.5)  # 50 mM → 5 mS/cm = 0.5 S/m
    with pytest.raises(ChromSimError, match="salt_component"):
        predict_conductivity(linear_result(), project, calibration, salt_component="NaCl")
