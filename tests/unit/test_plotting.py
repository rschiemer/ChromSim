"""Plot helper builds the expected traces from a result (no simulation needed)."""

import numpy as np
import pytest

from chromsim.backends.base import SimulationResult
from chromsim.core.method import ResolvedStep
from chromsim.plotting import plot_chromatogram


def fake_result() -> SimulationResult:
    time_s = np.linspace(0.0, 100.0, 11)
    step = ResolvedStep("run", 0.0, 100.0, 1e-8, (50.0, 0.0), (50.0, 0.0))
    return SimulationResult(
        component_names=["salt", "mAb"],
        time_s=time_s,
        outlet_mol_per_m3=np.column_stack([np.full(11, 50.0), np.linspace(0.0, 1e-3, 11)]),
        volume_m3=time_s * 1e-8,
        column_volume_m3=1e-6,
        steps=[step],
        molar_masses_kg_per_mol=[None, 148.0],
    )


def test_proteins_left_axis_in_g_per_l_and_salt_right_axis_in_mm():
    figure = plot_chromatogram(fake_result())
    traces = {trace.name: trace for trace in figure.data}
    assert set(traces) == {"salt / mM", "mAb / g/L"}
    assert traces["mAb / g/L"].y[-1] == pytest.approx(1e-3 * 148.0)
    assert traces["salt / mM"].yaxis == "y2"
    assert traces["mAb / g/L"].x[-1] == pytest.approx(1.0)  # 1e-6 m³ pumped = 1 CV


def test_time_axis_in_minutes():
    figure = plot_chromatogram(fake_result(), x="min")
    assert figure.data[0].x[-1] == pytest.approx(100.0 / 60.0)


def test_unknown_axis_is_rejected():
    with pytest.raises(ValueError, match="CV, mL, min"):
        plot_chromatogram(fake_result(), x="hours")
