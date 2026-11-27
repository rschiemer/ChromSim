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


def test_overlay_draws_measured_solid_and_simulated_dashed_with_fractions():
    from chromsim.data.experiment import Experiment, Fraction, Signal
    from chromsim.plotting import plot_overlay

    x = np.linspace(0.0, 10e-6, 11)
    uv = Signal("UV", "uv", "volume", x, np.linspace(0, 0.1, 11), "AU", "mAU")
    cond = Signal("Cond", "conductivity", "volume", x, np.full(11, 1.0), "S/m", "mS/cm")
    experiment = Experiment("e", {"UV": uv, "Cond": cond}, fractions=[Fraction("A1", 2e-6, 4e-6)])
    simulated = Signal("UV (sim)", "uv", "volume", x, np.linspace(0, 0.1, 11), "AU", "mAU")
    figure = plot_overlay(experiment, simulated=[simulated])
    traces = {trace.name: trace for trace in figure.data}
    assert traces["UV / mAU"].y[-1] == pytest.approx(100.0)  # AU → mAU
    assert traces["Cond / mS/cm"].y[0] == pytest.approx(10.0)  # S/m → mS/cm
    assert traces["UV (sim) / mAU"].line.dash == "dash"
    assert traces["UV / mAU"].x[-1] == pytest.approx(10.0)  # mL
    assert len(figure.layout.shapes) == 1  # one shaded fraction
    with pytest.raises(ValueError, match="needs column_volume_m3"):
        plot_overlay(experiment, x="CV")


def test_overlay_axis_titles_list_only_drawn_signals():
    from chromsim.data.experiment import Signal
    from chromsim.plotting import plot_overlay

    x = np.linspace(0.0, 10e-6, 11)
    uv = Signal("UV (sim)", "uv", "volume", x, np.zeros(11), "AU", "mAU")
    cond = Signal("Cond (sim)", "conductivity", "volume", x, np.ones(11), "S/m", "mS/cm")
    ph = Signal("pH (sim)", "ph", "volume", x, np.full(11, 7.0), "pH", "pH")
    figure = plot_overlay(simulated=[uv, cond, ph], kinds=("conductivity", "ph"))
    assert [trace.name for trace in figure.data] == ["Cond (sim) / mS/cm", "pH (sim) / pH"]
    assert figure.layout.yaxis.title.text == "UV"  # no UV curve drawn, so no UV unit
    assert figure.layout.yaxis2.title.text == "mS/cm, pH"


def test_fraction_comparison_plot():
    import pandas as pd

    from chromsim.analysis.fractionation import AssayComparison
    from chromsim.plotting import plot_fraction_comparison

    table = pd.DataFrame(
        {
            "fraction": ["A1", "A2"],
            "analyte": ["monomer", "monomer"],
            "unit": ["g/L", "g/L"],
            "measured": [0.5, np.nan],
            "sd": [0.02, 0.01],
            "below_loq": [False, True],
            "simulated": [0.45, 0.01],
            "start_m3": [1e-6, 2e-6],
            "end_m3": [2e-6, 3e-6],
        }
    )
    figure = plot_fraction_comparison(AssayComparison(table), "monomer")
    names = [trace.name for trace in figure.data]
    assert names == ["monomer measured", "below LOQ", "monomer simulated (fraction average)"]
    assert list(figure.data[0].x) == pytest.approx([1.5])  # fraction middle in mL
    with pytest.raises(ValueError, match="No comparison rows"):
        plot_fraction_comparison(AssayComparison(table), "hcp")


def test_yamamoto_and_correlation_plots():
    import pandas as pd

    from chromsim.estimation.uncertainty import UncertaintyReport
    from chromsim.estimation.workflows.yamamoto import LgeRun, yamamoto_fit
    from chromsim.plotting import plot_correlation, plot_yamamoto

    runs = [LgeRun("a", 1e7, 300.0), LgeRun("b", 2e7, 330.0), LgeRun("c", 4e7, 365.0)]
    figure = plot_yamamoto(yamamoto_fit(runs, 3.9e-6, 0.9, 1200.0))
    assert [trace.mode for trace in figure.data] == ["markers", "lines", "markers"]
    report = UncertaintyReport(
        table=pd.DataFrame(),
        correlation=pd.DataFrame(
            [[1.0, -0.99], [-0.99, 1.0]], index=["a", "b"], columns=["a", "b"]
        ),
    )
    heatmap = plot_correlation(report).data[0]
    assert list(heatmap.x) == ["a", "b"] and heatmap.z[0][1] == pytest.approx(-0.99)


def test_chromatogram_pools_and_pareto_plot():
    import pandas as pd

    from chromsim.analysis.fractionation import Pool
    from chromsim.plotting import plot_pareto

    figure = plot_chromatogram(fake_result(), x="min", pools=[Pool("mAb", 30.0, 60.0)])
    assert len(figure.layout.shapes) >= 2  # step line + pool shading
    shading = [shape for shape in figure.layout.shapes if shape.type == "rect"]
    assert shading[0].x0 == pytest.approx(0.5) and shading[0].x1 == pytest.approx(1.0)  # minutes
    table = pd.DataFrame({"recovery": [0.9, 0.8], "productivity": [1.0, 2.0], "gradient": [30, 10]})
    pareto = plot_pareto(table, "recovery", "productivity", color="gradient")
    assert list(pareto.data[0].x) == [0.8, 0.9]  # sorted along x
