"""Plotly figures shared by notebooks and the UI (SPEC §4.1).

Plotting only reads results; it never computes science beyond unit conversion for display.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from chromsim import units
from chromsim.backends.base import SimulationResult

X_AXES = {
    "CV": ("Volume / CV", lambda result: result.column_volumes),
    "mL": ("Volume / mL", lambda result: units.from_si(1.0, "mL") * result.volume_m3),
    "min": ("Time / min", lambda result: result.time_s / 60.0),
}


def plot_chromatogram(
    result: SimulationResult, x: str = "CV", protein_unit: str = "g/L"
) -> go.Figure:
    """Outlet concentrations over volume or time, with method steps marked.

    Components with a molecular weight (proteins) are drawn on the left axis in
    ``protein_unit``; components without one (salt, pH, tracers) on the right axis in mM.

    Parameters
    ----------
    result : SimulationResult
        A simulation result.
    x : {"CV", "mL", "min"}, optional
        Horizontal axis: column volumes (default), millilitres or minutes.
    protein_unit : str, optional
        Concentration unit for components with a molecular weight, default ``"g/L"``.

    Returns
    -------
    plotly.graph_objects.Figure
        The chromatogram.
    """
    if x not in X_AXES:
        raise ValueError(f"x must be one of {', '.join(X_AXES)}, got '{x}'.")
    x_title, x_values_of = X_AXES[x]
    x_values = x_values_of(result)

    figure = make_subplots(specs=[[{"secondary_y": True}]])
    for name, molar_mass in zip(
        result.component_names, result.molar_masses_kg_per_mol, strict=True
    ):
        if molar_mass is not None:
            figure.add_trace(
                go.Scatter(
                    x=x_values, y=result.outlet(name, protein_unit), name=f"{name} / {protein_unit}"
                ),
                secondary_y=False,
            )
        else:
            figure.add_trace(
                go.Scatter(
                    x=x_values,
                    y=result.outlet(name, "mM"),
                    name=f"{name} / mM",
                    line={"dash": "dot"},
                ),
                secondary_y=True,
            )

    for step in result.steps:
        start = float(np.interp(step.start_s, result.time_s, x_values))
        figure.add_vline(x=start, line={"width": 1, "color": "lightgray"})
        figure.add_annotation(
            x=start,
            y=1.0,
            yref="paper",
            text=step.name,
            showarrow=False,
            textangle=-90,  # vertical, so labels of short steps do not overlap
            xanchor="left",
            yanchor="top",
            font={"size": 10, "color": "gray"},
        )

    figure.update_xaxes(title_text=x_title)
    figure.update_yaxes(title_text=f"Concentration / {protein_unit}", secondary_y=False)
    figure.update_yaxes(
        title_text="Concentration / mM",
        secondary_y=True,
        rangemode="tozero",
        tickmode="auto",
        showgrid=False,
    )
    figure.update_yaxes(rangemode="tozero", secondary_y=False)
    figure.update_layout(template="simple_white", legend={"orientation": "h", "y": -0.2})
    return figure


def _x_values(x_m3, x: str, column_volume_m3: float | None):
    if x == "mL":
        return x_m3 * 1e6
    if x == "CV":
        if column_volume_m3 is None:
            raise ValueError("x='CV' needs column_volume_m3 (e.g. project.column.volume_m3).")
        return x_m3 / column_volume_m3
    raise ValueError(f"x must be 'mL' or 'CV', got '{x}'.")


def plot_overlay(
    experiment=None,
    simulated=(),
    kinds=("uv", "conductivity"),
    x: str = "mL",
    column_volume_m3: float | None = None,
    show_fractions: bool = True,
) -> go.Figure:
    """Measured and simulated detector signals on one volume axis.

    UV is drawn on the left axis, all other kinds on the right axis, each in the unit of the
    original export (e.g. mAU, mS/cm). Measured curves are solid, simulated curves dashed.

    Parameters
    ----------
    experiment : Experiment, optional
        Imported run; its signals of the requested kinds and its fractions are drawn.
    simulated : sequence of Signal
        Predicted signals, e.g. from ``predict_uv`` and ``predict_conductivity``.
    kinds : sequence of str
        Signal kinds to draw.
    x : {"mL", "CV"}
        Horizontal axis.
    column_volume_m3 : float, optional
        Needed for ``x="CV"``.
    show_fractions : bool
        Mark collected fractions.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    figure = make_subplots(specs=[[{"secondary_y": True}]])
    measured = []
    if experiment is not None:
        measured = [
            s for s in experiment.signals.values() if s.kind in kinds and s.x_axis == "volume"
        ]
    for signal, dash in [(s, "solid") for s in measured] + [
        (s, "dash") for s in simulated if s.kind in kinds
    ]:
        figure.add_trace(
            go.Scatter(
                x=_x_values(signal.x, x, column_volume_m3),
                y=signal.values_in(signal.original_unit),
                name=f"{signal.name} / {signal.original_unit}",
                line={"dash": dash},
            ),
            secondary_y=signal.kind != "uv",
        )
    if experiment is not None and show_fractions and experiment.fractions:
        # Shade fractions alternately; label at most ~12 of them so labels stay readable.
        label_every = max(1, int(np.ceil(len(experiment.fractions) / 12)))
        for index, fraction in enumerate(experiment.fractions):
            start = float(_x_values(np.array(fraction.start_m3), x, column_volume_m3))
            end = float(_x_values(np.array(fraction.end_m3), x, column_volume_m3))
            figure.add_vrect(
                x0=start,
                x1=end,
                fillcolor="gray",
                opacity=0.06 if index % 2 else 0.12,
                line_width=0,
            )
            if index % label_every == 0:
                figure.add_annotation(
                    x=start,
                    y=1.0,
                    yref="paper",
                    text=fraction.label,
                    showarrow=False,
                    xanchor="left",
                    yanchor="bottom",
                    font={"size": 9, "color": "gray"},
                )
    left_units = sorted({s.original_unit for s in measured + list(simulated) if s.kind == "uv"})
    right_units = sorted(
        {s.original_unit for s in measured + list(simulated) if s.kind != "uv" and s.kind in kinds}
    )
    figure.update_xaxes(title_text="Volume / mL" if x == "mL" else "Volume / CV")
    figure.update_yaxes(
        title_text="UV / " + ", ".join(left_units) if left_units else "UV",
        secondary_y=False,
        rangemode="tozero",
    )
    figure.update_yaxes(
        title_text=" , ".join(right_units), secondary_y=True, showgrid=False, tickmode="auto"
    )
    figure.update_layout(template="simple_white", legend={"orientation": "h", "y": -0.2})
    return figure


def plot_fraction_comparison(
    comparison, analyte: str, result: SimulationResult | None = None, components=None
) -> go.Figure:
    """Measured vs simulated fraction-averaged values of one readout, as bars over the fractions.

    Parameters
    ----------
    comparison : AssayComparison
        From :func:`chromsim.analysis.fractionation.compare_with_assay`.
    analyte : str
        Readout to show.
    result : SimulationResult, optional
        Also draw the simulated outlet curve of the readout ...
    components : list of str, optional
        ... summed over these components (e.g. ``assay.readouts[analyte]``).

    Returns
    -------
    plotly.graph_objects.Figure
        Measured bars with SD error bars (values below LOQ as open markers at 0), simulated
        fraction averages as outlined bars, optional continuous simulated curve.
    """
    table = comparison.table[comparison.table["analyte"] == analyte]
    if table.empty:
        raise ValueError(f"No comparison rows for analyte '{analyte}'.")
    unit = table["unit"].iloc[0]
    middle_ml = 0.5 * (table["start_m3"] + table["end_m3"]) * 1e6
    width_ml = (table["end_m3"] - table["start_m3"]) * 1e6
    quantified = ~table["below_loq"].astype(bool)

    figure = go.Figure()
    figure.add_trace(
        go.Bar(
            x=middle_ml[quantified],
            y=table["measured"][quantified],
            width=width_ml[quantified],
            error_y={"type": "data", "array": table["sd"][quantified]},
            name=f"{analyte} measured",
            marker={"color": "rgba(31, 119, 180, 0.5)"},
            text=table["fraction"][quantified],
        )
    )
    if (~quantified).any():
        figure.add_trace(
            go.Scatter(
                x=middle_ml[~quantified],
                y=np.zeros((~quantified).sum()),
                mode="markers",
                marker={"symbol": "circle-open", "size": 9},
                name="below LOQ",
                text=table["fraction"][~quantified],
            )
        )
    figure.add_trace(
        go.Bar(
            x=middle_ml,
            y=table["simulated"],
            width=width_ml,
            name=f"{analyte} simulated (fraction average)",
            marker={"color": "rgba(0,0,0,0)", "line": {"color": "rgb(214, 39, 40)", "width": 2}},
        )
    )
    if result is not None and components:
        curve = sum(result.outlet(name, unit) for name in components)
        figure.add_trace(
            go.Scatter(
                x=result.volume_m3 * 1e6,
                y=curve,
                name=f"{analyte} simulated",
                line={"color": "rgb(214, 39, 40)"},
            )
        )
        figure.update_xaxes(
            range=[float(table["start_m3"].min() * 1e6) - 2, float(table["end_m3"].max() * 1e6) + 2]
        )
    figure.update_layout(
        barmode="overlay",
        template="simple_white",
        xaxis_title="Volume / mL",
        yaxis_title=f"{analyte} / {unit}",
        legend={"orientation": "h", "y": -0.2},
    )
    return figure


def plot_yamamoto(result) -> go.Figure:
    """Yamamoto plot: log10 GH against log10 I_R with the regression line, and residuals below.

    Parameters
    ----------
    result : YamamotoResult
        From :func:`chromsim.estimation.workflows.yamamoto.yamamoto_fit`.
    """
    table = result.table
    figure = make_subplots(
        rows=2, cols=1, shared_xaxes=True, row_heights=[0.75, 0.25], vertical_spacing=0.05
    )
    figure.add_trace(
        go.Scatter(
            x=table["log10_elution_salt"],
            y=table["log10_gh"],
            mode="markers",
            name="runs",
            text=table["name"],
            marker={"size": 9},
        ),
        row=1,
        col=1,
    )
    line_x = np.linspace(
        table["log10_elution_salt"].min() - 0.02, table["log10_elution_salt"].max() + 0.02, 50
    )
    figure.add_trace(
        go.Scatter(
            x=line_x,
            y=result.slope * line_x + result.intercept,
            mode="lines",
            name=f"fit: ν = {result.nu:.2f}, k_eq = {result.k_eq:.3g}, R² = {result.r_squared:.4f}",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=table["log10_elution_salt"], y=table["residual"], mode="markers", showlegend=False
        ),
        row=2,
        col=1,
    )
    figure.add_hline(y=0.0, line={"width": 1, "color": "gray"}, row=2, col=1)
    figure.update_yaxes(title_text="log10 GH / mM", row=1, col=1)
    figure.update_yaxes(title_text="residual", row=2, col=1)
    figure.update_xaxes(title_text="log10 I_R / mM", row=2, col=1)
    figure.update_layout(template="simple_white", legend={"orientation": "h", "y": -0.25})
    return figure


def plot_correlation(report) -> go.Figure:
    """Heat map of the parameter correlation matrix of an :class:`UncertaintyReport`."""
    correlation = report.correlation
    figure = go.Figure(
        go.Heatmap(
            z=correlation.to_numpy(),
            x=list(correlation.columns),
            y=list(correlation.index),
            zmin=-1,
            zmax=1,
            colorscale="RdBu",
            reversescale=True,
            text=np.round(correlation.to_numpy(), 2),
            texttemplate="%{text}",
        )
    )
    figure.update_layout(
        template="simple_white", title="Parameter correlation", yaxis={"autorange": "reversed"}
    )
    return figure
