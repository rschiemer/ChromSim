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
