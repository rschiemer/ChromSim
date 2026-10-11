"""Fraction-averaged simulation results and comparison with fraction assays (SPEC §7.2b).

The fraction-averaged concentration of component i in a fraction collected between pumped
volumes V_a and V_b is

.. math:: \\bar c_i = \\frac{1}{V_b - V_a} \\int_{V_a}^{V_b} c_i(V - V_{delay}) \\, dV

where V_delay is the volume between the simulated outlet and the fraction collector. Cut points,
pooling and KPIs (yield, purity, productivity) follow in Phase 6.
"""

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from chromsim import units
from chromsim.backends.base import SimulationResult
from chromsim.data.assays import FractionAssay
from chromsim.data.experiment import Fraction
from chromsim.errors import ChromSimError

logger = logging.getLogger(__name__)


def interval_average(x, y, start: float, end: float) -> float:
    """Average of the piecewise-linear function y(x) over [start, end] (exact for linear pieces)."""
    if end <= start:
        raise ChromSimError(f"Fraction end {end:g} must be after its start {start:g}.")
    inside = (x > start) & (x < end)
    x_part = np.concatenate([[start], x[inside], [end]])
    y_part = np.interp(x_part, x, y)
    return float(np.trapezoid(y_part, x_part) / (end - start))


def fraction_averages(
    result: SimulationResult, fractions: list[Fraction], delay_volume_m3: float = 0.0
) -> pd.DataFrame:
    """Fraction-averaged outlet concentrations of every component.

    Parameters
    ----------
    result : SimulationResult
        The simulation.
    fractions : list of Fraction
        Fractions on the run's volume axis (m³).
    delay_volume_m3 : float
        Volume between simulated outlet and fraction collector, m³ (default 0).

    Returns
    -------
    pandas.DataFrame
        Columns ``fraction``, ``component``, ``concentration_mol_per_m3``.
    """
    collector_volume = result.volume_m3 + delay_volume_m3
    end_simulated = float(collector_volume[-1])
    rows = []
    for fraction in fractions:
        if fraction.end_m3 > end_simulated + 1e-12:
            raise ChromSimError(
                f"Fraction '{fraction.label}' ends at {fraction.end_m3 * 1e6:.2f} mL, after the "
                f"simulation ({end_simulated * 1e6:.2f} mL). Extend the method or check alignment."
            )
        for index, name in enumerate(result.component_names):
            rows.append(
                {
                    "fraction": fraction.label,
                    "component": name,
                    "concentration_mol_per_m3": interval_average(
                        collector_volume,
                        result.outlet_mol_per_m3[:, index],
                        fraction.start_m3,
                        fraction.end_m3,
                    ),
                }
            )
    return pd.DataFrame(rows, columns=["fraction", "component", "concentration_mol_per_m3"])


@dataclass
class AssayComparison:
    """Measured vs simulated fraction-averaged readouts.

    Attributes
    ----------
    table : pandas.DataFrame
        Columns ``fraction``, ``analyte``, ``unit``, ``measured``, ``sd``, ``below_loq``,
        ``simulated``, ``start_m3``, ``end_m3`` (measured and simulated in ``unit``).
    assay_fractions_not_collected : list of str
        Fraction ids in the assay without matching fraction (no volumes) – reported, not fixed.
    collected_fractions_without_assay : list of str
        Collected fractions that the assay does not report.
    """

    table: pd.DataFrame
    assay_fractions_not_collected: list[str] = field(default_factory=list)
    collected_fractions_without_assay: list[str] = field(default_factory=list)


def compare_with_assay(
    result: SimulationResult,
    assay: FractionAssay,
    fractions: list[Fraction] | None = None,
    delay_volume_m3: float = 0.0,
    warn: bool = True,
) -> AssayComparison:
    """Compare assay results with fraction-averaged simulated concentrations.

    Fraction volumes come from the assay file if it has them, otherwise from ``fractions``
    (normally ``experiment.fractions`` from the ÄKTA export), matched by label. Mismatches are
    reported in the result and logged, never silently fixed (SPEC §7.2b).

    Parameters
    ----------
    result : SimulationResult
        Simulation to compare.
    assay : FractionAssay
        Measured results; analytes are mapped to components through ``assay.readouts``.
    fractions : list of Fraction, optional
        Collected fractions, used when the assay has no volumes.
    delay_volume_m3 : float
        Volume between simulated outlet and fraction collector, m³.
    warn : bool
        Log fraction-id mismatches (they are always listed in the result).

    Returns
    -------
    AssayComparison
        Table of measured vs simulated values plus mismatch lists.
    """
    if assay.fraction_volumes:
        available = {
            label: Fraction(label, start, end)
            for label, (start, end) in assay.fraction_volumes.items()
        }
    elif fractions is not None:
        available = {fraction.label: fraction for fraction in fractions}
    else:
        raise ChromSimError(
            f"Assay '{assay.name}' has no fraction volumes; pass fractions=experiment.fractions."
        )

    assay_labels = assay.fractions
    not_collected = [label for label in assay_labels if label not in available]
    without_assay = [label for label in available if label not in assay_labels]
    if not_collected and warn:
        logger.warning(
            "Assay '%s': fraction(s) %s have no matching collected fraction and are left out.",
            assay.name,
            ", ".join(not_collected),
        )
    matched = [available[label] for label in assay_labels if label in available]
    averages = fraction_averages(result, matched, delay_volume_m3)
    average_of = {
        (row.fraction, row.component): row.concentration_mol_per_m3 for row in averages.itertuples()
    }

    molar_masses = dict(zip(result.component_names, result.molar_masses_kg_per_mol, strict=True))
    rows = []
    for row in assay.results.itertuples():
        if row.fraction not in available or row.analyte not in assay.readouts:
            continue
        components = assay.readouts[row.analyte]
        unknown = [c for c in components if c not in molar_masses]
        if unknown:
            raise ChromSimError(
                f"Readout '{row.analyte}' refers to unknown component(s) {', '.join(unknown)}; "
                f"simulated components are {', '.join(result.component_names)}."
            )
        simulated = sum(
            units.concentration_from_si(average_of[(row.fraction, c)], row.unit, molar_masses[c])
            for c in components
        )
        fraction = available[row.fraction]
        rows.append(
            {
                "fraction": row.fraction,
                "analyte": row.analyte,
                "unit": row.unit,
                "measured": row.value,
                "sd": row.sd,
                "below_loq": row.below_loq,
                "simulated": float(simulated),
                "start_m3": fraction.start_m3,
                "end_m3": fraction.end_m3,
            }
        )
    columns = [
        "fraction",
        "analyte",
        "unit",
        "measured",
        "sd",
        "below_loq",
        "simulated",
        "start_m3",
        "end_m3",
    ]
    return AssayComparison(
        table=pd.DataFrame(rows, columns=columns),
        assay_fractions_not_collected=not_collected,
        collected_fractions_without_assay=without_assay,
    )
