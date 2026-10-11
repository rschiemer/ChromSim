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


# --- pooling and KPIs (SPEC §7.4) ---------------------------------------------------------------


@dataclass(frozen=True)
class Pool:
    """A product pool: everything leaving the outlet between two cut times.

    Attributes
    ----------
    target : str
        Product component.
    start_s, end_s : float
        Cut points, s from the start of the method.
    """

    target: str
    start_s: float
    end_s: float


@dataclass
class PoolPerformance:
    """Key performance indicators of a pool.

    Attributes
    ----------
    pool : Pool
    amounts_mol : dict[str, float]
        Amount of every component in the pool, mol.
    volume_m3 : float
        Pool volume, m³.
    purity : float
        Target / (target + impurities), by amount of substance (by mass when every compared
        component has a molecular weight), dimensionless.
    recovery : float
        Yield: target in the pool / target fed in one cycle, dimensionless.
    concentration_mol_per_m3 : float
        Target concentration in the pool, mol/m³.
    productivity_mol_per_m3_s : float
        Target in the pool per cycle time and packed-bed volume, mol/(m³ s).
    buffer_consumption_m3_per_mol : float
        Total pumped volume of a cycle per target in the pool, m³/mol (inf if nothing pooled).
    target_molar_mass_kg_per_mol : float or None
        For the mass-based properties below.
    cycle_time_s : float
        Method duration, s.
    """

    pool: Pool
    amounts_mol: dict[str, float]
    volume_m3: float
    purity: float
    recovery: float
    concentration_mol_per_m3: float
    productivity_mol_per_m3_s: float
    buffer_consumption_m3_per_mol: float
    target_molar_mass_kg_per_mol: float | None
    cycle_time_s: float

    @property
    def productivity_kg_per_m3_s(self) -> float | None:
        """Productivity by mass, kg/(m³ s) (= g/(L s)); ``None`` without molecular weight."""
        if self.target_molar_mass_kg_per_mol is None:
            return None
        return self.productivity_mol_per_m3_s * self.target_molar_mass_kg_per_mol

    @property
    def buffer_consumption_m3_per_kg(self) -> float | None:
        """Pumped volume per product mass, m³/kg (= L/g); ``None`` without molecular weight."""
        if self.target_molar_mass_kg_per_mol is None:
            return None
        return self.buffer_consumption_m3_per_mol / self.target_molar_mass_kg_per_mol


def fed_amounts(result: SimulationResult) -> dict[str, float]:
    """Amount of every component pumped into the system during the method, mol.

    Each step contributes Q Δt (c_start + c_end)/2: its inlet is linear in time.
    """
    totals = np.zeros(len(result.component_names))
    for step in result.steps:
        average = 0.5 * (
            np.array(step.inlet_start_mol_per_m3) + np.array(step.inlet_end_mol_per_m3)
        )
        totals += step.flow_rate_m3_per_s * step.duration_s * average
    return dict(zip(result.component_names, totals.tolist(), strict=True))


def default_impurities(project, target: str) -> list[str]:
    """All components except the target and the isotherm's modifiers (salt, pH)."""
    from chromsim.isotherms.registry import get_isotherm  # registry imports nothing from analysis

    modifiers = get_isotherm(project.isotherm.model).modifier_components(project.component_names)
    return [name for name in project.component_names if name != target and name not in modifiers]


def _segment_amounts(result: SimulationResult, component: str) -> np.ndarray:
    """Amount leaving the outlet between consecutive output times, mol (trapezoid in volume)."""
    c = _outlet(result, component)
    return 0.5 * (c[1:] + c[:-1]) * np.diff(result.volume_m3)


def _outlet(result: SimulationResult, component: str) -> np.ndarray:
    """Outlet concentration with numerical undershoots (tiny negative values) set to zero.

    Solvers return values like −1e-12 mol/m³ where nothing elutes; summed over a pool they would
    count as negative impurity and inflate purity.
    """
    return np.maximum(result.outlet_mol_per_m3[:, result.component_names.index(component)], 0.0)


def pool_performance(
    result: SimulationResult,
    pool: Pool,
    impurities: list[str],
    column_volume_m3: float | None = None,
) -> PoolPerformance:
    """KPIs of a pool.

    Parameters
    ----------
    result : SimulationResult
        Simulated cycle.
    pool : Pool
        Cut points and target.
    impurities : list of str
        Components counted against purity (see :func:`default_impurities`).
    column_volume_m3 : float, optional
        Packed-bed volume for productivity; default ``result.column_volume_m3``.

    Returns
    -------
    PoolPerformance

    Notes
    -----
    Definitions (Schmidt-Traub, Schulte & Seidel-Morgenstern (eds.), *Preparative
    Chromatography*, 2nd ed., Wiley-VCH 2012, ch. 7): recovery = m_pool/m_fed; purity =
    m_target/(m_target + Σ m_impurities); productivity = m_pool/(t_cycle V_c); specific buffer
    consumption = V_pumped/m_pool.
    """
    start_v, end_v = np.interp([pool.start_s, pool.end_s], result.time_s, result.volume_m3)
    amounts = {}
    for name in result.component_names:
        if end_v > start_v:
            average = interval_average(result.volume_m3, _outlet(result, name), start_v, end_v)
            amounts[name] = float(average * (end_v - start_v))
        else:
            amounts[name] = 0.0
    target = amounts[pool.target]
    impurity = sum(amounts[name] for name in impurities)
    molar_masses = dict(zip(result.component_names, result.molar_masses_kg_per_mol, strict=True))
    compared = [pool.target, *impurities]
    known = {name: m for name, m in molar_masses.items() if m is not None}
    if all(name in known for name in compared):  # purity by mass
        target_mass = target * known[pool.target]
        total_mass = target_mass + sum(amounts[n] * known[n] for n in impurities)
        purity = target_mass / total_mass if total_mass > 0 else 0.0
    else:
        purity = target / (target + impurity) if target + impurity > 0 else 0.0
    cycle_time = result.steps[-1].end_s - result.steps[0].start_s
    pumped = float(result.volume_m3[-1] - result.volume_m3[0])
    volume = float(end_v - start_v)
    fed = fed_amounts(result)[pool.target]
    column_volume = column_volume_m3 or result.column_volume_m3
    return PoolPerformance(
        pool=pool,
        amounts_mol=amounts,
        volume_m3=volume,
        purity=float(purity),
        recovery=float(target / fed) if fed > 0 else 0.0,
        concentration_mol_per_m3=float(target / volume) if volume > 0 else 0.0,
        productivity_mol_per_m3_s=float(target / (cycle_time * column_volume)),
        buffer_consumption_m3_per_mol=float(pumped / target) if target > 0 else float("inf"),
        target_molar_mass_kg_per_mol=molar_masses[pool.target],
        cycle_time_s=float(cycle_time),
    )


def weighted_segments(
    result: SimulationResult, target: str, impurities: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Target and total impurity leaving the outlet per output interval.

    By mass when every compared component has a molecular weight, else by amount of substance
    (the same basis as :attr:`PoolPerformance.purity`).
    """
    molar_masses = dict(zip(result.component_names, result.molar_masses_kg_per_mol, strict=True))
    by_mass = all(molar_masses[n] is not None for n in [target, *impurities])

    def amount(name):
        weight = molar_masses[name] if by_mass else 1.0
        return _segment_amounts(result, name) * weight

    a = amount(target)
    b = sum((amount(name) for name in impurities), np.zeros_like(a))
    return a, b


def best_interval_purity(result: SimulationResult, target: str, impurities: list[str]) -> float:
    """Highest purity of any output interval that carries ≥ 0.1 % of the largest target amount.

    An upper bound for what any pool can reach; 0 if no target leaves the column.
    """
    a, b = weighted_segments(result, target, impurities)
    if a.max() <= 0:
        return 0.0
    carrying = a > 1e-3 * a.max()
    return float(np.max(a[carrying] / (a[carrying] + b[carrying])))


def optimal_pool(
    result: SimulationResult,
    target: str,
    impurities: list[str],
    purity: float,
    pool_steps: tuple[str, ...] | None = None,
) -> Pool | None:
    """Cut points that maximise recovery subject to purity ≥ ``purity`` (exact on the output grid).

    For one contiguous pool between output times t_i < t_j, with a_k the target amount and b_k the
    impurity amount (mass when all have molecular weights) leaving in segment k,

    .. math:: \\frac{\\sum a_k}{\\sum a_k + \\sum b_k} \\ge P \\iff
              \\sum_{k=i}^{j-1} w_k \\ge 0, \\qquad w_k = (1-P)\\,a_k - P\\,b_k,

    so the best pool maximises A_j − A_i subject to W_j − W_i ≥ 0 (prefix sums A, W). All pairs
    are checked, which finds the global optimum on the grid.

    ``pool_steps`` restricts collection to method steps by name (e.g. ``("elute",)``): otherwise
    the highest recovery may come from collecting traces during load or wash, which dilutes the
    pool.

    Returns
    -------
    Pool or None
        ``None`` if no pool reaches the purity.
    """
    if not 0.0 < purity < 1.0:
        raise ChromSimError(f"Required purity must be between 0 and 1, got {purity}.")
    a, b = weighted_segments(result, target, impurities)
    w = (1.0 - purity) * a - purity * b
    prefix_a = np.concatenate([[0.0], np.cumsum(a)])
    prefix_w = np.concatenate([[0.0], np.cumsum(w)])

    # Only segments carrying a noticeable amount of target, inside the allowed steps, can bound an
    # optimal pool.
    allowed = np.ones(a.size, dtype=bool)
    if pool_steps is not None:
        chosen = [step for step in result.steps if step.name in pool_steps]
        if not chosen:
            names = ", ".join(step.name for step in result.steps)
            raise ChromSimError(f"pool_steps {pool_steps} match no method step; steps are {names}.")
        start = min(step.start_s for step in chosen)
        end = max(step.end_s for step in chosen)
        allowed = (result.time_s[:-1] >= start) & (result.time_s[1:] <= end)
    carrying = (
        np.nonzero(allowed & (a > 1e-6 * a.max()))[0] if a.max() > 0 else np.array([], dtype=int)
    )
    if carrying.size == 0:
        return None
    starts, ends = carrying, carrying + 1
    best_gain, best = 0.0, None
    for i in starts:
        candidates = ends[ends > i]
        feasible = candidates[prefix_w[candidates] - prefix_w[i] >= 0.0]
        if feasible.size == 0:
            continue
        gains = prefix_a[feasible] - prefix_a[i]
        k = int(np.argmax(gains))
        if gains[k] > best_gain:
            best_gain, best = gains[k], (i, feasible[k])
    # A pool holding a negligible share of the target is no product pool.
    if best is None or best_gain < 1e-6 * prefix_a[-1]:
        return None
    i, j = best
    return Pool(target=target, start_s=float(result.time_s[i]), end_s=float(result.time_s[j]))
