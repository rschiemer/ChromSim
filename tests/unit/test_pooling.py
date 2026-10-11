"""Pools, KPIs and optimal cut points on synthetic chromatograms (exact expectations)."""

import numpy as np
import pytest

from chromsim.analysis.fractionation import (
    Pool,
    fed_amounts,
    optimal_pool,
    pool_performance,
)
from chromsim.backends.base import SimulationResult
from chromsim.core.method import ResolvedStep
from chromsim.errors import ChromSimError

Q = 1e-8  # m³/s
FEED_S = 10.0  # 10 s load of 1 mol/m³ of A and B → 1e-7 mol each


def two_peaks(separation_s=60.0, molar_masses=(None, None), sigma=8.0):
    """A at 100 s and B at 100 + separation s, each Gaussian with exactly the fed amount."""
    t = np.arange(0.0, 400.0, 0.5)
    fed = Q * FEED_S * 1.0

    def peak(centre):
        shape = np.exp(-0.5 * ((t - centre) / sigma) ** 2)
        return shape * fed / np.trapezoid(shape * Q, t)  # ∫ c Q dt = fed

    steps = [
        ResolvedStep("load", 0.0, FEED_S, Q, (1.0, 1.0), (1.0, 1.0)),
        ResolvedStep("elute", FEED_S, 400.0, Q, (0.0, 0.0), (0.0, 0.0)),
    ]
    return SimulationResult(
        component_names=["A", "B"],
        time_s=t,
        outlet_mol_per_m3=np.column_stack([peak(100.0), peak(100.0 + separation_s)]),
        volume_m3=Q * t,
        column_volume_m3=1e-6,
        steps=steps,
        molar_masses_kg_per_mol=list(molar_masses),
    )


def test_fed_amounts_with_gradient():
    result = two_peaks()
    result.steps[1] = ResolvedStep("elute", FEED_S, 400.0, Q, (0.0, 0.0), (2.0, 0.0))
    # load: Q·10 s·1; gradient 0 → 2 over 390 s: Q·390·1 (average 1)
    assert fed_amounts(result)["A"] == pytest.approx(Q * 10.0 + Q * 390.0)


def test_kpis_of_pools():
    result = two_peaks()
    everything = pool_performance(result, Pool("A", 0.0, 399.5), impurities=["B"])
    assert everything.recovery == pytest.approx(1.0, rel=1e-6)
    assert everything.purity == pytest.approx(0.5, rel=1e-6)  # equal amounts, by mol
    only_a = pool_performance(result, Pool("A", 60.0, 130.0), impurities=["B"])
    assert only_a.purity > 0.999 and only_a.recovery > 0.999
    # productivity = n_pool / (t_cycle V_c); buffer consumption = V_pumped / n_pool
    n = only_a.amounts_mol["A"]
    assert only_a.productivity_mol_per_m3_s == pytest.approx(n / (400.0 * 1e-6))
    assert only_a.buffer_consumption_m3_per_mol == pytest.approx(Q * 399.5 / n)
    assert only_a.volume_m3 == pytest.approx(Q * 70.0)
    assert only_a.concentration_mol_per_m3 == pytest.approx(n / (Q * 70.0))
    assert only_a.productivity_kg_per_m3_s is None


def test_purity_by_mass_when_molecular_weights_are_known():
    # B is twice as heavy: by mass, purity of the full pool is 1 / (1 + 2) = 1/3.
    result = two_peaks(molar_masses=(100.0, 200.0))
    performance = pool_performance(result, Pool("A", 0.0, 399.5), impurities=["B"])
    assert performance.purity == pytest.approx(1 / 3, rel=1e-6)
    assert performance.productivity_kg_per_m3_s == pytest.approx(
        performance.productivity_mol_per_m3_s * 100.0
    )


def brute_force(result, purity):
    """Best (recovery, i, j) over every pair of grid points (slow reference)."""
    best = (0.0, None)
    for i in range(0, len(result.time_s), 4):
        for j in range(i + 4, len(result.time_s), 4):
            kpi = pool_performance(result, Pool("A", result.time_s[i], result.time_s[j]), ["B"])
            if kpi.purity >= purity - 1e-12 and kpi.recovery > best[0]:
                best = (kpi.recovery, (i, j))
    return best[0]


@pytest.mark.parametrize("required", [0.9, 0.99])
def test_optimal_pool_matches_brute_force_on_overlapping_peaks(required):
    result = two_peaks(separation_s=25.0)  # strongly overlapping
    pool = optimal_pool(result, "A", ["B"], required)
    kpi = pool_performance(result, pool, ["B"])
    assert kpi.purity >= required - 1e-9
    # the exact search uses every grid point, the brute force every 4th: it is at least as good
    assert kpi.recovery >= brute_force(result, required) - 1e-12
    assert kpi.recovery < 1.0  # purity costs yield when peaks overlap


def test_optimal_pool_when_purity_is_impossible():
    result = two_peaks(separation_s=0.0)  # identical peaks: purity 50 % everywhere
    assert optimal_pool(result, "A", ["B"], 0.6) is None
    with pytest.raises(ChromSimError, match="between 0 and 1"):
        optimal_pool(result, "A", ["B"], 1.0)


def test_numerical_undershoots_do_not_fake_purity():
    # Tiny negative impurity values (solver undershoot) must not make a pool "purer than pure".
    result = two_peaks(separation_s=25.0)
    result.outlet_mol_per_m3[:20, 1] = -1e-12
    result.outlet_mol_per_m3[:20, 0] = 1e-15
    performance = pool_performance(result, Pool("A", 0.0, 9.5), impurities=["B"])
    assert performance.purity <= 1.0
    pool = optimal_pool(result, "A", ["B"], 0.995)
    assert pool is None or pool_performance(result, pool, ["B"]).recovery > 1e-3


def test_pool_steps_restrict_collection():
    result = two_peaks(separation_s=25.0)
    unrestricted = optimal_pool(result, "A", ["B"], 0.9)
    restricted = optimal_pool(result, "A", ["B"], 0.9, pool_steps=("elute",))
    assert restricted.start_s >= FEED_S
    assert (
        pool_performance(result, restricted, ["B"]).recovery
        <= pool_performance(result, unrestricted, ["B"]).recovery
    )
    with pytest.raises(ChromSimError, match="match no method step; steps are load, elute"):
        optimal_pool(result, "A", ["B"], 0.9, pool_steps=("strip",))
