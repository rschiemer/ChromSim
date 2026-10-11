"""Process optimization on a fast two-protein column (SPEC §7.4 workflows tier).

A target A and a more strongly retained impurity B: a shallower gradient separates them better,
so at 99 % purity the recovery rises with gradient length (measured: 0.18 at 5 min, 0.41 at
10 min, 0.66 at 30 min, 0.76 at 120 min) and the productivity falls.
"""

import numpy as np
import pytest
from project_examples import make_binding_project

from chromsim import api
from chromsim.estimation.parameters import FitParameter

pytestmark = pytest.mark.slow

PURITY = 0.99


def separation_project():
    components = [{"name": "salt"}, {"name": "A", "mw_kDa": 150.0}, {"name": "B", "mw_kDa": 150.0}]
    isotherm = {
        "model": "SMA",
        "parameters": {
            "A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0},
            "B": {"k_eq": 3e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0},
        },
    }
    project = make_binding_project(
        components,
        isotherm,
        {"salt": "50 mM", "A": "0.1 mM", "B": "0.02 mM"},
        elute={"salt": "1000 mM"},
        n_col=40,
    )
    project.method[2].duration = "10 min"  # a short, steep gradient to start from
    return project


GRADIENT = FitParameter(
    path="method.2.duration", lower=5.0, upper=120.0, unit="min", name="gradient_min"
)


def test_single_objective_raises_recovery_at_fixed_purity():
    project = separation_project()
    start = api.simulate(project)
    impurities = api.default_impurities(project, "A")
    assert impurities == ["B"]
    start_pool = api.optimal_pool(start, "A", impurities, PURITY, pool_steps=("elute",))
    start_recovery = (
        0.0 if start_pool is None else api.pool_performance(start, start_pool, impurities).recovery
    )

    result = api.optimize_process(
        project, [GRADIENT], "A", PURITY, pool_steps=("elute",), max_evaluations=25
    )
    best = result.best.performance
    assert best is not None
    assert best.purity >= PURITY - 1e-9  # exact cut points guarantee the constraint
    # From 0.41 at the 10 min start, a clearly better method (≥ 0.6) by a shallower gradient.
    assert start_recovery == pytest.approx(0.41, abs=0.02)
    assert best.recovery > 0.6
    assert result.optimization.x[0] > 25.0


def test_pareto_front_trades_recovery_against_productivity():
    result = api.optimize_process(
        separation_project(),
        [GRADIENT],
        "A",
        PURITY,
        objectives=("recovery", "productivity"),
        pool_steps=("elute",),
        pop_size=10,
        n_gen=4,
    )
    front = result.table.sort_values("recovery")
    assert len(front) >= 3
    # Non-dominated: along increasing recovery, productivity must not increase.
    assert np.all(np.diff(front["productivity"].to_numpy()) <= 1e-15)
    assert front["recovery"].iloc[-1] > front["recovery"].iloc[0]
