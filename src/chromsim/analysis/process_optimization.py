"""Process optimization: choose method settings for the best pooled product (SPEC §7.4, U6).

Decision variables are project entries addressed by path, as in fitting (e.g. the gradient
length ``method.3.duration`` in CV and the end salt ``method.3.inlet_end.salt`` in mM). For every
candidate the process is simulated and the **cut points are chosen exactly** by
:func:`~chromsim.analysis.fractionation.optimal_pool`: the pool with the highest recovery that
meets the purity requirement. The outer optimizer therefore only searches the method settings,
and every reported point satisfies the purity constraint by construction.

Objectives (all from :class:`~chromsim.analysis.fractionation.PoolPerformance`):

- ``recovery`` (maximised), ``productivity`` (maximised, per packed-bed volume),
- ``buffer_consumption`` (minimised), ``concentration`` (maximised, product in the pool).

If no pool reaches the purity for a candidate, all objectives get their worst value (0
recovery …) and the constraint g = P_required − P_best_sample > 0 reports by how much it failed,
so constraint-handling optimizers steer back.
"""

import logging
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np
import pandas as pd

from chromsim.analysis.fractionation import (
    Pool,
    PoolPerformance,
    best_interval_purity,
    default_impurities,
    optimal_pool,
    pool_performance,
)
from chromsim.backends import get_backend
from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.errors import ChromSimError
from chromsim.estimation.parameters import FitParameter, get_value, set_values
from chromsim.optimizers import (
    Constraint,
    OptimizationProblem,
    OptimizationResult,
    Variable,
    get_optimizer,
)

logger = logging.getLogger(__name__)

# Objective name → (PoolPerformance attribute, sign); optimizers minimise sign × value.
OBJECTIVES = MappingProxyType(
    {
        "recovery": ("recovery", -1.0),
        "productivity": ("productivity_mol_per_m3_s", -1.0),
        "buffer_consumption": ("buffer_consumption_m3_per_mol", 1.0),
        "concentration": ("concentration_mol_per_m3", -1.0),
    }
)
# Objective values for candidates without any pool reaching the purity (worst possible).
INFEASIBLE = MappingProxyType(
    {
        "recovery": 0.0,
        "productivity": 0.0,
        "buffer_consumption": 1e30,
        "concentration": 0.0,
    }
)


@dataclass
class ProcessEvaluation:
    """Outcome of one candidate: project, simulation, best pool and its KPIs (or why none)."""

    project: Project
    result: SimulationResult
    pool: Pool | None
    performance: PoolPerformance | None
    best_sample_purity: float


class ProcessEvaluator:
    """Simulates a candidate and pools it optimally; picklable for parallel optimizers.

    The last evaluation is cached, so objectives and the purity constraint of one candidate
    share one simulation.
    """

    def __init__(
        self,
        project: Project,
        parameters: list[FitParameter],
        target: str,
        purity: float,
        objectives: list[str],
        impurities: list[str] | None = None,
        pool_steps: tuple[str, ...] | None = None,
    ):
        unknown = [name for name in objectives if name not in OBJECTIVES]
        if unknown:
            raise ChromSimError(
                f"Unknown objective(s) {', '.join(unknown)}; choose from {', '.join(OBJECTIVES)}."
            )
        if target not in project.component_names:
            raise ChromSimError(
                f"Target '{target}' is not a component ({', '.join(project.component_names)})."
            )
        self.project = project
        self.parameters = parameters
        self.target = target
        self.purity = purity
        self.objectives = list(objectives)
        self.impurities = (
            impurities if impurities is not None else default_impurities(project, target)
        )
        self.pool_steps = pool_steps
        self._cache: tuple[bytes, ProcessEvaluation] | None = None

    def evaluate(self, x) -> ProcessEvaluation:
        """Simulate the candidate with physical values ``x`` and pool it optimally."""
        x = np.asarray(x, dtype=float)
        if self._cache is not None and self._cache[0] == x.tobytes():
            return self._cache[1]
        project = set_values(
            self.project,
            [(p.path, float(v), p.unit) for p, v in zip(self.parameters, x, strict=True)],
        )
        result = get_backend(project.solver.backend).simulate(project)
        pool = optimal_pool(result, self.target, self.impurities, self.purity, self.pool_steps)
        performance = None if pool is None else pool_performance(result, pool, self.impurities)
        best_purity = best_interval_purity(result, self.target, self.impurities)
        evaluation = ProcessEvaluation(project, result, pool, performance, best_purity)
        self._cache = (x.tobytes(), evaluation)
        return evaluation

    def objective_values(self, x) -> np.ndarray:
        """Values to minimise, in the order of ``objectives``."""
        evaluation = self.evaluate(x)
        values = []
        for name in self.objectives:
            attribute, sign = OBJECTIVES[name]
            raw = (
                INFEASIBLE[name]
                if evaluation.performance is None
                else getattr(evaluation.performance, attribute)
            )
            values.append(sign * raw)
        return np.array(values)

    def purity_constraint(self, x) -> float:
        """g ≤ 0 when a pool reaches the purity; else the purity shortfall of the best interval."""
        evaluation = self.evaluate(x)
        if evaluation.pool is not None:
            return 0.0
        return self.purity - evaluation.best_sample_purity


@dataclass
class ProcessOptimizationResult:
    """Outcome of :func:`optimize_process`.

    Attributes
    ----------
    optimization : OptimizationResult
        Raw result (history, Pareto set …).
    best : ProcessEvaluation
        The chosen point (for several objectives: the Pareto point with the highest first
        objective), re-simulated in the main process.
    table : pandas.DataFrame
        One row per evaluated (single objective) or Pareto (several objectives) point: variable
        values and objective values in their natural sign.
    """

    optimization: OptimizationResult
    best: ProcessEvaluation
    table: pd.DataFrame = field(default_factory=pd.DataFrame)


def optimize_process(
    project: Project,
    parameters: list[FitParameter],
    target: str,
    purity: float,
    objectives: tuple[str, ...] = ("recovery",),
    impurities: list[str] | None = None,
    pool_steps: tuple[str, ...] | None = None,
    optimizer: str | None = None,
    parallel: int = 1,
    callback=None,
    **options,
) -> ProcessOptimizationResult:
    """Optimise method settings for the best product pool at a required purity.

    Parameters
    ----------
    project : Project
        Calibrated project; the starting point is its current values (or ``initial``).
    parameters : list of FitParameter
        Decision variables (paths, bounds, units, lin/log scale).
    target : str
        Product component.
    purity : float
        Required purity of the pool (fraction, e.g. 0.99).
    objectives : tuple of str
        One or more of ``recovery``, ``productivity``, ``buffer_consumption``, ``concentration``.
    impurities : list of str, optional
        Components counted against purity; default all except target, salt and pH.
    pool_steps : tuple of str, optional
        Method steps during which the pool may be collected, e.g. ``("elute",)``.
    optimizer : str, optional
        Default ``nelder_mead`` for one objective, ``nsga2`` for several.
    parallel : int
        Processes for simultaneous simulations.
    callback : callable, optional
        ``callback(evaluation_number, objectives)``.
    **options
        Optimizer options (e.g. ``max_evaluations=40``, ``pop_size=16, n_gen=8``).

    Returns
    -------
    ProcessOptimizationResult
    """
    if not 0.0 < purity < 1.0:
        raise ChromSimError("purity must be a fraction between 0 and 1, e.g. 0.99.")
    evaluator = ProcessEvaluator(
        project, parameters, target, purity, list(objectives), impurities, pool_steps
    )
    variables = []
    for parameter in parameters:
        start = (
            parameter.initial
            if parameter.initial is not None
            else get_value(project, parameter.path, parameter.unit)
        )
        variables.append(
            Variable(
                parameter.label,
                parameter.lower,
                parameter.upper,
                log_scale=parameter.scale == "log",
                initial=min(max(start, parameter.lower), parameter.upper),
            )
        )
    problem = OptimizationProblem(
        variables,
        objectives=evaluator.objective_values,
        n_objectives=len(objectives),
        constraints=[Constraint("purity", evaluator.purity_constraint)],
        parallel=parallel,
    )
    name = optimizer or ("nelder_mead" if len(objectives) == 1 else "nsga2")
    chosen = get_optimizer(name, **options)
    if "constraints" not in chosen.capabilities:
        # Cut points already enforce purity; infeasible candidates get the worst objective values.
        problem.constraints = []
    outcome = chosen.optimize(problem, callback=callback)

    def natural(f):
        return {o: OBJECTIVES[o][1] * value for o, value in zip(objectives, f, strict=True)}

    if outcome.pareto_x is not None and len(objectives) > 1:
        rows = [
            dict(zip(outcome.variable_names, x, strict=True)) | natural(f)
            for x, f in zip(outcome.pareto_x, outcome.pareto_f, strict=True)
        ]
    else:
        rows = [
            dict(zip(outcome.variable_names, e.x, strict=True)) | natural(e.f)
            for e in outcome.history
        ]
    best = ProcessEvaluator(
        project, parameters, target, purity, list(objectives), impurities, pool_steps
    ).evaluate(outcome.x)
    if best.pool is None:
        logger.warning("No settings found that reach purity %.3g for '%s'.", purity, target)
    return ProcessOptimizationResult(optimization=outcome, best=best, table=pd.DataFrame(rows))


def process_variable_candidates(project: Project) -> list[FitParameter]:
    """Method settings that are commonly optimised, with default bounds (half to double now).

    Every step duration (in its own unit) and, for gradient steps, the end concentration of each
    modifier (salt, pH) that changes.
    """
    from chromsim.units import Quantity

    candidates = []
    for index, step in enumerate(project.method):
        duration = Quantity.parse(step.duration)
        candidates.append(
            FitParameter(
                path=f"method.{index}.duration",
                lower=duration.magnitude / 2,
                upper=duration.magnitude * 2,
                unit=duration.unit,
                name=f"{step.name} duration ({duration.unit})",
            )
        )
        for component, end in (step.inlet_end or {}).items():
            quantity = Quantity.parse(end)
            if quantity.magnitude <= 0:
                continue
            candidates.append(
                FitParameter(
                    path=f"method.{index}.inlet_end.{component}",
                    lower=quantity.magnitude / 2,
                    upper=quantity.magnitude * 2,
                    unit=quantity.unit,
                    name=f"{step.name} end {component} ({quantity.unit})",
                )
            )
    return candidates
