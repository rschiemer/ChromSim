"""One optimizer interface for parameter estimation and process optimization (SPEC §6).

A problem is described once (:class:`OptimizationProblem`): variables with bounds, and either
objectives (one or several, all minimised) or residuals (least squares), plus optional
constraints. Every optimizer plugin receives it through :class:`ScaledProblem`, which

- maps each variable to x ∈ [0, 1] (linear, or log10 for ``log_scale`` variables), so optimizers
  never see units or decades of magnitude;
- evaluates batches of points, in parallel processes when ``problem.parallel > 1``;
- records every evaluation (x, f, time) and calls the progress callback.

Objectives, residuals and constraints receive the **physical** variable values and must be
picklable (module-level functions or objects of module-level classes) for parallel runs.

Constraint convention: g(x) ≤ 0 is feasible.
"""

import logging
import multiprocessing
import time
from abc import ABC, abstractmethod
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np

from chromsim.errors import ChromSimError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Variable:
    """A decision variable.

    Attributes
    ----------
    name : str
        Label for reports.
    lower, upper : float
        Bounds (physical units of the problem).
    log_scale : bool
        Search in log10 (for rate constants and other positive quantities spanning decades).
    initial : float, optional
        Start value for local optimizers; default the middle of the scaled range.
    """

    name: str
    lower: float
    upper: float
    log_scale: bool = False
    initial: float | None = None

    def __post_init__(self):
        if not self.lower < self.upper:
            raise ChromSimError(f"Variable '{self.name}': lower bound must be below upper bound.")
        if self.log_scale and self.lower <= 0:
            raise ChromSimError(f"Variable '{self.name}': log scale needs a positive lower bound.")

    def to_scaled(self, value: float) -> float:
        """Physical value → x ∈ [0, 1]."""
        if self.log_scale:
            return float(
                (np.log10(value) - np.log10(self.lower))
                / (np.log10(self.upper) - np.log10(self.lower))
            )
        return float((value - self.lower) / (self.upper - self.lower))

    def from_scaled(self, x: float) -> float:
        """x ∈ [0, 1] → physical value."""
        if self.log_scale:
            return float(
                10.0 ** (np.log10(self.lower) + x * (np.log10(self.upper) - np.log10(self.lower)))
            )
        return float(self.lower + x * (self.upper - self.lower))


@dataclass(frozen=True)
class Constraint:
    """Inequality constraint g(x) ≤ 0 on the physical variables (feasible when ≤ 0)."""

    name: str
    function: Callable[[np.ndarray], float]


@dataclass
class OptimizationProblem:
    """What to optimise.

    Attributes
    ----------
    variables : list of Variable
    objectives : callable, optional
        x (physical) → array of ``n_objectives`` values to minimise.
    n_objectives : int
    constraints : list of Constraint
    residuals : callable, optional
        x (physical) → residual vector; least-squares optimizers minimise ½‖r‖². Single-objective
        optimizers use ½‖r‖² as objective when no ``objectives`` are given.
    parallel : int
        Number of processes for batch evaluations (1 = serial).
    """

    variables: list[Variable]
    objectives: Callable[[np.ndarray], np.ndarray] | None = None
    n_objectives: int = 1
    constraints: list[Constraint] = field(default_factory=list)
    residuals: Callable[[np.ndarray], np.ndarray] | None = None
    parallel: int = 1

    def __post_init__(self):
        if self.objectives is None and self.residuals is None:
            raise ChromSimError("An OptimizationProblem needs objectives or residuals.")
        if self.n_objectives > 1 and self.objectives is None:
            raise ChromSimError("Multi-objective problems need an objectives function.")


@dataclass
class Evaluation:
    """One evaluation: physical x, objective values (½‖r‖² for residual problems), constraint
    values, and seconds since the start."""

    x: np.ndarray
    f: np.ndarray
    g: np.ndarray
    time_s: float


@dataclass
class OptimizationResult:
    """Outcome of an optimization.

    Attributes
    ----------
    optimizer : str
    variable_names : list of str
    x : numpy.ndarray
        Best point (physical). For multi-objective problems: one Pareto point (lowest first
        objective); see ``pareto_x``.
    f : numpy.ndarray
        Objectives at ``x``.
    success : bool
    message : str
    history : list of Evaluation
        Every evaluation in order.
    pareto_x, pareto_f : numpy.ndarray, optional
        Non-dominated points (multi-objective optimizers).
    residuals, jacobian : numpy.ndarray, optional
        At the optimum, for least-squares optimizers; the Jacobian is dr/dx in the scaled space.
    """

    optimizer: str
    variable_names: list[str]
    x: np.ndarray
    f: np.ndarray
    success: bool
    message: str
    history: list[Evaluation] = field(default_factory=list)
    pareto_x: np.ndarray | None = None
    pareto_f: np.ndarray | None = None
    residuals: np.ndarray | None = None
    jacobian: np.ndarray | None = None

    @property
    def n_evaluations(self) -> int:
        """Number of evaluations."""
        return len(self.history)


def _evaluate(problem: OptimizationProblem, x: np.ndarray):
    """Objectives (or ½‖r‖²), constraint values and residuals at physical x. Runs in workers."""
    residuals = None
    if problem.residuals is not None:
        residuals = np.asarray(problem.residuals(x), dtype=float)
    if problem.objectives is not None:
        f = np.atleast_1d(np.asarray(problem.objectives(x), dtype=float))
    else:
        assert residuals is not None  # guaranteed by OptimizationProblem.__post_init__
        f = np.array([0.5 * residuals @ residuals])
    g = np.array([constraint.function(x) for constraint in problem.constraints], dtype=float)
    return f, g, residuals


class ScaledProblem:
    """The problem as optimizers see it: x ∈ [0, 1]ⁿ, with history, callback and parallelism."""

    def __init__(
        self,
        problem: OptimizationProblem,
        callback: Callable[[int, np.ndarray], None] | None = None,
    ):
        self.problem = problem
        self.callback = callback
        self.history: list[Evaluation] = []
        self.start = time.perf_counter()
        self._executor: ProcessPoolExecutor | None = None
        self._cache: dict[bytes, tuple] = {}

    @property
    def n_variables(self) -> int:
        return len(self.problem.variables)

    def to_physical(self, x_scaled) -> np.ndarray:
        return np.array(
            [v.from_scaled(xi) for v, xi in zip(self.problem.variables, x_scaled, strict=True)]
        )

    def to_scaled(self, x_physical) -> np.ndarray:
        return np.array(
            [v.to_scaled(xi) for v, xi in zip(self.problem.variables, x_physical, strict=True)]
        )

    def start_point(self) -> np.ndarray:
        """Scaled start point from the variables' ``initial`` values (middle of range otherwise)."""
        return np.array(
            [0.5 if v.initial is None else v.to_scaled(v.initial) for v in self.problem.variables]
        )

    def evaluate_batch(
        self, points_scaled
    ) -> list[tuple[np.ndarray, np.ndarray, np.ndarray | None]]:
        """Evaluate several scaled points (in parallel if configured) and record them.

        Points evaluated before are taken from a cache, so an optimizer that asks for objectives
        and constraints separately does not simulate twice.
        """
        points = [np.clip(np.asarray(p, dtype=float), 0.0, 1.0) for p in points_scaled]
        keys = [p.tobytes() for p in points]
        new = [i for i, key in enumerate(keys) if key not in self._cache]
        physical = [self.to_physical(points[i]) for i in new]
        if self.problem.parallel > 1 and len(new) > 1:
            if self._executor is None:
                context = multiprocessing.get_context("spawn")  # safe on Windows and macOS
                self._executor = ProcessPoolExecutor(self.problem.parallel, mp_context=context)
            outputs = list(self._executor.map(_evaluate, [self.problem] * len(new), physical))
        else:
            outputs = [_evaluate(self.problem, x) for x in physical]
        for i, x, output in zip(new, physical, outputs, strict=True):
            self._cache[keys[i]] = output
            self.history.append(
                Evaluation(x=x, f=output[0], g=output[1], time_s=time.perf_counter() - self.start)
            )
            if self.callback is not None:
                self.callback(len(self.history), output[0])
        return [self._cache[key] for key in keys]

    def objectives(self, x_scaled) -> np.ndarray:
        return self.evaluate_batch([x_scaled])[0][0]

    def residuals(self, x_scaled) -> np.ndarray:
        output = self.evaluate_batch([x_scaled])[0]
        if output[2] is None:
            raise ChromSimError("This optimizer needs a problem with residuals (least squares).")
        return output[2]

    def close(self) -> None:
        """Shut down worker processes."""
        if self._executor is not None:
            self._executor.shutdown()
            self._executor = None


class OptimizerBase(ABC):
    """Base class of optimizer plugins.

    Subclasses set ``name`` and ``capabilities`` and implement :meth:`optimize`, using
    :meth:`prepare` to get the :class:`ScaledProblem`. Capabilities: ``"bounds"``,
    ``"least_squares"``, ``"multi_objective"``, ``"constraints"``, ``"parallel"``, ``"global"``.
    Options are given at construction, e.g. ``get_optimizer("nsga2", pop_size=40, n_gen=50)``.
    """

    name: ClassVar[str]
    capabilities: ClassVar[set[str]]
    default_options: ClassVar[dict] = {}

    def __init__(self, **options):
        unknown = set(options) - set(self.default_options)
        if unknown:
            raise ChromSimError(
                f"Optimizer '{self.name}' has no option(s) {', '.join(sorted(unknown))}; "
                f"options: {', '.join(self.default_options) or 'none'}."
            )
        self.options = dict(self.default_options) | options

    def prepare(self, problem: OptimizationProblem, callback=None) -> ScaledProblem:
        """Check the problem against the capabilities and wrap it in the scaled space."""
        if problem.n_objectives > 1 and "multi_objective" not in self.capabilities:
            raise ChromSimError(
                f"Optimizer '{self.name}' handles one objective only; use e.g. 'nsga2'."
            )
        if problem.constraints and "constraints" not in self.capabilities:
            raise ChromSimError(
                f"Optimizer '{self.name}' cannot handle constraints; "
                "use e.g. 'differential_evolution' or 'nsga2'."
            )
        if "least_squares" in self.capabilities and problem.residuals is None:
            raise ChromSimError(f"Optimizer '{self.name}' needs a problem with residuals.")
        if problem.parallel > 1 and "parallel" not in self.capabilities:
            logger.warning(
                "Optimizer '%s' evaluates serially; parallel=%d ignored.",
                self.name,
                problem.parallel,
            )
        return ScaledProblem(problem, callback)

    def result(
        self, scaled: ScaledProblem, x_scaled, success: bool, message: str, **extra
    ) -> OptimizationResult:
        """Build the result for the best scaled point and close worker processes."""
        scaled.close()
        x = scaled.to_physical(x_scaled)
        best = min(scaled.history, key=lambda e: float(np.sum(np.abs(e.x - x))))
        return OptimizationResult(
            optimizer=self.name,
            variable_names=[v.name for v in scaled.problem.variables],
            x=x,
            f=best.f,
            success=bool(success),
            message=str(message),
            history=scaled.history,
            **extra,
        )

    @abstractmethod
    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        """Solve the problem. ``callback(evaluation_number, f)`` is called after each evaluation."""
