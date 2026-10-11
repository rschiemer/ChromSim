"""Inverse fit of model parameters to one or more experiments (SPEC §7.3, use case U4).

The fit adjusts the chosen :class:`FitParameter` values until simulated detector signals and/or
fraction-assay results match the measurements, in the weighted least-squares sense:

.. math:: \\min_x \\; \\tfrac12 \\sum_{e} \\sum_{t \\in e} \\lVert r_{e,t}(x) \\rVert^2

Parameters are either shared by all experiments (e.g. isotherm parameters) or estimated per
experiment (e.g. a load concentration, ``per_experiment=True``). The search runs in the scaled
space x ∈ [0, 1] with scipy's trust-region-reflective least squares (finite-difference Jacobian).
``multi_start`` adds random starting points to reduce the risk of a local minimum.

Phase 5 routes this through the optimizer plugin interface (SPEC §6); the call signature stays.

Reference: Branch, Coleman & Li, SIAM J. Sci. Comput. 21 (1999) 1 (TRF algorithm).
"""

import logging
import time
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from chromsim.backends import get_backend
from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.errors import ChromSimError
from chromsim.estimation.objectives import FitExperiment
from chromsim.estimation.parameters import FitParameter, get_value, set_values

logger = logging.getLogger(__name__)

# Residual assigned to every point when a simulation fails, so the optimiser backs away.
FAILED_SIMULATION_RESIDUAL = 10.0
# Finite-difference step in the scaled space: 0.1 % of each parameter's range (log range for log
# parameters). Large against CADET's relative tolerance (1e-6), small against curvature.
DIFF_STEP = 1e-3


@dataclass
class FitResult:
    """Outcome of :func:`inverse_fit`.

    Attributes
    ----------
    parameters : list of FitParameter
        As given.
    experiment_names : list of str
        Experiments in order.
    values : dict[str, float]
        Fitted values by label; per-experiment parameters as ``"label[experiment]"``.
    x : numpy.ndarray
        Fitted point in the scaled space.
    cost : float
        ½ Σ r² at the optimum.
    residuals : numpy.ndarray
        All residuals at the optimum.
    jacobian : numpy.ndarray
        dr/dx at the optimum (scaled space), shape (n_residuals, n_x).
    projects, results : dict
        Fitted project and its simulation, per experiment.
    n_evaluations : int
        Residual evaluations (each simulates every experiment).
    success : bool
    message : str
    history : list of (float, float)
        (time since start in s, cost) per evaluation.
    starts : list of (numpy.ndarray, float)
        (start point, final cost) of every start.
    """

    parameters: list[FitParameter]
    experiment_names: list[str]
    values: dict[str, float]
    x: np.ndarray
    cost: float
    residuals: np.ndarray
    jacobian: np.ndarray
    projects: dict[str, Project]
    results: dict[str, SimulationResult]
    n_evaluations: int
    success: bool
    message: str
    history: list[tuple[float, float]] = field(default_factory=list)
    starts: list[tuple[np.ndarray, float]] = field(default_factory=list)

    @property
    def x_labels(self) -> list[str]:
        """Label of every entry of ``x``."""
        return list(self.values)

    def x_parameters(self) -> list[FitParameter]:
        """The FitParameter behind every entry of ``x``."""
        return _expand(self.parameters, self.experiment_names)[1]


def _expand(parameters, experiment_names):
    """Labels, parameters and owning experiment (None if shared) of every entry of x."""
    labels, owners, experiment_of = [], [], []
    for parameter in parameters:
        if parameter.per_experiment:
            for name in experiment_names:
                labels.append(f"{parameter.label}[{name}]")
                owners.append(parameter)
                experiment_of.append(name)
        else:
            labels.append(parameter.label)
            owners.append(parameter)
            experiment_of.append(None)
    return labels, owners, experiment_of


class _Objective:
    """Residual function over x; simulates every experiment with the parameters applied."""

    def __init__(self, experiments, parameters, callback):
        self.experiments = experiments
        self.parameters = parameters
        self.names = [experiment.name for experiment in experiments]
        self.labels, self.owners, self.experiment_of = _expand(parameters, self.names)
        self.callback = callback
        self.start = time.perf_counter()
        self.history: list[tuple[float, float]] = []
        self.n_evaluations = 0
        self.size: int | None = None

    def projects(self, x) -> dict[str, Project]:
        projects = {}
        for experiment in self.experiments:
            assignments = []
            for x_index, parameter in enumerate(self.owners):
                if self.experiment_of[x_index] not in (None, experiment.name):
                    continue
                assignments.append(
                    (parameter.path, parameter.from_scaled(x[x_index]), parameter.unit)
                )
            projects[experiment.name] = set_values(experiment.project, assignments)
        return projects

    def evaluate(self, x):
        projects = self.projects(x)
        results, pieces = {}, []
        for experiment in self.experiments:
            project = projects[experiment.name]
            result = get_backend(project.solver.backend).simulate(project)
            results[experiment.name] = result
            pieces += [target.residuals(result, project) for target in experiment.targets]
        return np.concatenate(pieces), projects, results

    def __call__(self, x):
        self.n_evaluations += 1
        try:
            residuals = self.evaluate(x)[0]
            self.size = residuals.size
        except ChromSimError as error:
            if self.size is None:
                raise
            logger.warning("Simulation failed at x=%s (%s); penalised.", np.round(x, 4), error)
            residuals = np.full(self.size, FAILED_SIMULATION_RESIDUAL)
        cost = 0.5 * float(residuals @ residuals)
        self.history.append((time.perf_counter() - self.start, cost))
        if self.callback is not None:
            self.callback(self.n_evaluations, cost)
        return residuals


def inverse_fit(
    experiments: list[FitExperiment],
    parameters: list[FitParameter],
    multi_start: int = 0,
    seed: int = 0,
    max_evaluations: int = 200,
    callback: Callable[[int, float], None] | None = None,
) -> FitResult:
    """Estimate parameters from one or more experiments.

    Parameters
    ----------
    experiments : list of FitExperiment
        Each with its project and targets (``SignalTarget``, ``AssayTarget``).
    parameters : list of FitParameter
        Parameters to estimate; start values default to the first experiment's project.
    multi_start : int
        Number of additional random starting points (uniform in the scaled space).
    seed : int
        Random seed for the extra starts.
    max_evaluations : int
        Maximum residual evaluations per start.
    callback : callable, optional
        Called as ``callback(evaluation_number, cost)`` after every evaluation (progress bars).

    Returns
    -------
    FitResult
        Best fit over all starts.
    """
    if not experiments or not parameters:
        raise ChromSimError("inverse_fit needs at least one experiment and one parameter.")
    names = [experiment.name for experiment in experiments]
    if len(set(names)) != len(names):
        raise ChromSimError("Experiment names must be unique.")

    objective = _Objective(experiments, parameters, callback)
    x0 = []
    for parameter, label, owner in zip(
        objective.owners, objective.labels, objective.experiment_of, strict=True
    ):
        if parameter.initial is not None:
            value = parameter.initial
        else:
            experiment = experiments[names.index(owner)] if owner is not None else experiments[0]
            value = get_value(experiment.project, parameter.path, parameter.unit)
        if not parameter.lower <= value <= parameter.upper:
            raise ChromSimError(
                f"Start value {value:g} of '{label}' is outside its bounds "
                f"[{parameter.lower:g}, {parameter.upper:g}]; set initial= or widen the bounds."
            )
        x0.append(parameter.to_scaled(value))

    rng = np.random.default_rng(seed)
    starting_points = [np.array(x0)] + [rng.uniform(0.0, 1.0, len(x0)) for _ in range(multi_start)]
    best, starts = None, []
    for start in starting_points:
        solution = least_squares(
            objective,
            start,
            bounds=(0.0, 1.0),
            method="trf",
            diff_step=DIFF_STEP,
            max_nfev=max_evaluations,
        )
        starts.append((start, float(solution.cost)))
        logger.info(
            "Start %s → cost %.4g (%s)", np.round(start, 3), solution.cost, solution.message
        )
        if best is None or solution.cost < best.cost:
            best = solution

    assert best is not None  # there is always at least one start
    residuals, projects, results = objective.evaluate(best.x)
    values = {
        label: float(parameter.from_scaled(x))
        for label, parameter, x in zip(objective.labels, objective.owners, best.x, strict=True)
    }
    return FitResult(
        parameters=list(parameters),
        experiment_names=names,
        values=values,
        x=best.x,
        cost=0.5 * float(residuals @ residuals),
        residuals=residuals,
        jacobian=best.jac,
        projects=projects,
        results=results,
        n_evaluations=objective.n_evaluations,
        success=bool(best.success),
        message=str(best.message),
        history=objective.history,
        starts=starts,
    )
