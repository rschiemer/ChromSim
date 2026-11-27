"""Inverse fit of model parameters to one or more experiments (SPEC §7.3, use case U4).

The fit adjusts the chosen :class:`FitParameter` values until simulated detector signals and/or
fraction-assay results match the measurements, in the weighted least-squares sense:

.. math:: \\min_x \\; \\tfrac12 \\sum_{e} \\sum_{t \\in e} \\lVert r_{e,t}(x) \\rVert^2

Parameters are either shared by all experiments (e.g. isotherm parameters) or estimated per
experiment (e.g. a load concentration, ``per_experiment=True``). The fit is an
:class:`~chromsim.optimizers.OptimizationProblem` solved by any registered optimizer (default
``least_squares``: scipy's trust-region-reflective method), in the scaled space x ∈ [0, 1]. With
``parallel=N`` the finite-difference Jacobian columns are simulated in N processes.
``multi_start`` adds random starting points to reduce the risk of a local minimum.

Reference: Branch, Coleman & Li, SIAM J. Sci. Comput. 21 (1999) 1 (TRF algorithm).
"""

import logging
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from chromsim.backends import get_backend
from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.errors import ChromSimError
from chromsim.estimation.objectives import FitExperiment
from chromsim.estimation.parameters import FitParameter, get_value, set_values
from chromsim.optimizers import OptimizationProblem, Variable, get_optimizer
from chromsim.optimizers.base import ScaledProblem
from chromsim.optimizers.scipy_opt import finite_difference_jacobian

logger = logging.getLogger(__name__)

# Residual assigned to every point when a simulation fails, so the optimiser backs away.
FAILED_SIMULATION_RESIDUAL = 10.0


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


class FitResiduals:
    """Residuals of all experiments for physical parameter values (picklable for workers)."""

    def __init__(self, experiments, owners, experiment_of):
        self.experiments = experiments
        self.owners = owners
        self.experiment_of = experiment_of
        self.size: int | None = None  # set after the first successful evaluation

    def projects(self, theta) -> dict[str, Project]:
        projects = {}
        for experiment in self.experiments:
            assignments = [
                (parameter.path, float(value), parameter.unit)
                for parameter, owner, value in zip(
                    self.owners, self.experiment_of, theta, strict=True
                )
                if owner in (None, experiment.name)
            ]
            projects[experiment.name] = set_values(experiment.project, assignments)
        return projects

    def evaluate(self, theta):
        """Residuals, fitted projects and simulation results at physical parameter values."""
        projects = self.projects(theta)
        results, pieces = {}, []
        for experiment in self.experiments:
            project = projects[experiment.name]
            result = get_backend(project.solver.backend).simulate(project)
            results[experiment.name] = result
            pieces += [target.residuals(result, project) for target in experiment.targets]
        return np.concatenate(pieces), projects, results

    def __call__(self, theta) -> np.ndarray:
        try:
            residuals = self.evaluate(theta)[0]
        except ChromSimError as error:
            if self.size is None:
                raise
            logger.warning("Simulation failed at %s (%s); penalised.", np.round(theta, 6), error)
            return np.full(self.size, FAILED_SIMULATION_RESIDUAL)
        self.size = residuals.size
        return residuals


def inverse_fit(
    experiments: list[FitExperiment],
    parameters: list[FitParameter],
    multi_start: int = 0,
    seed: int = 0,
    max_evaluations: int = 200,
    callback: Callable[[int, float], None] | None = None,
    optimizer: str = "least_squares",
    parallel: int = 1,
) -> FitResult:
    """Estimate parameters from one or more experiments.

    Parameters
    ----------
    experiments : list of FitExperiment
        Each with its project and targets (``SignalTarget``, ``AssayTarget``).
    parameters : list of FitParameter
        Parameters to estimate; start values default to the experiments' projects.
    multi_start : int
        Number of additional random starting points (uniform in the scaled space).
    seed : int
        Random seed for the extra starts.
    max_evaluations : int
        Maximum residual evaluations per start (passed to the optimizer).
    callback : callable, optional
        Called as ``callback(evaluation_number, cost)`` after every evaluation (progress bars).
    optimizer : str
        Registered single-objective optimizer (``least_squares``, ``nelder_mead``, ``lbfgsb``,
        ``differential_evolution`` …).
    parallel : int
        Number of processes for simultaneous simulations (Jacobian columns, populations).

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

    labels, owners, experiment_of = _expand(parameters, names)
    start_values = []
    for parameter, label, owner in zip(owners, labels, experiment_of, strict=True):
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
        start_values.append(value)

    residuals = FitResiduals(experiments, owners, experiment_of)
    residuals(np.array(start_values))  # fixes the residual size for failure penalties

    rng = np.random.default_rng(seed)
    first = np.array([p.to_scaled(v) for p, v in zip(owners, start_values, strict=True)])
    starts = [first] + [rng.uniform(0.0, 1.0, len(first)) for _ in range(multi_start)]
    options = {"max_evaluations": max_evaluations} if optimizer != "differential_evolution" else {}
    history, start_costs, best = [], [], None
    evaluations_so_far = 0
    for start in starts:
        variables = [
            Variable(label, p.lower, p.upper, log_scale=p.scale == "log", initial=p.from_scaled(x))
            for label, p, x in zip(labels, owners, start, strict=True)
        ]
        problem = OptimizationProblem(variables, residuals=residuals, parallel=parallel)

        def progress(n, f, offset=evaluations_so_far):
            if callback is not None:
                callback(offset + n, float(f[0]))

        outcome = get_optimizer(optimizer, **options).optimize(problem, callback=progress)
        evaluations_so_far += outcome.n_evaluations
        history += [(e.time_s, float(e.f[0])) for e in outcome.history]
        start_costs.append((start, float(outcome.f[0])))
        logger.info("Start %s → cost %.4g (%s)", np.round(start, 3), outcome.f[0], outcome.message)
        if best is None or outcome.f[0] < best[1].f[0]:
            best = (problem, outcome)

    assert best is not None  # there is always at least one start
    problem, outcome = best
    final_residuals, projects, results = residuals.evaluate(outcome.x)
    x_scaled = np.array(
        [v.to_scaled(value) for v, value in zip(problem.variables, outcome.x, strict=True)]
    )
    jacobian = outcome.jacobian
    if jacobian is None:  # optimizer without a Jacobian: compute it at the optimum for uncertainty
        scaled = ScaledProblem(problem)
        try:
            jacobian = finite_difference_jacobian(scaled, x_scaled, final_residuals)
        finally:
            scaled.close()  # stops the worker processes of parallel > 1
    return FitResult(
        parameters=list(parameters),
        experiment_names=names,
        values={label: float(value) for label, value in zip(labels, outcome.x, strict=True)},
        x=x_scaled,
        cost=0.5 * float(final_residuals @ final_residuals),
        residuals=final_residuals,
        jacobian=jacobian,
        projects=projects,
        results=results,
        n_evaluations=evaluations_so_far,
        success=outcome.success,
        message=outcome.message,
        history=history,
        starts=start_costs,
    )
