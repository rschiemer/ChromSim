"""Long-running workspace jobs and their serialisable setups (SPEC §8.1).

Each setup is a pydantic model that names its inputs in the workspace storage (project names,
uploaded file names), so it can be stored, sent to another process, or exposed through a REST API
later. The ``run_*`` functions do the work and store every output; the ``*_task`` functions wrap
them for a :class:`~chromsim.jobs.base.JobRunner`, reporting progress and honouring cancels.

The same ``run_*`` functions can be called directly from a notebook (no job runner needed).
"""

import logging

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from chromsim.analysis.process_optimization import optimize_process
from chromsim.backends import get_backend
from chromsim.data.assays import read_assay
from chromsim.data.importers import read_chromatogram
from chromsim.data.preprocessing import subtract_baseline
from chromsim.errors import ChromSimError
from chromsim.estimation.objectives import AssayTarget, FitExperiment, SignalTarget
from chromsim.estimation.parameters import FitParameter
from chromsim.estimation.uncertainty import parameter_uncertainty
from chromsim.estimation.workflows.inverse_fit import inverse_fit
from chromsim.jobs.base import JobContext
from chromsim.storage.file_storage import FileStorage, stamp
from chromsim.units import Quantity, to_si

logger = logging.getLogger(__name__)


def _ml(value: float) -> float:
    return to_si(Quantity(magnitude=value, unit="mL"), "volume")


# --- simulation -----------------------------------------------------------------------------------


def run_simulation(storage: FileStorage, project_name: str, result_name: str) -> str:
    """Simulate a stored project and store the result under ``result_name``."""
    project = storage.load_project(project_name)
    result = get_backend(project.solver.backend).simulate(project)
    storage.save_simulation_result(result_name, result)
    storage.save_record(
        "simulations",
        result_name,
        {"project": project_name, "solver_stats": result.solver_stats} | stamp(),
    )
    return result_name


def simulation_task(context: JobContext, project_name: str, result_name: str) -> str:
    context.report(0.1, "simulating")
    return run_simulation(context.storage, project_name, result_name)


# --- inverse fit ----------------------------------------------------------------------------------


class AssayInput(BaseModel):
    """An uploaded assay file and the uploaded mapping profile to read it with."""

    model_config = ConfigDict(extra="forbid")

    file: str
    profile: str


class FitSetup(BaseModel):
    """An inverse fit described by names in the workspace (see :func:`run_fit`).

    Attributes
    ----------
    project : str
        Stored project to start from.
    result_name : str
        Outputs: project ``<result_name>-fitted``, simulation ``<result_name>``, record
        ``fits/<result_name>``.
    chromatogram : str, optional
        Uploaded chromatogram; its UV trace is fitted (if ``fit_uv``) and its fractions are used
        for assays without volumes.
    fit_uv : bool
        Include the UV trace as a target.
    uv_baseline_ml : list of (float, float)
        One or two protein-free volume windows (mL) for the UV baseline.
    uv_window_ml : (float, float), optional
        Volume range (mL) of the UV trace to fit.
    assays : list of AssayInput
    parameters : list of FitParameter
    optimizer : str
    max_evaluations : int
    parallel : int
    multi_start : int
    """

    model_config = ConfigDict(extra="forbid")

    project: str
    result_name: str
    chromatogram: str | None = None
    fit_uv: bool = True
    uv_baseline_ml: list[tuple[float, float]] = Field(default_factory=list)
    uv_window_ml: tuple[float, float] | None = None
    assays: list[AssayInput] = Field(default_factory=list)
    parameters: list[FitParameter]
    optimizer: str = "least_squares"
    max_evaluations: int = 100
    parallel: int = 1
    multi_start: int = 0

    def expected_evaluations(self) -> int:
        """Rough number of residual evaluations over all starts (for progress bars)."""
        return self.max_evaluations * (self.multi_start + 1)


def build_fit_experiment(storage: FileStorage, setup: FitSetup) -> FitExperiment:
    """Turn a setup into a FitExperiment with UV and assay targets."""
    project = storage.load_project(setup.project)
    targets: list = []
    fractions = None
    if setup.chromatogram is not None:
        experiment = read_chromatogram(storage.raw_file_path(setup.chromatogram))
        fractions = experiment.fractions or None
        if setup.fit_uv:
            uv = experiment.signal("uv")
            if setup.uv_baseline_ml:
                uv = subtract_baseline(uv, *[(_ml(a), _ml(b)) for a, b in setup.uv_baseline_ml])
            window = (
                None
                if setup.uv_window_ml is None
                else (_ml(setup.uv_window_ml[0]), _ml(setup.uv_window_ml[1]))
            )
            targets.append(SignalTarget(uv, window=window))
    for assay_input in setup.assays:
        assay = read_assay(
            storage.raw_file_path(assay_input.file), storage.raw_file_path(assay_input.profile)
        )
        targets.append(AssayTarget(assay, fractions=fractions))
    if not targets:
        raise ChromSimError("The fit has no data: choose a chromatogram (with UV) and/or assays.")
    return FitExperiment(setup.result_name, project, targets)


def run_fit(storage: FileStorage, setup: FitSetup, callback=None) -> str:
    """Run an inverse fit and store the fitted project, its simulation and a summary record."""
    experiment = build_fit_experiment(storage, setup)
    fit = inverse_fit(
        [experiment],
        setup.parameters,
        optimizer=setup.optimizer,
        max_evaluations=setup.max_evaluations,
        parallel=setup.parallel,
        multi_start=setup.multi_start,
        callback=callback,
    )
    report = parameter_uncertainty(fit)
    fitted_project = fit.projects[setup.result_name]
    storage.save_project(f"{setup.result_name}-fitted", fitted_project)
    storage.save_simulation_result(setup.result_name, fit.results[setup.result_name])
    storage.save_record(
        "fits",
        setup.result_name,
        {
            "setup": setup.model_dump(mode="json"),
            "values": fit.values,
            "cost": fit.cost,
            "success": fit.success,
            "message": fit.message,
            "n_evaluations": fit.n_evaluations,
            "uncertainty": report.table.to_dict(orient="records"),
            "correlation": {
                "labels": list(report.correlation.columns),
                "matrix": np.asarray(report.correlation).tolist(),
            },
            "high_correlations": report.high_correlations,
        }
        | stamp(),
    )
    return setup.result_name


def fit_task(context: JobContext, setup_json: dict) -> str:
    setup = FitSetup.model_validate(setup_json)
    expected = setup.expected_evaluations()

    def progress(evaluation, cost):
        context.report(
            min(evaluation / expected, 0.99), f"evaluation {evaluation}: cost {cost:.4g}"
        )

    context.report(0.0, "preparing")
    return run_fit(context.storage, setup, callback=progress)


# --- process optimization -------------------------------------------------------------------------


class ProcessSetup(BaseModel):
    """A process optimization described by names in the workspace.

    See :func:`run_process_optimization`.

    Outputs: project and simulation ``<result_name>-best``, optimization history
    ``optimizations/<result_name>``, record ``processes/<result_name>`` with the result table.
    """

    model_config = ConfigDict(extra="forbid")

    project: str
    result_name: str
    target: str
    purity: float
    variables: list[FitParameter]
    objectives: list[str] = Field(default_factory=lambda: ["recovery"])
    pool_steps: list[str] | None = None
    optimizer: str | None = None
    options: dict = Field(default_factory=dict)
    parallel: int = 1

    def expected_evaluations(self) -> int:
        """Rough number of simulations (for progress bars)."""
        if "max_evaluations" in self.options:
            return int(self.options["max_evaluations"])
        if "n_gen" in self.options or len(self.objectives) > 1:
            return int(self.options.get("pop_size", 40)) * int(self.options.get("n_gen", 50))
        return 400


def run_process_optimization(storage: FileStorage, setup: ProcessSetup, callback=None) -> str:
    """Optimise the process and store the best project, its simulation and the result table."""
    project = storage.load_project(setup.project)
    outcome = optimize_process(
        project,
        setup.variables,
        setup.target,
        setup.purity,
        objectives=tuple(setup.objectives),
        pool_steps=tuple(setup.pool_steps) if setup.pool_steps else None,
        optimizer=setup.optimizer,
        parallel=setup.parallel,
        callback=callback,
        **setup.options,
    )
    best = outcome.best
    storage.save_optimization_result(setup.result_name, outcome.optimization)
    storage.save_project(f"{setup.result_name}-best", best.project)
    storage.save_simulation_result(f"{setup.result_name}-best", best.result)
    performance = best.performance
    storage.save_record(
        "processes",
        setup.result_name,
        {
            "setup": setup.model_dump(mode="json"),
            "table": outcome.table.to_dict(orient="records"),
            "best_variables": dict(
                zip(
                    outcome.optimization.variable_names,
                    outcome.optimization.x.tolist(),
                    strict=True,
                )
            ),
            "best_pool": None
            if best.pool is None
            else {"start_s": best.pool.start_s, "end_s": best.pool.end_s},
            "best_kpis": None
            if performance is None
            else {
                "recovery": performance.recovery,
                "purity": performance.purity,
                "productivity_kg_per_m3_s": performance.productivity_kg_per_m3_s,
                "buffer_consumption_m3_per_kg": performance.buffer_consumption_m3_per_kg,
            },
        }
        | stamp(),
    )
    return setup.result_name


def process_task(context: JobContext, setup_json: dict) -> str:
    setup = ProcessSetup.model_validate(setup_json)
    expected = setup.expected_evaluations()

    def progress(evaluation, objectives):
        context.report(min(evaluation / expected, 0.99), f"{evaluation} methods simulated")

    context.report(0.0, "preparing")
    return run_process_optimization(context.storage, setup, callback=progress)
