"""Public API of ChromSim.

This is the only module the CLI and the web UI import from (SPEC §4, §8). Functions here take and
return serialisable objects (dataclasses, pydantic models, file paths) so they can later be exposed
as a REST API without redesign (SPEC §8.1).
"""

from chromsim import codegen
from chromsim.analysis.detectors import predict_conductivity, predict_uv
from chromsim.analysis.fractionation import (
    AssayComparison,
    Pool,
    PoolPerformance,
    compare_with_assay,
    default_impurities,
    fed_amounts,
    fraction_averages,
    optimal_pool,
    pool_performance,
)
from chromsim.analysis.peaks import (
    PeakMoments,
    asymmetry_factor,
    hetp,
    peak_maximum,
    peak_moments,
    plate_number,
)
from chromsim.analysis.process_optimization import (
    ProcessOptimizationResult,
    optimize_process,
    process_variable_candidates,
)
from chromsim.backends import SimulationResult, get_backend
from chromsim.backends.checks import DiscretisationCheck, check_discretisation
from chromsim.core.column import Column, ColumnModel
from chromsim.core.components import Component
from chromsim.core.method import ResolvedStep, Step
from chromsim.core.project import IsothermConfig, Project, SolverSettings
from chromsim.core.resin import Resin
from chromsim.core.system import System
from chromsim.data.assays import (
    AssayMappingProfile,
    FractionAssay,
    load_mapping_profile,
    read_assay,
)
from chromsim.data.experiment import Experiment, Fraction, Signal
from chromsim.data.importers import list_importers, read_chromatogram
from chromsim.data.preprocessing import (
    ConductivityCalibration,
    conductivity_to_salt,
    load_conductivity_calibration,
    resample,
    subtract_baseline,
    uv_to_concentration,
    volume_to_time,
)
from chromsim.diagnostics import DoctorReport, SmokeTestResult, run_doctor
from chromsim.errors import (
    ChromSimError,
    EnvironmentCheckError,
    IsothermConfigError,
    JobError,
    ProjectFileError,
    SimulationError,
    SimulationSetupError,
    UnitError,
)
from chromsim.estimation.objectives import AssayTarget, FitExperiment, SignalTarget
from chromsim.estimation.parameters import (
    FitParameter,
    fit_parameter_candidates,
    get_value,
    set_values,
)
from chromsim.estimation.uncertainty import UncertaintyReport, parameter_uncertainty
from chromsim.estimation.workflows.inverse_fit import FitResult, inverse_fit
from chromsim.estimation.workflows.tracer import (
    TracerAnalysis,
    TracerRun,
    analyse_tracers,
    tracer_run_from_file,
    tracer_summary,
)
from chromsim.estimation.workflows.yamamoto import (
    LgeRun,
    YamamotoResult,
    elution_salt_at_peak,
    gradient_slope,
    lge_run_from_file,
    yamamoto_fit,
)
from chromsim.io.project_file import (
    load_project,
    project_from_dict,
    project_from_yaml,
    project_to_yaml,
    save_project,
)
from chromsim.io.results import (
    read_optimization_result,
    read_simulation_result,
    write_optimization_result,
    write_simulation_result,
)
from chromsim.isotherms.registry import isotherm_parameter_table, list_isotherms
from chromsim.jobs import (
    AssayInput,
    FitSetup,
    JobCancelled,
    JobStatus,
    LocalProcessRunner,
    ProcessSetup,
    run_fit,
    run_process_optimization,
    run_simulation,
)
from chromsim.optimizers import (
    Constraint,
    OptimizationProblem,
    OptimizationResult,
    Variable,
    get_optimizer,
    list_optimizers,
)
from chromsim.plotting import (
    plot_chromatogram,
    plot_correlation,
    plot_correlation_matrix,
    plot_fraction_comparison,
    plot_overlay,
    plot_pareto,
    plot_yamamoto,
)
from chromsim.storage import FileStorage
from chromsim.units import Quantity, field_dimensions, unit_choices
from chromsim.workspace import Workspace

__all__ = [
    # errors
    "ChromSimError",
    "EnvironmentCheckError",
    "IsothermConfigError",
    "ProjectFileError",
    "SimulationError",
    "SimulationSetupError",
    "UnitError",
    # environment
    "DoctorReport",
    "SmokeTestResult",
    "doctor",
    # domain model
    "Column",
    "ColumnModel",
    "Component",
    "IsothermConfig",
    "Project",
    "Quantity",
    "ResolvedStep",
    "Resin",
    "SolverSettings",
    "Step",
    "System",
    # project files
    "load_project",
    "project_from_dict",
    "project_from_yaml",
    "project_to_yaml",
    "save_project",
    # isotherms
    "list_isotherms",
    # simulation
    "DiscretisationCheck",
    "SimulationResult",
    "check_discretisation",
    "plot_chromatogram",
    "simulate",
    # experimental data
    "Experiment",
    "Fraction",
    "Signal",
    "list_importers",
    "read_chromatogram",
    "AssayMappingProfile",
    "FractionAssay",
    "load_mapping_profile",
    "read_assay",
    "ConductivityCalibration",
    "load_conductivity_calibration",
    "conductivity_to_salt",
    "resample",
    "subtract_baseline",
    "uv_to_concentration",
    "volume_to_time",
    # comparing simulation and experiment
    "AssayComparison",
    "compare_with_assay",
    "fraction_averages",
    "plot_fraction_comparison",
    "plot_overlay",
    "predict_conductivity",
    "predict_uv",
    # peak analysis
    "PeakMoments",
    "asymmetry_factor",
    "hetp",
    "peak_maximum",
    "peak_moments",
    "plate_number",
    # calibration
    "TracerAnalysis",
    "TracerRun",
    "analyse_tracers",
    "LgeRun",
    "YamamotoResult",
    "elution_salt_at_peak",
    "gradient_slope",
    "yamamoto_fit",
    "plot_yamamoto",
    "AssayTarget",
    "FitExperiment",
    "FitParameter",
    "FitResult",
    "SignalTarget",
    "get_value",
    "inverse_fit",
    "set_values",
    "UncertaintyReport",
    "parameter_uncertainty",
    "plot_correlation",
    "plot_correlation_matrix",
    # optimization
    "Constraint",
    "OptimizationProblem",
    "OptimizationResult",
    "Variable",
    "get_optimizer",
    "list_optimizers",
    "optimize",
    # result files
    "read_optimization_result",
    "read_simulation_result",
    "write_optimization_result",
    "write_simulation_result",
    # fractionation and process optimization
    "Pool",
    "PoolPerformance",
    "default_impurities",
    "fed_amounts",
    "optimal_pool",
    "pool_performance",
    "ProcessOptimizationResult",
    "optimize_process",
    "plot_pareto",
    # workspaces, jobs and the UI
    "AssayInput",
    "FileStorage",
    "FitSetup",
    "JobCancelled",
    "JobError",
    "JobStatus",
    "LocalProcessRunner",
    "ProcessSetup",
    "UIUnavailableError",
    "Workspace",
    "codegen",
    "field_dimensions",
    "isotherm_parameter_table",
    "run_fit",
    "run_process_optimization",
    "run_simulation",
    "run_ui",
    "ui_available",
    "unit_choices",
    "fit_parameter_candidates",
    "lge_run_from_file",
    "process_variable_candidates",
    "tracer_run_from_file",
    "tracer_summary",
]


def doctor() -> DoctorReport:
    """Check that ChromSim and CADET are installed and working.

    Collects the versions of Python, ChromSim, CADET-Process, CADET-Python and CADET-Core, and runs
    a sub-second CADET test simulation (tracer pulse through a non-binding column, checked for mass
    balance).

    Returns
    -------
    DoctorReport
        Versions, platform, smoke-test result and a list of problems; ``report.ok`` is ``True``
        if everything works.

    Examples
    --------
    >>> from chromsim import api
    >>> report = api.doctor()
    >>> report.ok
    True
    """
    return run_doctor()


def simulate(project: Project) -> SimulationResult:
    """Simulate a project with the backend named in ``project.solver.backend``.

    Parameters
    ----------
    project : Project
        A complete project (all transport parameters set for its column model).

    Returns
    -------
    SimulationResult
        Outlet concentrations over time and volume, with provenance.

    Raises
    ------
    SimulationSetupError
        If parameters are missing; the message lists all of them.
    SimulationError
        If the solver fails.
    """
    return get_backend(project.solver.backend).simulate(project)


def optimize(
    problem: OptimizationProblem, optimizer: str = "least_squares", callback=None, **options
) -> OptimizationResult:
    """Solve an optimization problem with a registered optimizer (SPEC §6).

    Parameters
    ----------
    problem : OptimizationProblem
        Variables, objectives and/or residuals, constraints, parallelism.
    optimizer : str
        Name from :func:`list_optimizers` (e.g. ``"least_squares"``, ``"nsga2"``).
    callback : callable, optional
        ``callback(evaluation_number, objectives)`` after every evaluation.
    **options
        Optimizer options, e.g. ``n_gen=100`` for ``nsga2``.

    Returns
    -------
    OptimizationResult
        Best point, objectives, full evaluation history (and Pareto set for multi-objective
        optimizers).
    """
    return get_optimizer(optimizer, **options).optimize(problem, callback=callback)


class UIUnavailableError(ChromSimError):
    """The optional web UI is not installed."""


def ui_available() -> bool:
    """``True`` if the optional web UI (NiceGUI) is installed."""
    import importlib.util

    return importlib.util.find_spec("nicegui") is not None


def run_ui(
    workspace_dir="chromsim-workspace",
    port: int = 8080,
    native: bool = False,
    show: bool = True,
    host: str = "127.0.0.1",
) -> None:
    """Start the local web UI (SPEC §8) on a workspace folder.

    The UI has no login (v1), so by default it only accepts connections from this computer
    (``host="127.0.0.1"``). ``host="0.0.0.0"`` makes it reachable from the network: everyone who
    can reach the port can read and change the workspace.

    Raises
    ------
    UIUnavailableError
        If NiceGUI is not installed; the message says how to install it.
    """
    if not ui_available():
        raise UIUnavailableError(
            "The web UI is optional and not installed. Install it with `pip install "
            "'chromsim[ui]'`, or use the pixi default environment (`pixi run ui`). Everything "
            "the UI does is also available from Python (chromsim.api)."
        )
    from chromsim.ui.app import run  # imported only here: chromsim itself never imports NiceGUI

    run(workspace_dir=workspace_dir, port=port, native=native, show=show, host=host)
