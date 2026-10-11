"""CADET backend: Project → CADET-Process ``Process`` → CADET-Core → :class:`SimulationResult`.

Flow sheet (units with zero volume are left out)::

    inlet → mixer (Cstr) → tubing_in (TubularReactor) → column
          → tubing_out (TubularReactor) → outlet

The tubing volume is split equally before and after the column.

Method → events: every step sets the inlet flow rate and the inlet concentrations at its start
time. Inlet concentrations are first-order polynomials in the time since the step started,
c(t) = c_start + (c_end − c_start)/Δt · t, so a gradient is linear within its step.

Initial state: a clean column equilibrated with the first step's modifiers. Mobile-phase
modifiers the isotherm needs (salt for SMA) start at their first-step inlet value in the column,
mixer and tubing; all other components start at zero and nothing is bound (for SMA all sites carry
salt, q_0 = Λ). Tracer pulses and breakthrough experiments can therefore start at t = 0.

All CADET-Process names used here were read from the installed CADET-Process (0.12) and are pinned
by ``tests/backends/test_cadet_names.py`` (CLAUDE.md rule 4). ``Cadet()`` always gets its own
temporary folder as a ``Path`` (see ADR 0001).

CADET-Process is imported inside the functions that need it, so ``import chromsim`` stays fast and
works for code that never simulates.
"""

import logging
import tempfile
import time
from pathlib import Path

import numpy as np

from chromsim import units
from chromsim.backends.base import SimulationResult, make_provenance, pumped_volume
from chromsim.core.column import ColumnModel
from chromsim.core.project import Project
from chromsim.diagnostics import collect_versions
from chromsim.errors import SimulationError, SimulationSetupError
from chromsim.isotherms.cadet_native import CadetNativeIsotherm
from chromsim.isotherms.registry import get_isotherm

logger = logging.getLogger(__name__)

# Axial cells for each tubing segment. Tubing is short and nearly plug flow; 50 cells resolve
# its residence time well, and its broadening is small next to the mixer and the column.
TUBING_AXIAL_CELLS = 50


def check_simulatable(project: Project) -> None:
    """Check that a project has everything the chosen column model needs.

    Raises
    ------
    SimulationSetupError
        Listing every missing transport parameter.
    """
    column = project.column
    names = project.component_names
    problems = []

    if column.axial_dispersion is None:
        problems.append(
            "column.axial_dispersion is not set; add e.g. `axial_dispersion: 1e-7 m^2/s`"
        )
    if column.model in (ColumnModel.LRMP, ColumnModel.GRM):
        missing = [name for name in names if name not in (column.film_diffusion or {})]
        if missing:
            example = ", ".join(f"{name}: 1e-5 m/s" for name in missing)
            problems.append(
                f"column.film_diffusion is missing for {', '.join(missing)}; "
                f"add e.g. `film_diffusion: {{{example}}}`"
            )
    if column.model is ColumnModel.GRM:
        missing = [name for name in names if name not in (column.pore_diffusion or {})]
        if missing:
            example = ", ".join(f"{name}: 1e-11 m^2/s" for name in missing)
            problems.append(
                f"column.pore_diffusion is missing for {', '.join(missing)}; "
                f"add e.g. `pore_diffusion: {{{example}}}`"
            )

    if problems:
        raise SimulationSetupError(
            f"Project '{project.name}' cannot be simulated with the {column.model} yet:\n"
            + "\n".join(f"  - {problem}" for problem in problems)
        )


def _cadet_isotherm(project: Project) -> CadetNativeIsotherm:
    """The project's isotherm, which must be one CADET implements."""
    isotherm = get_isotherm(project.isotherm.model)
    if not isinstance(isotherm, CadetNativeIsotherm):
        raise SimulationSetupError(
            f"Isotherm '{isotherm.name}' is not built into CADET; simulate it with "
            "`solver: {backend: reference}`."
        )
    return isotherm


def _binding_model(project: Project, component_system):
    """Create and fill the CADET-Process binding model for the project's isotherm."""
    import CADETProcess.processModel as process_model

    isotherm = _cadet_isotherm(project)
    parameters = project.isotherm.parameters
    names = project.component_names
    if not parameters and not isotherm.requires_salt:
        return process_model.NoBinding(component_system, name="binding")

    binding_class = getattr(process_model, isotherm.cadet_class)
    binding = binding_class(component_system, name="binding")
    binding.bound_states = isotherm.bound_states(names, parameters)
    binding.is_kinetic = isotherm.is_kinetic(parameters)
    values = isotherm.cadet_parameters(
        names,
        parameters,
        project.resin.ionic_capacity_mol_per_m3,
        project.isotherm.global_parameters or {},
    )
    for attribute, value in values.items():
        setattr(binding, attribute, value)
    return binding


def _column_unit(project: Project, component_system, c_init: list[float]):
    """Create the column unit (LRM, LRMP or GRM) with geometry, transport and initial state."""
    import CADETProcess.processModel as process_model

    column, resin, names = project.column, project.resin, project.component_names
    unit_classes = {
        ColumnModel.LRM: process_model.LumpedRateModelWithoutPores,
        ColumnModel.LRMP: process_model.LumpedRateModelWithPores,
        ColumnModel.GRM: process_model.GeneralRateModel,
    }
    # Presence of these parameters is guaranteed by check_simulatable().
    assert column.axial_dispersion is not None
    film_diffusion = column.film_diffusion or {}
    pore_diffusion = column.pore_diffusion or {}

    unit = unit_classes[column.model](component_system, name="column")
    unit.length = column.length_m
    unit.diameter = column.diameter_m
    unit.axial_dispersion = units.to_si(column.axial_dispersion, "diffusivity")
    unit.discretization.ncol = project.solver.n_col

    if column.model is ColumnModel.LRM:
        unit.total_porosity = resin.total_porosity
    else:
        unit.bed_porosity = resin.col_porosity
        unit.particle_porosity = resin.par_porosity
        unit.particle_radius = resin.particle_radius_m
        unit.film_diffusion = [units.to_si(film_diffusion[name], "velocity") for name in names]
    if column.model is ColumnModel.GRM:
        unit.pore_diffusion = [units.to_si(pore_diffusion[name], "diffusivity") for name in names]
        unit.discretization.npar = project.solver.n_par

    unit.binding_model = _binding_model(project, component_system)
    unit.c = c_init
    if column.model is not ColumnModel.LRM:
        unit.cp = c_init
    isotherm = _cadet_isotherm(project)
    q_init = isotherm.initial_bound_state(
        names, project.isotherm.parameters, resin.ionic_capacity_mol_per_m3
    )
    if q_init:
        unit.q = q_init
    return unit


def build_process(project: Project):
    """Build the CADET-Process ``Process`` for a project.

    Parameters
    ----------
    project : Project
        A project that passes :func:`check_simulatable`.

    Returns
    -------
    CADETProcess.processModel.Process
        Flow sheet with events for every method step; ready for ``Cadet().simulate``.
    """
    import CADETProcess.processModel as process_model

    check_simulatable(project)
    names = project.component_names
    steps = project.resolved_method()
    # Clean column, equilibrated with the first step's modifiers (e.g. salt); see module docstring.
    modifiers = _cadet_isotherm(project).modifier_components(names)
    c_init = [
        concentration if name in modifiers else 0.0
        for name, concentration in zip(names, steps[0].inlet_start_mol_per_m3, strict=True)
    ]
    system = project.system
    component_system = process_model.ComponentSystem(names)
    flow_sheet = process_model.FlowSheet(component_system)

    inlet = process_model.Inlet(component_system, name="inlet")
    inlet.flow_rate = steps[0].flow_rate_m3_per_s
    chain = [inlet]

    if system.mixer_volume_m3 > 0:
        mixer = process_model.Cstr(component_system, name="mixer")
        mixer.init_liquid_volume = system.mixer_volume_m3
        mixer.c = c_init
        chain.append(mixer)

    tubing_units = []
    if system.tubing_volume_m3 > 0:
        tubing_cross_section = np.pi * (system.tubing_diameter_m / 2.0) ** 2
        for unit_name in ("tubing_in", "tubing_out"):
            tubing = process_model.TubularReactor(component_system, name=unit_name)
            tubing.length = (system.tubing_volume_m3 / 2.0) / tubing_cross_section
            tubing.diameter = system.tubing_diameter_m
            tubing.axial_dispersion = system.tubing_axial_dispersion_m2_per_s
            tubing.discretization.ncol = TUBING_AXIAL_CELLS
            tubing.c = c_init
            tubing_units.append(tubing)
        chain.append(tubing_units[0])

    chain.append(_column_unit(project, component_system, c_init))
    if tubing_units:
        chain.append(tubing_units[1])
    chain.append(process_model.Outlet(component_system, name="outlet"))

    for unit in chain:
        flow_sheet.add_unit(unit)
    for upstream, downstream in zip(chain[:-1], chain[1:], strict=True):
        flow_sheet.add_connection(upstream, downstream)

    process = process_model.Process(flow_sheet, project.name)
    process.cycle_time = steps[-1].end_s
    for index, step in enumerate(steps):
        slopes = [
            (end - start) / step.duration_s
            for start, end in zip(
                step.inlet_start_mol_per_m3, step.inlet_end_mol_per_m3, strict=True
            )
        ]
        concentrations = [
            [start, slope] for start, slope in zip(step.inlet_start_mol_per_m3, slopes, strict=True)
        ]
        process.add_event(
            f"step{index}_{step.name}_c", "flow_sheet.inlet.c", concentrations, step.start_s
        )
        process.add_event(
            f"step{index}_{step.name}_flow",
            "flow_sheet.inlet.flow_rate",
            step.flow_rate_m3_per_s,
            step.start_s,
        )
    return process


class CadetBackend:
    """Simulate projects with CADET-Core through CADET-Process."""

    name = "cadet"

    def simulate(self, project: Project) -> SimulationResult:
        """Simulate the project's method.

        Parameters
        ----------
        project : Project
            Project with all transport parameters for its column model.

        Returns
        -------
        SimulationResult
            Outlet concentrations on a 1 s time grid, with provenance.

        Raises
        ------
        SimulationSetupError
            If the project is incomplete (see :func:`check_simulatable`).
        SimulationError
            If CADET fails.
        """
        from CADETProcess.simulator import Cadet

        process = build_process(project)
        steps = project.resolved_method()

        start = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory(prefix="chromsim-cadet-") as temp_dir:
                simulator = Cadet(temp_dir=Path(temp_dir))
                simulator.time_integrator_parameters.abstol = project.solver.abstol
                cadet_core_version = str(simulator.version)
                results = simulator.simulate(process)
        except Exception as error:  # CADET-Process raises plain exceptions on solver failure
            raise SimulationError(
                f"CADET failed to simulate project '{project.name}': {error}. Try a finer grid "
                "(solver.n_col) or check parameter magnitudes."
            ) from error
        wall_time_s = time.perf_counter() - start

        if results.exit_flag != 0:
            raise SimulationError(
                f"CADET ended with exit flag {results.exit_flag} for project '{project.name}': "
                f"{results.exit_message}"
            )
        logger.info("Simulated '%s' with CADET in %.2f s", project.name, wall_time_s)

        outlet = results.solution.outlet.outlet
        time_s = np.asarray(outlet.time, dtype=float)
        versions = collect_versions() | {"CADET-Core": cadet_core_version}
        return SimulationResult(
            component_names=project.component_names,
            time_s=time_s,
            outlet_mol_per_m3=np.asarray(outlet.solution, dtype=float),
            volume_m3=pumped_volume(time_s, steps),
            column_volume_m3=project.column.volume_m3,
            steps=steps,
            molar_masses_kg_per_mol=[c.molar_mass_kg_per_mol for c in project.components],
            solver_stats={
                "wall_time_s": wall_time_s,
                "exit_flag": int(results.exit_flag),
                "exit_message": str(results.exit_message),
            },
            provenance=make_provenance(project, self.name, versions),
        )
