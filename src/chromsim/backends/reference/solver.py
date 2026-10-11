"""Assemble and integrate the reference model: inlet → mixer → tubing → column → tubing → outlet.

Same flow sheet, initial state and method handling as the CADET backend (cadet_backend.py), so the
two can be compared directly. Column models, with cell averages c (bulk), c_p (pores, LRMP) and q
(bound, per solid volume):

LRM, u = Q / (A ε_t), F_t = (1 − ε_t)/ε_t:

.. math:: \\frac{\\partial c}{\\partial t} = -u\\frac{\\partial c}{\\partial z}
          + D_{ax}\\frac{\\partial^2 c}{\\partial z^2} - F_t \\frac{\\partial q}{\\partial t}

LRMP, u = Q / (A ε_c), F_c = (1 − ε_c)/ε_c, F_p = (1 − ε_p)/ε_p, j = 3 k_f (c − c_p)/r_p:

.. math:: \\frac{\\partial c}{\\partial t} = -u\\frac{\\partial c}{\\partial z}
          + D_{ax}\\frac{\\partial^2 c}{\\partial z^2} - F_c\\, j

.. math:: \\frac{\\partial c_p}{\\partial t} = \\frac{j}{\\varepsilon_p}
          - F_p \\frac{\\partial q}{\\partial t}

Mixer (ideal stirred tank of volume V_m): dc/dt = Q (c_in − c)/V_m. Tubing: dispersed plug flow
(LRM without binding, ε = 1).

Each method step is integrated separately with scipy's BDF method and a sparse Jacobian pattern
(the inlet changes discontinuously between steps).

Reference: Leweke & von Lieres, Comput. Chem. Eng. 113 (2018) 274 (model equations as in CADET);
Shampine & Reichelt, SIAM J. Sci. Comput. 18 (1997) 1 (BDF in scipy).
"""

import logging
import time
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

from chromsim import units
from chromsim.backends.base import SimulationResult, make_provenance, pumped_volume
from chromsim.backends.cadet_backend import TUBING_AXIAL_CELLS, check_simulatable
from chromsim.backends.reference.binding import binding_rate, check_reference_capable
from chromsim.backends.reference.convection import convection
from chromsim.backends.reference.dispersion import dispersion
from chromsim.backends.reference.film_transfer import film_transfer
from chromsim.core.column import ColumnModel
from chromsim.core.project import Project
from chromsim.diagnostics import collect_versions
from chromsim.errors import SimulationError, SimulationSetupError
from chromsim.isotherms.python_isotherm import binding_parameters
from chromsim.isotherms.registry import get_isotherm

logger = logging.getLogger(__name__)

RELATIVE_TOLERANCE = 1e-6  # same as CADET's default


@dataclass
class _Field:
    """Where one concentration field (n_components × n_cells) lives in the state vector."""

    start: int
    n_cells: int
    n_components: int

    @property
    def size(self) -> int:
        return self.n_cells * self.n_components

    def get(self, y: np.ndarray) -> np.ndarray:
        return y[self.start : self.start + self.size].reshape(self.n_components, self.n_cells)

    def index(self, component: int, cell: int | np.ndarray) -> int | np.ndarray:
        return self.start + component * self.n_cells + cell


class _Model:
    """State layout, right-hand side and Jacobian pattern for one project."""

    def __init__(self, project: Project):
        check_simulatable(project)
        column = project.column
        if column.model not in (ColumnModel.LRM, ColumnModel.LRMP):
            raise SimulationSetupError(
                f"The reference backend supports LRM and LRMP columns, not {column.model}; use "
                "`column: {model: LRMP}` or the CADET backend."
            )
        self.project = project
        self.isotherm = get_isotherm(project.isotherm.model)
        names = project.component_names
        self.n_comp = len(names)
        self.binds = bool(project.isotherm.parameters) or self.isotherm.requires_salt
        if self.binds:
            check_reference_capable(
                self.isotherm, self.isotherm.is_kinetic(project.isotherm.parameters)
            )
        self.parameters = binding_parameters(
            names,
            project.isotherm.parameters,
            project.isotherm.global_parameters,
            project.resin.ionic_capacity_mol_per_m3,
        )

        system, resin = project.system, project.resin
        self.lrmp = column.model is ColumnModel.LRMP
        self.n_col = project.solver.n_col
        self.dz = column.length_m / self.n_col
        self.area = column.cross_section_m2
        self.porosity = resin.col_porosity if self.lrmp else resin.total_porosity
        assert column.axial_dispersion is not None  # guaranteed by check_simulatable()
        self.axial_dispersion = units.to_si(column.axial_dispersion, "diffusivity")
        if self.lrmp:
            film = column.film_diffusion or {}
            self.film = np.array([units.to_si(film[name], "velocity") for name in names])
            self.particle_radius = resin.particle_radius_m
            self.par_porosity = resin.par_porosity

        self.mixer_volume = system.mixer_volume_m3
        self.tubing_area = np.pi * (system.tubing_diameter_m / 2.0) ** 2
        self.tubing_dz = (system.tubing_volume_m3 / 2.0) / self.tubing_area / TUBING_AXIAL_CELLS
        self.tubing_dispersion = system.tubing_axial_dispersion_m2_per_s

        # State layout.
        start = 0
        self.mixer = None
        if self.mixer_volume > 0:
            self.mixer = _Field(start, 1, self.n_comp)
            start += self.mixer.size
        self.tubing_in = self.tubing_out = None
        if system.tubing_volume_m3 > 0:
            self.tubing_in = _Field(start, TUBING_AXIAL_CELLS, self.n_comp)
            start += self.tubing_in.size
        self.c = _Field(start, self.n_col, self.n_comp)
        start += self.c.size
        self.cp = None
        if self.lrmp:
            self.cp = _Field(start, self.n_col, self.n_comp)
            start += self.cp.size
        # Bound phase only for components with a bound state: a state that influences nothing
        # makes the finite-difference Jacobian enlarge its step until it overflows.
        bound_states = (
            self.isotherm.bound_states(names, project.isotherm.parameters)
            if self.binds
            else [0] * self.n_comp
        )
        self.bound = np.array([n > 0 for n in bound_states], dtype=bool)
        self.q = _Field(start, self.n_col, int(self.bound.sum()))
        start += self.q.size
        if system.tubing_volume_m3 > 0:
            self.tubing_out = _Field(start, TUBING_AXIAL_CELLS, self.n_comp)
            start += self.tubing_out.size
        self.size = start

    # --- initial state ----------------------------------------------------------------------

    def initial_state(self, first_inlet: np.ndarray) -> np.ndarray:
        """Clean column, equilibrated with the first step's modifiers (as in the CADET backend)."""
        names = self.project.component_names
        modifiers = self.isotherm.modifier_components(names)
        c0 = np.array(
            [first_inlet[i] if names[i] in modifiers else 0.0 for i in range(self.n_comp)]
        )
        y = np.zeros(self.size)
        for field in (self.mixer, self.tubing_in, self.c, self.cp, self.tubing_out):
            if field is not None:
                field.get(y)[:] = c0[:, np.newaxis]
        if self.binds:
            initial = self.isotherm.initial_bound_state(
                names, self.project.isotherm.parameters, self.parameters.ionic_capacity
            )
            self.q.get(y)[:] = np.array(initial, dtype=float)[:, np.newaxis]
        return y

    def _binding(self, liquid: np.ndarray, q_bound: np.ndarray) -> np.ndarray:
        """dq/dt for every component (zero rows where nothing binds) from the bound-state rows."""
        q = np.zeros_like(liquid)
        q[self.bound] = q_bound
        if not self.binds:
            return q
        return binding_rate(self.isotherm, liquid, q, self.parameters)

    # --- right-hand side --------------------------------------------------------------------

    def rhs(
        self, t: float, y: np.ndarray, inlet_start, inlet_slope, step_start: float, flow: float
    ) -> np.ndarray:
        dydt = np.zeros_like(y)
        c_in = inlet_start + inlet_slope * (t - step_start)

        if self.mixer is not None:
            c_m = self.mixer.get(y)[:, 0]
            self.mixer.get(dydt)[:, 0] = flow * (c_in - c_m) / self.mixer_volume
            c_in = c_m
        if self.tubing_in is not None:
            c_in = self._tubing(y, dydt, self.tubing_in, c_in, flow)

        c = self.c.get(y)
        velocity = flow / (self.area * self.porosity)
        dc = convection(c, c_in, velocity, self.dz) + dispersion(c, self.axial_dispersion, self.dz)
        if self.cp is not None:  # LRMP
            c_p = self.cp.get(y)
            dq = self._binding(c_p, self.q.get(y))
            transfer = film_transfer(c, c_p, self.film, self.particle_radius)
            dc -= (1.0 - self.porosity) / self.porosity * transfer
            self.cp.get(dydt)[:] = (
                transfer / self.par_porosity - (1.0 - self.par_porosity) / self.par_porosity * dq
            )
        else:
            dq = self._binding(c, self.q.get(y))
            dc -= (1.0 - self.porosity) / self.porosity * dq
        self.c.get(dydt)[:] = dc
        self.q.get(dydt)[:] = dq[self.bound]

        if self.tubing_out is not None:
            self._tubing(y, dydt, self.tubing_out, c[:, -1], flow)
        return dydt

    def _tubing(self, y, dydt, field: _Field, c_in, flow) -> np.ndarray:
        c = field.get(y)
        velocity = flow / self.tubing_area
        field.get(dydt)[:] = convection(c, c_in, velocity, self.tubing_dz) + dispersion(
            c, self.tubing_dispersion, self.tubing_dz
        )
        return c[:, -1]

    def outlet(self, y: np.ndarray) -> np.ndarray:
        """Concentration leaving the system (last tubing cell, else last column cell)."""
        last = self.tubing_out if self.tubing_out is not None else self.c
        return last.get(y)[:, -1]

    # --- Jacobian sparsity ------------------------------------------------------------------

    def sparsity(self):
        """Which unknowns each equation depends on (for the sparse finite-difference Jacobian)."""
        pattern = lil_matrix((self.size, self.size), dtype=bool)
        transport_fields = [f for f in (self.tubing_in, self.c, self.tubing_out) if f is not None]
        for field in transport_fields:
            for k in range(self.n_comp):
                for cell in range(field.n_cells):
                    for offset in (-2, -1, 0, 1):  # MUSCL stencil + dispersion
                        neighbour = cell + offset
                        if 0 <= neighbour < field.n_cells:
                            pattern[field.index(k, cell), field.index(k, neighbour)] = True
        # Inlet coupling: the first two cells of a field see the upstream outlet.
        chain = [f for f in (self.mixer, self.tubing_in, self.c) if f is not None]
        for upstream, downstream in zip(chain[:-1], chain[1:], strict=True):
            for k in range(self.n_comp):
                for cell in (0, 1):
                    if cell < downstream.n_cells:
                        pattern[
                            downstream.index(k, cell), upstream.index(k, upstream.n_cells - 1)
                        ] = True
        if self.tubing_out is not None:
            for k in range(self.n_comp):
                for cell in (0, 1):
                    for column_cell in (self.n_col - 1,):
                        pattern[self.tubing_out.index(k, cell), self.c.index(k, column_cell)] = True
        if self.mixer is not None:
            for k in range(self.n_comp):
                pattern[self.mixer.index(k, 0), self.mixer.index(k, 0)] = True
        # Binding and film transfer couple all components within one cell.
        liquid = self.cp if self.lrmp else self.c
        for cell in range(self.n_col):
            local_liquid = [liquid.index(k, cell) for k in range(self.n_comp)]
            local_q = [self.q.index(k, cell) for k in range(self.q.n_components)]
            for row in local_liquid + local_q:
                for column in local_liquid + local_q:
                    pattern[row, column] = True
            if self.lrmp:
                for k in range(self.n_comp):
                    pattern[self.c.index(k, cell), self.cp.index(k, cell)] = True
                    pattern[self.cp.index(k, cell), self.c.index(k, cell)] = True
        return pattern.tocsr()


class ReferenceBackend:
    """Simulate LRM/LRMP projects with the pure-Python reference model."""

    name = "reference"

    def simulate(self, project: Project) -> SimulationResult:
        """Simulate the project's method with the reference model.

        Raises
        ------
        SimulationSetupError
            For GRM columns, isotherms without a Python flux, rapid equilibrium, or missing
            transport parameters.
        SimulationError
            If the integrator fails.
        """
        model = _Model(project)
        steps = project.resolved_method()
        y = model.initial_state(np.array(steps[0].inlet_start_mol_per_m3))
        sparsity = model.sparsity()

        times, outlets = [], []
        start_time = time.perf_counter()
        n_evaluations = 0
        for index, step in enumerate(steps):
            inlet_start = np.array(step.inlet_start_mol_per_m3)
            inlet_slope = (np.array(step.inlet_end_mol_per_m3) - inlet_start) / step.duration_s
            first = np.ceil(step.start_s) if index == 0 else np.floor(step.start_s) + 1
            t_eval = np.arange(first, np.floor(step.end_s) + 1)
            t_eval = t_eval[(t_eval >= step.start_s) & (t_eval <= step.end_s)]
            solution = solve_ivp(
                model.rhs,
                (step.start_s, step.end_s),
                y,
                method="BDF",
                t_eval=t_eval,
                args=(inlet_start, inlet_slope, step.start_s, step.flow_rate_m3_per_s),
                jac_sparsity=sparsity,
                rtol=RELATIVE_TOLERANCE,
                atol=project.solver.abstol,
            )
            if not solution.success:
                raise SimulationError(
                    f"Reference backend failed in step '{step.name}': {solution.message}. Try a "
                    "finer grid (solver.n_col) or smaller rate constants."
                )
            n_evaluations += solution.nfev
            y = solution.y[:, -1]
            times.append(solution.t)
            outlets.append(np.array([model.outlet(column) for column in solution.y.T]))
        wall_time_s = time.perf_counter() - start_time
        logger.info("Reference backend simulated '%s' in %.2f s", project.name, wall_time_s)

        time_s = np.concatenate(times)
        return SimulationResult(
            component_names=project.component_names,
            time_s=time_s,
            outlet_mol_per_m3=np.concatenate(outlets),
            volume_m3=pumped_volume(time_s, steps),
            column_volume_m3=project.column.volume_m3,
            steps=steps,
            molar_masses_kg_per_mol=[c.molar_mass_kg_per_mol for c in project.components],
            solver_stats={"wall_time_s": wall_time_s, "rhs_evaluations": int(n_evaluations)},
            provenance=make_provenance(project, self.name, collect_versions()),
        )
