"""Backend-independent simulation interface: the ``SimulatorBackend`` protocol and its result.

Every backend (CADET, later the reference solver) turns a :class:`~chromsim.core.project.Project`
into a :class:`SimulationResult`. Code that analyses or plots results only uses this module, so it
works the same for every backend.
"""

import datetime
import getpass
from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np

from chromsim import units
from chromsim.core.method import ResolvedStep
from chromsim.core.project import Project
from chromsim.io.project_file import project_to_yaml


@dataclass
class SimulationResult:
    """Outlet concentrations of a simulated run, with time/volume axes and provenance.

    Attributes
    ----------
    component_names : list of str
        Components, in project order (columns of ``outlet_mol_per_m3``).
    time_s : numpy.ndarray
        Time since the start of the method, s. Shape (n_time,).
    outlet_mol_per_m3 : numpy.ndarray
        Concentration at the system outlet, mol/m³ (= mM). Shape (n_time, n_components).
    volume_m3 : numpy.ndarray
        Volume pumped since the start of the method, m³. Shape (n_time,).
    column_volume_m3 : float
        Packed-bed volume (1 CV), m³.
    steps : list of ResolvedStep
        The method on the SI time axis.
    molar_masses_kg_per_mol : list of float or None
        Molar mass per component, kg/mol (``None`` if ``mw_kDa`` is not set).
    solver_stats : dict
        Backend-specific information, e.g. wall time and exit message.
    provenance : dict
        Project YAML, package versions, backend, solver settings, creation time and user
        (SPEC §2.5, §8.1).
    """

    component_names: list[str]
    time_s: np.ndarray
    outlet_mol_per_m3: np.ndarray
    volume_m3: np.ndarray
    column_volume_m3: float
    steps: list[ResolvedStep]
    molar_masses_kg_per_mol: list[float | None]
    solver_stats: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def column_volumes(self) -> np.ndarray:
        """Volume pumped since the start of the method, in column volumes (CV)."""
        return self.volume_m3 / self.column_volume_m3

    def outlet(self, component: str, unit: str = "mM") -> np.ndarray:
        """Outlet concentration of one component.

        Parameters
        ----------
        component : str
            Component name.
        unit : str, optional
            Molar (``"mM"``, default) or mass (``"g/L"``) concentration unit. Mass units need
            the component's ``mw_kDa``.

        Returns
        -------
        numpy.ndarray
            Concentration over time in ``unit``. Shape (n_time,).
        """
        index = self.component_names.index(component)
        return units.concentration_from_si(
            self.outlet_mol_per_m3[:, index], unit, self.molar_masses_kg_per_mol[index]
        )


class SimulatorBackend(Protocol):
    """What every simulation backend provides."""

    name: str

    def simulate(self, project: Project) -> SimulationResult:
        """Simulate the project's method and return the outlet concentrations."""
        ...


def pumped_volume(time_s: np.ndarray, steps: list[ResolvedStep]) -> np.ndarray:
    """Volume pumped since t = 0, m³, for piecewise-constant flow rates.

    Parameters
    ----------
    time_s : numpy.ndarray
        Times, s.
    steps : list of ResolvedStep
        Method steps with start/end times (s) and flow rates (m³/s).

    Returns
    -------
    numpy.ndarray
        V(t) = Σ Q_k · (time spent in step k up to t), m³.
    """
    knot_times = [steps[0].start_s] + [step.end_s for step in steps]
    knot_volumes = [0.0]
    for step in steps:
        knot_volumes.append(knot_volumes[-1] + step.flow_rate_m3_per_s * step.duration_s)
    return np.interp(time_s, knot_times, knot_volumes)


def make_provenance(project: Project, backend_name: str, versions: dict[str, str]) -> dict:
    """Everything needed to reproduce a result (SPEC §2.5).

    Parameters
    ----------
    project : Project
        The simulated project.
    backend_name : str
        Backend used.
    versions : dict
        Package name → version (ChromSim, CADET-Process, CADET-Core …).
    """
    return {
        "project_yaml": project_to_yaml(project),
        "versions": dict(versions),
        "backend": backend_name,
        "solver": project.solver.model_dump(mode="json"),
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "created_by": getpass.getuser(),
    }
