"""Numerical checks on simulations, e.g. whether the grid is fine enough (SPEC §7.1)."""

from dataclasses import dataclass

import numpy as np

from chromsim.backends import get_backend
from chromsim.core.project import Project


def nrmse(simulated: np.ndarray, reference: np.ndarray) -> float:
    """Normalised root-mean-square error, dimensionless.

    .. math:: \\mathrm{NRMSE} = \\frac{\\sqrt{\\overline{(y - y_{ref})^2}}}
                                     {\\max y_{ref} - \\min y_{ref}}

    Returns 0 when both signals are identical, also if the reference is constant.
    """
    rmse = float(np.sqrt(np.mean((simulated - reference) ** 2)))
    signal_range = float(np.max(reference) - np.min(reference))
    if rmse == 0.0:
        return 0.0
    if signal_range == 0.0:
        return float("inf")
    return rmse / signal_range


@dataclass
class DiscretisationCheck:
    """Change in the outlet profile when the grid is doubled.

    Attributes
    ----------
    n_col, n_par : int
        Grid of the project.
    refined_n_col, refined_n_par : int
        Doubled grid.
    nrmse : dict[str, float]
        Per component: NRMSE between the outlet on the project grid and on the doubled grid.
    """

    n_col: int
    n_par: int
    refined_n_col: int
    refined_n_par: int
    nrmse: dict[str, float]

    @property
    def max_nrmse(self) -> float:
        """Largest NRMSE over all components."""
        return max(self.nrmse.values())


def check_discretisation(project: Project) -> DiscretisationCheck:
    """Simulate on the project grid and on a doubled grid and compare the outlet profiles.

    An NRMSE well below 1 % means the project grid is fine enough; otherwise increase
    ``solver.n_col`` (and ``n_par`` for the GRM).

    Parameters
    ----------
    project : Project
        Project to check. Not modified.

    Returns
    -------
    DiscretisationCheck
        NRMSE per component between the two grids.
    """
    refined = project.model_copy(deep=True)
    refined.solver.n_col = 2 * project.solver.n_col
    refined.solver.n_par = 2 * project.solver.n_par

    backend = get_backend(project.solver.backend)
    coarse_result = backend.simulate(project)
    fine_result = backend.simulate(refined)
    return DiscretisationCheck(
        n_col=project.solver.n_col,
        n_par=project.solver.n_par,
        refined_n_col=refined.solver.n_col,
        refined_n_par=refined.solver.n_par,
        nrmse={
            name: nrmse(coarse_result.outlet_mol_per_m3[:, i], fine_result.outlet_mol_per_m3[:, i])
            for i, name in enumerate(project.component_names)
        },
    )
