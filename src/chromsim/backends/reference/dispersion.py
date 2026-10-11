"""Axial dispersion D_ax ∂²c/∂z², finite volumes with closed (Danckwerts) boundaries.

Interior faces carry the diffusive flux −D_ax (c_{i+1} − c_i)/Δz. At the inlet the total flux is
u·c_in (Danckwerts), which convection.py supplies, so the dispersive flux there is zero; at the
outlet ∂c/∂z = 0. Together these conserve mass exactly and give the closed–closed first moment.

Reference: Danckwerts, Chem. Eng. Sci. 2 (1953) 1; Leweke & von Lieres, Comput. Chem. Eng. 113
(2018) 274 (boundary conditions as in CADET).
"""

import numpy as np


def dispersion(c: np.ndarray, axial_dispersion: float, dz: float) -> np.ndarray:
    """Rate of change of the cell averages due to axial dispersion.

    Parameters
    ----------
    c : numpy.ndarray
        Cell averages, shape (n_components, n_cells), mol/m³.
    axial_dispersion : float
        D_ax, m²/s.
    dz : float
        Cell length, m.

    Returns
    -------
    numpy.ndarray
        dc/dt from dispersion, shape like ``c``, mol/(m³ s).
    """
    face_flux = -axial_dispersion * np.diff(c, axis=1) / dz  # interior faces, n_cells − 1
    zeros = np.zeros((c.shape[0], 1))
    flux_in = np.concatenate([zeros, face_flux], axis=1)
    flux_out = np.concatenate([face_flux, zeros], axis=1)
    return (flux_in - flux_out) / dz
