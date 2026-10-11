"""Film mass transfer between bulk liquid and lumped particle pores (LRMP).

.. math:: j = \\frac{3}{r_p} k_f \\,(c - c_p)

is the transfer rate per particle volume, mol/(m³ s). In the LRMP it enters the bulk balance as
−(1 − ε_c)/ε_c · j and the particle balance as j/ε_p (see solver.py).

Reference: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
Chromatography*, 2nd ed. 2006, ch. 6; Leweke & von Lieres, Comput. Chem. Eng. 113 (2018) 274.
"""

import numpy as np


def film_transfer(
    c: np.ndarray, c_p: np.ndarray, film_coefficient: np.ndarray, particle_radius: float
) -> np.ndarray:
    """Transfer rate into the particles, 3 k_f (c − c_p)/r_p.

    Parameters
    ----------
    c, c_p : numpy.ndarray
        Bulk and pore concentrations, shape (n_components, n_cells), mol/m³.
    film_coefficient : numpy.ndarray
        k_f per component, shape (n_components,), m/s.
    particle_radius : float
        r_p, m.

    Returns
    -------
    numpy.ndarray
        mol/(m³ particle · s), shape like ``c``.
    """
    return 3.0 / particle_radius * film_coefficient[:, np.newaxis] * (c - c_p)
