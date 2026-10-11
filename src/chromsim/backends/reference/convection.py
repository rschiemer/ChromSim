"""Convection −u ∂c/∂z, finite volumes with a second-order van Leer (TVD) reconstruction.

Cell i has average c_i; the concentration on its downstream face is reconstructed as

.. math:: c_{i+1/2} = c_i + \\tfrac12 s_i, \\qquad
          s_i = \\frac{\\Delta^-|\\Delta^+| + |\\Delta^-|\\Delta^+}{|\\Delta^-| + |\\Delta^+|},
          \\quad \\Delta^- = c_i - c_{i-1}, \\quad \\Delta^+ = c_{i+1} - c_i

(the van Leer limiter: the harmonic mean of the two slopes, zero at extrema, so no new
oscillations). The inlet face carries u·c_in (Danckwerts: the dispersive part is zero there,
see dispersion.py); the outlet face is first order, c_{N−1/2} = c_{N−1}.

Reference: van Leer, J. Comput. Phys. 14 (1974) 361; LeVeque, *Finite Volume Methods for
Hyperbolic Problems*, Cambridge University Press 2002, ch. 6.
"""

import numpy as np

# Avoids 0/0 in the limiter where both slopes vanish (then the limited slope is 0 anyway).
_TINY = 1e-30


def convection(c: np.ndarray, c_in: np.ndarray, velocity: float, dz: float) -> np.ndarray:
    """Rate of change of the cell averages due to convection, −∂(u c)/∂z.

    Parameters
    ----------
    c : numpy.ndarray
        Cell averages, shape (n_components, n_cells), mol/m³.
    c_in : numpy.ndarray
        Inlet concentrations, shape (n_components,), mol/m³.
    velocity : float
        Interstitial velocity u > 0, m/s.
    dz : float
        Cell length, m.

    Returns
    -------
    numpy.ndarray
        dc/dt from convection, shape like ``c``, mol/(m³ s).
    """
    upstream = np.concatenate([c_in[:, np.newaxis], c[:, :-1]], axis=1)  # c_{i-1}, ghost = c_in
    downstream = np.concatenate([c[:, 1:], c[:, -1:]], axis=1)  # c_{i+1}, zero slope at outlet
    backward, forward = c - upstream, downstream - c
    slope = (backward * np.abs(forward) + np.abs(backward) * forward) / (
        np.abs(backward) + np.abs(forward) + _TINY
    )
    face = c + 0.5 * slope  # value on the downstream face of each cell
    face[:, -1] = c[:, -1]  # outlet face: first order
    flux_out = velocity * face
    flux_in = np.concatenate([velocity * c_in[:, np.newaxis], flux_out[:, :-1]], axis=1)
    return (flux_in - flux_out) / dz
