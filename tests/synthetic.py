"""Exact synthetic signals for calibration tests and fixtures (no CADET needed).

The chromatography system is modelled as plug flow (tubing, volume V_t) followed by an ideal
stirred tank (mixer, volume V_m). For an input that is piecewise linear in pumped volume V,
the stirred tank dc_out/dV = (c_in − c_out)/V_m is solved exactly step by step:

    c_out(V_k + h) = a + s h − s V_m + (c_out(V_k) − a + s V_m) exp(−h/V_m)

with a = c_in(V_k) and s the slope of c_in over the step. The first moment then shifts by
exactly V_t + V_m and the variance grows by exactly V_m² (Levenspiel 1999, ch. 11).
"""

import numpy as np


def stirred_tank(volume, c_in, mixer_volume):
    """Outlet of an ideal stirred tank for a piecewise-linear inlet on the grid ``volume``."""
    c_out = np.zeros_like(c_in, dtype=float)
    c_out[0] = c_in[0]
    for k in range(len(volume) - 1):
        h = volume[k + 1] - volume[k]
        a, s = c_in[k], (c_in[k + 1] - c_in[k]) / h
        decay = np.exp(-h / mixer_volume)
        c_out[k + 1] = a + s * h - s * mixer_volume + (c_out[k] - a + s * mixer_volume) * decay
    return c_out


def through_system(volume, c_in, mixer_volume, tubing_volume):
    """Pass a signal through tubing (pure delay) and then the mixer."""
    delayed = np.interp(volume - tubing_volume, volume, c_in, left=c_in[0])
    return stirred_tank(volume, delayed, mixer_volume) if mixer_volume > 0 else delayed


def rectangular_pulse(volume, start, width, height=1.0):
    """Rectangular injection, linear ramps of one grid step at the edges."""
    return height * np.interp(
        volume, [start - 1e-30, start, start + width, start + width + 1e-30], [0.0, 1.0, 1.0, 0.0]
    )
