"""Reference-backend operators: conservation, boundary conditions and the TVD property."""

import numpy as np
import pytest

from chromsim.backends.reference.convection import convection
from chromsim.backends.reference.dispersion import dispersion
from chromsim.backends.reference.film_transfer import film_transfer

DZ = 0.01


def test_convection_conserves_mass():
    # Σ dc/dt · Δz = u c_in − u c_out (flux in minus flux out)
    rng = np.random.default_rng(0)
    c = rng.uniform(0.0, 2.0, (2, 50))
    c_in = np.array([1.5, 0.2])
    rate = convection(c, c_in, velocity=1e-3, dz=DZ)
    assert rate.sum(axis=1) * DZ == pytest.approx(1e-3 * (c_in - c[:, -1]))


def test_convection_of_uniform_profile_equal_to_inlet_is_zero():
    c = np.full((1, 20), 3.0)
    assert convection(c, np.array([3.0]), 1e-3, DZ) == pytest.approx(0.0)


def test_convection_is_total_variation_diminishing_for_a_step():
    # A step must not create over- or undershoots (van Leer limiter).
    c = np.where(np.arange(40) < 20, 1.0, 0.0)[np.newaxis, :]
    for _ in range(200):  # explicit Euler at CFL 0.4
        c = c + 0.4 * DZ / 1e-3 * convection(c, np.array([1.0]), 1e-3, DZ)
    assert c.max() <= 1.0 + 1e-12 and c.min() >= -1e-12


def test_dispersion_conserves_mass_and_has_closed_ends():
    rng = np.random.default_rng(1)
    c = rng.uniform(0.0, 2.0, (3, 30))
    assert dispersion(c, 1e-7, DZ).sum(axis=1) == pytest.approx(0.0, abs=1e-18)
    linear = np.linspace(0.0, 1.0, 30)[np.newaxis, :]
    interior = dispersion(linear, 1e-7, DZ)[0, 1:-1]
    assert interior == pytest.approx(0.0, abs=1e-12)  # ∂²c/∂z² = 0 inside


def test_film_transfer_rate():
    # 3/r_p k_f (c − c_p) = 3 / 30e-6 m × 1e-5 m/s × (2 − 0.5) mol/m³ = 1.5 mol/(m³ s)
    rate = film_transfer(np.array([[2.0]]), np.array([[0.5]]), np.array([1e-5]), 30e-6)
    assert rate[0, 0] == pytest.approx(1.5)
