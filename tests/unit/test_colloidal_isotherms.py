"""Colloidal isotherms: the CPA building blocks against textbook limits, and both flux functions
against independent formulas.

Tolerances: identities that hold exactly are checked to round-off (rel 1e-9 to 1e-6); limits
that hold only asymptotically state the size of the neglected term next to the check.
"""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from chromsim.isotherms.custom.colloidal_particle_adsorption import (
    AVOGADRO,
    BOLTZMANN,
    ELEMENTARY_CHARGE,
    VACUUM_PERMITTIVITY,
    available_surface_function,
    debye_kappa,
    grahame_potential,
    henry_factor,
    lateral_energy,
    well_depth,
)
from chromsim.isotherms.python_isotherm import binding_parameters
from chromsim.isotherms.registry import get_isotherm

EPS = 78.4 * VACUUM_PERMITTIVITY
KT = BOLTZMANN * 298.15


# --- CPA building blocks ------------------------------------------------------------------


def test_debye_length_of_100_mM_salt_is_0_96_nm():
    # Israelachvili: 1/κ = 0.304 nm / √(I in M) for 1:1 salt at 25 °C → 0.961 nm at 0.1 M.
    assert 1.0 / debye_kappa(100.0) == pytest.approx(0.961e-9, rel=5e-3)
    assert debye_kappa(400.0) / debye_kappa(100.0) == pytest.approx(2.0, rel=1e-12)


def test_grahame_reduces_to_debye_hueckel_for_small_charge():
    kappa = debye_kappa(100.0)
    sigma = 1e-4  # C/m²: e ψ/k_B T ≈ 0.006, so the cubic correction is ~1e-5 relative
    assert grahame_potential(sigma, kappa) == pytest.approx(sigma / (EPS * kappa), rel=1e-4)
    big = 0.2
    assert abs(grahame_potential(big, kappa)) < big / (EPS * kappa)  # Grahame saturates
    assert grahame_potential(-big, kappa) == -grahame_potential(big, kappa)


def hhf_energy(z, psi_a, psi_i, radius, kappa):
    """Constant-charge sphere–plate energy (Usui 1973), written in z, J."""
    x = np.exp(-kappa * z)
    return (
        np.pi
        * EPS
        * radius
        * (2 * psi_a * psi_i * np.log((1 + x) / (1 - x)) - (psi_a**2 + psi_i**2) * np.log(1 - x**2))
    )


@pytest.mark.parametrize(("psi_a", "psi_i"), [(-0.12, 0.02), (-0.05, 0.04), (0.08, -0.01)])
def test_well_depth_matches_numerical_minimum(psi_a, psi_i):
    radius, kappa = 5e-9, debye_kappa(50.0)
    result = minimize_scalar(
        lambda z: hhf_energy(z, psi_a, psi_i, radius, kappa),
        bounds=(1e-12, 20 / kappa),
        method="bounded",
        options={"xatol": 1e-16},
    )
    assert well_depth(psi_a, psi_i, radius) == pytest.approx(-result.fun / KT, rel=1e-6)


def test_like_charges_have_no_well():
    assert well_depth(-0.1, -0.02, 5e-9) == 0.0
    assert well_depth(0.1, 0.0, 5e-9) == 0.0


def test_henry_factor_limits():
    assert henry_factor(0.0) == pytest.approx(1.0)
    assert henry_factor(1.0) == pytest.approx(np.e - 1.0, rel=1e-12)
    assert np.isfinite(henry_factor(1e4))


def test_available_surface_function_single_component_is_spt():
    radius = np.array([[5e-9]])
    coverage = np.array([0.01, 0.2, 0.5])
    density = coverage[np.newaxis, :] / (np.pi * AVOGADRO * radius**2)
    expected = (1 - coverage) * np.exp(
        -3 * coverage / (1 - coverage) - (coverage / (1 - coverage)) ** 2
    )
    assert available_surface_function(density, radius)[0] == pytest.approx(expected, rel=1e-12)
    # Low coverage: B ≈ 1 − 4Θ (second virial coefficient of hard discs); O(Θ²) ≈ 1e-4 here.
    low = np.array([[1e-3]]) / (np.pi * AVOGADRO * radius**2)
    assert available_surface_function(low, radius)[0, 0] == pytest.approx(1 - 4e-3, abs=2e-5)


def test_available_surface_function_identical_species_add_up():
    radius = np.array([[4e-9], [4e-9]])
    density = np.array([[1e-9], [2e-9]])
    together = available_surface_function(density.sum(axis=0, keepdims=True), radius[:1])
    assert available_surface_function(density, radius) == pytest.approx(
        np.vstack([together, together]), rel=1e-12
    )


def test_available_surface_function_large_proteins_are_excluded_more():
    radius = np.array([[2e-9], [6e-9]])
    density = np.array([[5e-9], [1e-9]])
    small, large = available_surface_function(density, radius)[:, 0]
    assert large < small < 1.0


def test_lateral_energy_sparse_limit_is_six_nearest_neighbours():
    # κD = 8: the lattice-sum factor 1/(1 − e^{−0.83 κD}) differs from 1 by 1.3e-3.
    radius, z_lat = np.array([[4e-9]]), np.array([[6.0]])
    kappa = debye_kappa(100.0)
    spacing = 8.0 / kappa
    density = np.array([[2.0 / (np.sqrt(3.0) * AVOGADRO * spacing**2)]])
    beta = (
        z_lat[0, 0] ** 2
        * ELEMENTARY_CHARGE**2
        / (4 * np.pi * EPS)
        * np.exp(2 * kappa * radius[0, 0])
        / (1 + kappa * radius[0, 0]) ** 2
    )
    six_neighbours = 6 * beta * np.exp(-kappa * spacing) / spacing / KT
    result = lateral_energy(density, radius, z_lat, np.array([kappa]))[0, 0]
    assert result == pytest.approx(six_neighbours, rel=2e-3)


def test_lateral_energy_vanishes_without_charge_or_neighbours():
    radius = np.array([[4e-9]])
    kappa = np.array([debye_kappa(50.0)])
    assert lateral_energy(np.array([[1e-8]]), radius, np.array([[0.0]]), kappa)[0, 0] == 0.0
    assert lateral_energy(np.array([[0.0]]), radius, np.array([[5.0]]), kappa)[0, 0] == 0.0


# --- CPA isotherm ---------------------------------------------------------------------------


def cpa_parameters(charge=30.0, lateral_charge=5.0):
    example = get_isotherm("CPA").example
    values = example["parameters"]["A"] | {"charge": charge, "lateral_charge": lateral_charge}
    return binding_parameters(["salt", "A"], {"A": values}, example["global_parameters"], None)


def test_cpa_henry_coefficient_falls_with_salt():
    cpa = get_isotherm("CPA")
    salt = np.array([10.0, 25.0, 50.0, 100.0, 200.0, 500.0, 1000.0])
    henry = cpa.equilibrium_factor(salt, np.zeros((2, salt.size)), cpa_parameters())[1]
    assert np.all(np.diff(henry) < 0)


def test_cpa_like_charged_protein_only_partitions():
    # A protein with the ligand's charge sign has no well: K = Δ (geometric partition only).
    cpa = get_isotherm("CPA")
    p = cpa_parameters(charge=-10.0)
    factor = cpa.equilibrium_factor(np.array([50.0]), np.zeros((2, 1)), p)[1, 0]
    assert factor == pytest.approx(p["Delta"][1, 0], rel=1e-12)


def test_cpa_isotherm_saturates_below_a_monolayer():
    cpa = get_isotherm("CPA")
    p = cpa_parameters()
    example = cpa.example
    monolayer = example["global_parameters"]["specific_surface"] / (
        np.pi * AVOGADRO * example["parameters"]["A"]["radius"] ** 2
    )  # q at Θ = 1, mol/m³ solid
    q = np.linspace(0.0, 0.75 * monolayer, 400)
    states = np.vstack([np.zeros_like(q), q])
    factor = cpa.equilibrium_factor(np.full_like(q, 50.0), states, p)[1]
    c_eq = q / factor  # liquid concentration in equilibrium with q
    assert np.all(np.diff(c_eq) > 0), "q*(c) must be monotonic"
    # Three quarters of a monolayer need > 1 M protein in the liquid: the isotherm has saturated.
    assert c_eq[-1] > 1e3


def test_cpa_flux_vanishes_at_equilibrium():
    cpa = get_isotherm("CPA")
    p = cpa_parameters()
    q = np.array([[0.0, 0.0], [0.5, 2.0]])
    salt = np.array([50.0, 80.0])
    c = np.vstack([salt, q[1] / cpa.equilibrium_factor(salt, q, p)[1]])
    assert cpa.flux(c, q, p) == pytest.approx(np.zeros_like(q), abs=1e-12)


# --- Colloidal (Xu & Lenhoff) flux ---------------------------------------------------------

COLLOIDAL_A = {
    "radius": 4.5e-9,
    "ln_ke_salt_factor": 40.0,
    "ln_ke_salt_exponent": 0.4,
    "ln_ke_exp_factor": 1.0,
    "ln_ke_exp_salt_coefficient": -1e-3,
    "b_pp_salt_factor": 5.0,
    "b_pp_salt_exponent": 0.2,
    "b_pp_exp_factor": 2.0,
    "b_pp_exp_salt_coefficient": -2e-3,
    "k_kin": 1.5,
}
COLLOIDAL_GLOBAL = {
    "phase_ratio": 5e7,
    "kappa_factor": 1e-10,
    "kappa_exponent": 0.3,
    "kappa_constant": 5e-10,
    "coordination_number": 6,
}


def colloidal_by_hand(c_s, c, q, v, g):
    """One component, written as in the CADET documentation (not as in the implementation)."""
    ln_ke = v["ln_ke_salt_factor"] * c_s ** (-v["ln_ke_salt_exponent"]) + v[
        "ln_ke_exp_factor"
    ] * np.exp(v["ln_ke_exp_salt_coefficient"] * c_s)
    b_pp = v["b_pp_salt_factor"] * c_s ** v["b_pp_salt_exponent"] + v["b_pp_exp_factor"] * np.exp(
        v["b_pp_exp_salt_coefficient"] * c_s
    )
    kappa = 1 / (g["kappa_factor"] * c_s ** g["kappa_exponent"] + g["kappa_constant"])
    big_r = np.sqrt(2 * g["phase_ratio"] / (AVOGADRO * np.sqrt(3) * q))
    r = v["radius"]
    n = g["coordination_number"]
    exponent = (
        n * (3 + kappa * big_r) / (4 * q)
        * q * b_pp * (2 * r) / (2 * big_r) * np.exp(-kappa * (big_r - 2 * r))
        - ln_ke
    )  # fmt: skip
    return v["k_kin"] * (c - q * np.exp(exponent))


@pytest.mark.parametrize(("c_s", "c", "q"), [(50.0, 0.02, 0.3), (300.0, 0.5, 1.0)])
def test_colloidal_flux_matches_documented_equation(c_s, c, q):
    colloidal = get_isotherm("Colloidal")
    p = binding_parameters(["salt", "A"], {"A": COLLOIDAL_A}, COLLOIDAL_GLOBAL, None)
    rate = colloidal.flux(np.array([[c_s], [c]]), np.array([[0.0], [q]]), p)
    assert rate[0, 0] == 0.0
    assert rate[1, 0] == pytest.approx(
        colloidal_by_hand(c_s, c, q, COLLOIDAL_A, COLLOIDAL_GLOBAL), rel=1e-10
    )


def test_colloidal_below_threshold_is_linear():
    colloidal = get_isotherm("Colloidal")
    p = binding_parameters(
        ["salt", "A"], {"A": COLLOIDAL_A}, COLLOIDAL_GLOBAL | {"linear_threshold": 1e-3}, None
    )
    c_s, c, q = 100.0, 1e-3, 5e-4
    ln_ke = 40.0 * c_s**-0.4 + np.exp(-1e-3 * c_s)
    expected = 1.5 * (c - q * np.exp(-ln_ke))
    rate = colloidal.flux(np.array([[c_s], [c]]), np.array([[0.0], [q]]), p)
    assert rate[1, 0] == pytest.approx(expected, rel=1e-12)


def test_colloidal_rapid_equilibrium_passes_unit_rate_to_cadet():
    colloidal = get_isotherm("Colloidal")
    values = {k: v for k, v in COLLOIDAL_A.items() if k != "k_kin"}
    assert not colloidal.is_kinetic({"A": values})
    mapped = colloidal.cadet_parameters(["salt", "A"], {"A": values}, None, COLLOIDAL_GLOBAL)
    assert mapped["kinetic_rate_constant"] == [0.0, 1.0]
    assert mapped["logkeq_power_factor"] == [0.0, 40.0]
    assert mapped["logkeq_power_exponent"] == [0.0, 0.4]


# --- SAS (self-association SMA) -------------------------------------------------------------


def sas_equilibrium_q(c, salt, values, capacity=500.0):
    """Single-component SAS equilibrium q(c) by bisection on q = K (1 + k_sa c) (q̄0/c_s)^ν c."""
    low, high = 0.0, capacity / (values["nu"] + values["sigma"])
    for _ in range(200):
        q = 0.5 * (low + high)
        free = capacity - (values["nu"] + values["sigma"]) * q
        target = (
            values["k_eq"]
            * (1 + values["self_association"] * c)
            * (free / salt) ** values["nu"]
            * c
        )
        low, high = (q, high) if target > q else (low, q)
    return q


def test_sas_flux_vanishes_at_its_equilibrium():
    sas = get_isotherm("SAS")
    values = {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "self_association": 10.0, "k_kin": 10.0}
    p = binding_parameters(["salt", "A"], {"A": values}, None, 500.0)
    c, salt = 0.2, 80.0
    q = sas_equilibrium_q(c, salt, values)
    bound_salt = 500.0 - values["nu"] * q
    rate = sas.flux(np.array([[salt], [c]]), np.array([[bound_salt], [q]]), p)
    assert abs(rate[1, 0]) < 1e-9 * values["k_kin"] * q


def test_sas_is_anti_langmuirian_at_low_loading_and_langmuirian_near_saturation():
    # q/c rises with c while (1 + k_sa c) wins, then falls when the steric term takes over.
    # Linearised at q → 0: d ln(q/c)/dc = k_sa − ν(ν+σ) K_H / Λ, with K_H = k_eq (Λ/c_s)^ν.
    # At 150 mM: ν(ν+σ) K_H / Λ = 125 · 4.1 / 500 ≈ 1.0 < k_sa = 10 → anti-Langmuirian start.
    values = {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "self_association": 10.0}
    c = np.geomspace(1e-4, 5.0, 60)
    affinity = np.array([sas_equilibrium_q(ci, 150.0, values) / ci for ci in c])
    peak = np.argmax(affinity)
    assert 0 < peak < c.size - 1
    assert np.all(np.diff(affinity[: peak + 1]) > 0) and np.all(np.diff(affinity[peak:]) < 0)
    # Without self-association (SMA) q/c only falls: Langmuirian everywhere.
    sma = values | {"self_association": 0.0}
    sma_affinity = np.array([sas_equilibrium_q(ci, 150.0, sma) / ci for ci in c])
    assert np.all(np.diff(sma_affinity) < 0)
