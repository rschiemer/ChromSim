"""Colloidal particle adsorption (CPA) model for ion exchange (Briskot et al. 2021), in Python.

The protein is a charged hard sphere; the adsorber a charged plane. Adsorption follows from
colloid physics instead of a stoichiometric exchange of counter-ions:

1. **Screening.** Ionic strength I sets the Debye parameter κ.
2. **Surface potentials.** Protein charge Z and ligand density Γ_L give surface charge densities;
   the Grahame equation turns them into surface potentials ψ.
3. **Protein–adsorber well.** The double-layer energy between sphere and plane (constant-charge
   Hogg–Healy–Fuerstenau form, Derjaguin approximation) has a minimum u_A(δ_m) < 0 for opposite
   charges. Its depth gives the Henry coefficient.
4. **Crowding.** Hard-disc exclusion on the surface (scaled particle theory: available surface
   function B) and screened repulsion between adsorbed proteins on a hexagonal lattice (u_lat)
   make the isotherm saturate.

Every step is a small function below, so each can be tested and reused (e.g. to plot K_H(I)).
Units: SI throughout; concentrations in mol/m³ (= mM).
"""

import numpy as np

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.python_isotherm import BindingParameters, PythonIsotherm
from chromsim.isotherms.registry import register_isotherm

ELEMENTARY_CHARGE = 1.602176634e-19  # C
AVOGADRO = 6.02214076e23  # 1/mol
BOLTZMANN = 1.380649e-23  # J/K
VACUUM_PERMITTIVITY = 8.8541878128e-12  # F/m

DEFAULT_TEMPERATURE = 298.15  # K
DEFAULT_RELATIVE_PERMITTIVITY = 78.4  # water at 25 °C
MAX_WELL_DEPTH = 700.0  # in k_B T; keeps exp() finite (e^700 ≈ 1e304)
MAX_COVERAGE = 1.0 - 1e-9  # B(Θ) → 0 as Θ → 1; keep the logarithm finite


def debye_kappa(
    ionic_strength, temperature=DEFAULT_TEMPERATURE, eps_r=DEFAULT_RELATIVE_PERMITTIVITY
):
    """Debye parameter κ (inverse Debye length).

    .. math:: \\kappa = \\sqrt{\\frac{2 e^2 N_A I}{\\varepsilon_r \\varepsilon_0 k_B T}}

    Parameters
    ----------
    ionic_strength : float or numpy.ndarray
        I, mol/m³ (= mM).
    temperature : float
        T, K.
    eps_r : float
        Relative permittivity of the solvent.

    Returns
    -------
    float or numpy.ndarray
        κ, 1/m. At 100 mM, 25 °C: 1/κ = 0.96 nm.

    Reference: Israelachvili, *Intermolecular and Surface Forces*, 3rd ed. (2011), eq. 14.36;
    Briskot et al., J. Chromatogr. A 1653 (2021) 462397, eq. 13.
    """
    return np.sqrt(
        2.0
        * ELEMENTARY_CHARGE**2
        * AVOGADRO
        * np.asarray(ionic_strength)
        / (eps_r * VACUUM_PERMITTIVITY * BOLTZMANN * temperature)
    )


def grahame_potential(
    charge_density, kappa, temperature=DEFAULT_TEMPERATURE, eps_r=DEFAULT_RELATIVE_PERMITTIVITY
):
    """Surface potential ψ from surface charge density σ (Grahame equation, 1:1 electrolyte).

    .. math:: \\psi = \\frac{2 k_B T}{e}\\,
              \\operatorname{arsinh}\\left(\\frac{e\\,\\sigma}{2 \\varepsilon_r\\varepsilon_0
              \\kappa k_B T}\\right)

    For small σ this is the Debye–Hückel result ψ = σ/(ε_r ε_0 κ).

    Parameters
    ----------
    charge_density : float or numpy.ndarray
        σ, C/m².
    kappa : float or numpy.ndarray
        κ, 1/m.
    temperature, eps_r : float
        T in K; relative permittivity.

    Returns
    -------
    float or numpy.ndarray
        ψ, V.

    Reference: Grahame, Chem. Rev. 41 (1947) 441; Briskot et al. (2021), eq. 12.
    """
    thermal_voltage = BOLTZMANN * temperature / ELEMENTARY_CHARGE  # k_B T / e, V
    eps = eps_r * VACUUM_PERMITTIVITY
    return (
        2.0 * thermal_voltage * np.arcsinh(charge_density / (2.0 * eps * kappa * thermal_voltage))
    )


def well_depth(
    psi_adsorber,
    psi_protein,
    radius,
    temperature=DEFAULT_TEMPERATURE,
    eps_r=DEFAULT_RELATIVE_PERMITTIVITY,
):
    """Depth of the protein–adsorber energy minimum, y = −u_A(δ_m)/(k_B T) ≥ 0.

    Sphere (radius a) and plane at constant surface charge, linearised Poisson–Boltzmann,
    Derjaguin approximation, with x = e^{−κz}:

    .. math:: u_A(z) = \\pi \\varepsilon_r\\varepsilon_0 a \\left[2\\psi_A\\psi_i
              \\ln\\frac{1+x}{1-x} - (\\psi_A^2 + \\psi_i^2)\\ln(1 - x^2)\\right]

    For opposite charges the minimum lies at x* = −2ψ_Aψ_i/(ψ_A² + ψ_i²) (from du/dx = 0), i.e.
    at δ_m = −ln(x*)/κ. Like charges give no well (y = 0). The depth does not depend on κ
    except through the potentials.

    Parameters
    ----------
    psi_adsorber, psi_protein : float or numpy.ndarray
        Surface potentials ψ_A and ψ_i, V.
    radius : float or numpy.ndarray
        Protein radius a, m.
    temperature, eps_r : float
        T in K; relative permittivity.

    Returns
    -------
    numpy.ndarray
        y = −u_A(δ_m)/(k_B T), dimensionless, ≥ 0.

    Reference: Usui, J. Colloid Interface Sci. 44 (1973) 107 (constant charge); Hogg, Healy &
    Fuerstenau, Trans. Faraday Soc. 62 (1966) 1638; Briskot et al. (2021), eqs. 18–19.
    """
    product = np.asarray(psi_adsorber * psi_protein)
    squares = psi_adsorber**2 + psi_protein**2
    attractive = product < 0.0
    x = np.where(attractive, -2.0 * product / np.where(squares > 0, squares, 1.0), 0.0)
    x = np.clip(x, 0.0, 1.0 - 1e-12)  # |ψ_A| = |ψ_i| puts the minimum at contact
    u_over_prefactor = 2.0 * product * np.log((1.0 + x) / (1.0 - x)) - squares * np.log(1.0 - x**2)
    u = np.pi * eps_r * VACUUM_PERMITTIVITY * radius * u_over_prefactor  # J
    return np.where(attractive, np.maximum(-u / (BOLTZMANN * temperature), 0.0), 0.0)


def henry_factor(y):
    """(e^y − 1)/y, the integral of the idealised linear well over its width (→ 1 for y → 0).

    The Henry coefficient is K_H = Δ (e^y − 1)/y, where Δ = A_s(δ* − δ_m) and y is the well
    depth in k_B T (Briskot et al. 2021, eq. 10, written for u_A < 0).
    """
    y = np.minimum(np.asarray(y, dtype=float), MAX_WELL_DEPTH)
    small = y < 1e-8
    return np.where(small, 1.0 + 0.5 * y, np.expm1(y) / np.where(small, 1.0, y))


def available_surface_function(surface_density, radius):
    """Available surface function B_i(Θ) of a hard-disc mixture (scaled particle theory).

    .. math:: \\ln B_i = \\ln(1-\\Theta) - \\frac{2\\pi a_i \\rho_1 + \\pi a_i^2 \\rho_0}{1-\\Theta}
              - \\frac{(\\pi a_i \\rho_1)^2}{(1-\\Theta)^2}

    with ρ_0 = N_A Σ_j Γ_j, ρ_1 = N_A Σ_j a_j Γ_j and Θ = π N_A Σ_j a_j² Γ_j. For one component
    this is B = (1−Θ) exp[−3Θ/(1−Θ) − Θ²/(1−Θ)²]; at low coverage B ≈ 1 − 4Θ.

    Parameters
    ----------
    surface_density : numpy.ndarray
        Γ_j, mol/m², shape (n_components, n_points); 0 for non-binding components.
    radius : numpy.ndarray
        a_j, m, shape (n_components, 1).

    Returns
    -------
    numpy.ndarray
        B_i, shape (n_components, n_points).

    Reference: Helfand, Frisch & Lebowitz, J. Chem. Phys. 34 (1961) 1037; Talbot et al.,
    Colloids Surf. A 165 (2000) 287; Briskot et al. (2021), eq. 28.
    """
    rho_0 = AVOGADRO * surface_density.sum(axis=0)
    rho_1 = AVOGADRO * (radius * surface_density).sum(axis=0)
    coverage = np.minimum(
        np.pi * AVOGADRO * (radius**2 * surface_density).sum(axis=0), MAX_COVERAGE
    )
    free = 1.0 - coverage
    ln_b = (
        np.log(free)
        - (2.0 * np.pi * radius * rho_1 + np.pi * radius**2 * rho_0) / free
        - (np.pi * radius * rho_1) ** 2 / free**2
    )
    return np.exp(ln_b)


def lateral_energy(
    surface_density,
    radius,
    lateral_charge,
    kappa,
    temperature=DEFAULT_TEMPERATURE,
    eps_r=DEFAULT_RELATIVE_PERMITTIVITY,
):
    """Screened repulsion of protein i with its adsorbed neighbours, u_lat,i/(k_B T).

    Proteins sit on a hexagonal lattice of spacing D (fixed by the total surface density) and
    interact through a Yukawa potential u_ij(r) = β_ij e^{−κr}/r:

    .. math:: \\beta_{ij} = \\frac{Z_{lat,i} Z_{lat,j}\\, e^2}{4\\pi\\varepsilon_r\\varepsilon_0}
              \\frac{e^{\\kappa(a_i+a_j)}}{(1+\\kappa a_i)(1+\\kappa a_j)}, \\qquad
              D^2 = \\frac{2}{\\sqrt{3}\\, N_A \\sum_j \\Gamma_j}

    .. math:: u_{lat,i} = 3\\sqrt{3}\\, D\\, N_A
              \\frac{e^{-\\kappa D}}{1 - \\exp\\left(-\\frac{3\\sqrt{3}}{2\\pi}\\kappa D\\right)}
              \\sum_j \\Gamma_j \\beta_{ij}

    For κD ≫ 1 this is the energy of six nearest neighbours at distance D.

    Parameters
    ----------
    surface_density : numpy.ndarray
        Γ_j, mol/m², shape (n_components, n_points).
    radius, lateral_charge : numpy.ndarray
        a_j (m) and Z_lat,j (dimensionless), shape (n_components, 1).
    kappa : numpy.ndarray
        κ, 1/m, shape (n_points,).
    temperature, eps_r : float
        T in K; relative permittivity.

    Returns
    -------
    numpy.ndarray
        u_lat,i/(k_B T), shape (n_components, n_points).

    Reference: Briskot et al., J. Chromatogr. A 1653 (2021) 462397, eqs. 20–24 (Yukawa
    potential with DLVO prefactor, e.g. Russel, Saville & Schowalter, *Colloidal Dispersions*,
    1989, ch. 4).
    """
    total = surface_density.sum(axis=0)
    occupied = total > 0.0
    spacing = np.sqrt(2.0 / (np.sqrt(3.0) * AVOGADRO * np.where(occupied, total, 1.0)))  # D, m
    bjerrum_energy = ELEMENTARY_CHARGE**2 / (4.0 * np.pi * eps_r * VACUUM_PERMITTIVITY)  # J·m
    screened = lateral_charge * np.exp(kappa * radius) / (1.0 + kappa * radius)  # (n, points)
    beta_sum = screened * (screened * surface_density).sum(axis=0) * bjerrum_energy  # Σ_j Γ_j β_ij
    kd = kappa * spacing
    lattice = 3.0 * np.sqrt(3.0) * spacing * AVOGADRO * np.exp(-kd)
    lattice /= -np.expm1(-3.0 * np.sqrt(3.0) / (2.0 * np.pi) * kd)
    u = np.where(occupied, lattice * beta_sum, 0.0)
    return u / (BOLTZMANN * temperature)


@register_isotherm
class ColloidalParticleAdsorption(PythonIsotherm):
    """Colloidal particle adsorption (CPA) model for ion exchange (Briskot et al. 2021).

    .. math::

        \\frac{dq_i}{dt} = k_{kin,i}\\left(K_{H,i}\\, B_i(\\Theta)\\,
        e^{-u_{lat,i}/k_B T}\\, c_{p,i} - q_i\\right), \\qquad
        K_{H,i} = \\Delta_i \\frac{e^{y_i} - 1}{y_i}, \\qquad
        y_i = -\\frac{u_{A,i}(\\delta_{m,i})}{k_B T}

    y_i is the depth of the protein–adsorber double-layer well (:func:`well_depth`), computed
    from the protein net charge Z_i, the ligand density Γ_L and the ionic strength
    I = c_{p,0} + I_0. B_i is the scaled-particle available surface function of hard discs
    (:func:`available_surface_function`) and u_lat,i the screened repulsion of adsorbed neighbours
    (:func:`lateral_energy`). Surface densities are Γ_j = q_j / A_s. q is per solid volume.

    Simplifications against the paper (see docs/theory/isotherms/cpa.md): strong (fully
    dissociated) ligand, no charge regulation; Z_i and Δ_i independent of pH; constant k_kin;
    one specific surface A_s for all components. Salt is component 0 (1:1 salt).

    References: Briskot, Hahn, Huuk & Hubbuch, J. Chromatogr. A 1653 (2021) 462397; Briskot et
    al., J. Chromatogr. A 1654 (2021) 462439 (kinetic form, eqs. 12–19); Briskot, Hahn, Huuk &
    Hubbuch, J. Chromatogr. A 1611 (2020) 460608.
    """

    name = "CPA"
    requires_salt = True
    n_modifiers = 1
    parameters = {
        # --- per binding component -------------------------------------------------------
        "radius": Parameter(unit="m", description="protein radius a_i (mAb ≈ 5e-9)"),
        "charge": Parameter(
            unit="dimensionless",
            description="net charge Z_i at the operating pH (opposite sign to the ligand binds)",
            bounds=(None, None),
        ),
        "Delta": Parameter(
            unit="dimensionless",
            description="Δ_i = A_s (δ* − δ_m), width of the adsorption layer × surface",
            log_scale=True,
        ),
        "lateral_charge": Parameter(
            unit="dimensionless",
            description="effective charge Z_lat,i for protein–protein repulsion (default 0)",
            required=False,
        ),
        "k_kin": Parameter(unit="1/s", description="desorption rate constant", log_scale=True),
        # --- whole isotherm ----------------------------------------------------------------
        "ligand_density": Parameter(
            unit="mol/m^2",
            description="charged ligands per adsorber surface Γ_L (≈ Λ / A_s)",
            log_scale=True,
            per_component=False,
        ),
        "ligand_charge": Parameter(
            unit="dimensionless",
            description="charge per ligand: −1 cation exchanger, +1 anion exchanger",
            bounds=(-3, 3),
            per_component=False,
        ),
        "specific_surface": Parameter(
            unit="m^2/m^3",
            description="adsorber surface per solid volume A_s",
            log_scale=True,
            per_component=False,
        ),
        "background_ionic_strength": Parameter(
            unit="mol/m^3",
            description="ionic strength of the buffer without salt I_0 (default 0)",
            per_component=False,
            required=False,
        ),
        "temperature": Parameter(
            unit="K",
            description="temperature (default 298.15)",
            per_component=False,
            required=False,
        ),
        "relative_permittivity": Parameter(
            unit="dimensionless",
            description="ε_r of the mobile phase (default 78.4)",
            per_component=False,
            required=False,
        ),
    }
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {
            "A": {"radius": 5e-9, "charge": 30.0, "Delta": 0.5, "lateral_charge": 5.0, "k_kin": 1.0}
        },
        "global_parameters": {
            "ligand_density": 2e-6,
            "ligand_charge": -1.0,
            "specific_surface": 2e8,
        },
        "load": {"salt": "50 mM", "A": "0.05 mM"},
        "elute": {"salt": "1000 mM"},
    }

    def bound_states(self, component_names, parameters):
        """Salt has no bound state."""
        return [0] + [1 if name in parameters else 0 for name in component_names[1:]]

    def modifier_components(self, component_names):
        """Salt (component 0) sets the ionic strength."""
        return [component_names[0]]

    def equilibrium_factor(
        self, c_salt: np.ndarray, q: np.ndarray, p: BindingParameters
    ) -> np.ndarray:
        """K_H,i B_i e^{−u_lat,i/k_BT}, so that q*_i = factor · c_i; shape (n_components, n_points).

        Parameters
        ----------
        c_salt : numpy.ndarray
            Salt concentration c_{p,0}, mol/m³, shape (n_points,).
        q : numpy.ndarray
            Bound concentrations, mol/m³ solid, shape (n_components, n_points).
        p : BindingParameters
            Isotherm parameters.
        """
        g = p.global_
        temperature = g.get("temperature", DEFAULT_TEMPERATURE)
        eps_r = g.get("relative_permittivity", DEFAULT_RELATIVE_PERMITTIVITY)
        binding = p.binding.copy()
        binding[0] = False

        ionic_strength = np.maximum(c_salt, 0.0) + g.get("background_ionic_strength", 0.0)
        kappa = debye_kappa(np.maximum(ionic_strength, 1e-6), temperature, eps_r)
        radius = np.where(binding[:, np.newaxis], p["radius"], 0.0)
        safe_radius = np.where(radius > 0, radius, 1.0)

        sigma_adsorber = g["ligand_charge"] * ELEMENTARY_CHARGE * AVOGADRO * g["ligand_density"]
        sigma_protein = p["charge"] * ELEMENTARY_CHARGE / (4.0 * np.pi * safe_radius**2)
        psi_adsorber = grahame_potential(sigma_adsorber, kappa, temperature, eps_r)
        psi_protein = grahame_potential(sigma_protein, kappa, temperature, eps_r)
        depth = well_depth(psi_adsorber, psi_protein, radius, temperature, eps_r)
        henry = p["Delta"] * henry_factor(depth)

        surface_density = (
            np.where(binding[:, np.newaxis], np.maximum(q, 0.0), 0.0) / g["specific_surface"]
        )
        crowding = available_surface_function(surface_density, radius)
        lateral_charge = (
            p["lateral_charge"] if "lateral_charge" in p.per_component else np.zeros_like(radius)
        )
        lateral = lateral_energy(surface_density, radius, lateral_charge, kappa, temperature, eps_r)
        return henry * crowding * np.exp(-np.minimum(lateral, MAX_WELL_DEPTH))

    def flux(self, c: np.ndarray, q: np.ndarray, p: BindingParameters) -> np.ndarray:
        """dq/dt = k_kin (K_H B e^{−u_lat/k_BT} c − q); zero rows for salt and non-binding ones."""
        binding = p.binding.copy()
        binding[0] = False
        rate = p["k_kin"] * (self.equilibrium_factor(c[0], q, p) * c - q)
        return np.where(binding[:, np.newaxis], rate, 0.0)
