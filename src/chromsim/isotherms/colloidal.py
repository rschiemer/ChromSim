"""Colloidal isotherm of Xu & Lenhoff (2009), as implemented in CADET-Core (SPEC §5).

Proteins are treated as charged colloidal spheres. The adsorbed layer is a hexagonal lattice whose
spacing R shrinks as the surface fills; neighbouring proteins repel each other through a screened
(Yukawa) potential. This gives a non-stoichiometric alternative to SMA for ion exchange and, with
a different salt/pH dependence of K_e, for Protein A (Bhoyar et al. 2024).

The flux below reproduces CADET-Core's ``ColloidalBinding.cpp`` (v5.1) line by line, so the
model also runs on the reference backend; ``tests/backends/test_reference_vs_cadet.py`` checks
both against each other.
"""

import numpy as np

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import CadetNativeIsotherm, per_component
from chromsim.isotherms.python_isotherm import BindingParameters
from chromsim.isotherms.registry import register_isotherm

# 2 / (N_A √3), m²·mol: R² = R_FACTOR · φ / q_tot (CADET-Core ColloidalBinding.cpp, `_rFactor`).
AVOGADRO = 6.02214076e23  # 1/mol
R_FACTOR = 2.0 / (AVOGADRO * np.sqrt(3.0))


@register_isotherm
class Colloidal(CadetNativeIsotherm):
    """Colloidal isotherm (Xu & Lenhoff 2009). Salt is component 0 and has no bound state.

    .. math::

        \\frac{dq_i}{dt} = k_{kin,i}\\left(c_{p,i} - c_{p,i}^\\star\\right), \\qquad
        c_{p,i}^\\star = q_i \\exp\\left[\\frac{n(3 + \\kappa R)}{4 q_{tot}}
        \\sum_j q_j \\sqrt{b_{pp,i} b_{pp,j}}\\, \\frac{r_i + r_j}{2R}\\,
        e^{-\\kappa (R - r_i - r_j)} - \\ln K_{e,i}\\right]

        R = \\sqrt{\\frac{2\\phi}{\\sqrt{3}\\, N_A\\, q_{tot}}}, \\qquad
        \\kappa = \\frac{1}{\\kappa_f c_{p,0}^{\\kappa_{ef}} + \\kappa_c}

        \\ln K_{e,i} = k_{a,i}\\, c_{p,0}^{-k_{b,i}} + k_{c,i}\\, e^{k_{d,i} c_{p,0}}, \\qquad
        b_{pp,i} = b_{a,i}\\, c_{p,0}^{b_{b,i}} + b_{c,i}\\, e^{b_{d,i} c_{p,0}}

    q_tot = Σ q_j over binding components. Below ``linear_threshold`` (q_tot ≤ c_ε) the lattice
    term is dropped: c*_i = q_i / K_e,i. At low loading the isotherm is linear with Henry constant
    K_e; the lattice term makes it saturate without a fixed capacity.

    **Kinetic constant.** Here k_kin multiplies the liquid-phase driving force c − c*, as in
    CADET; it is *not* the desorption rate of the other ChromSim isotherms (ADR 0011). Omit k_kin
    for rapid equilibrium. Concentrations are in mol/m³ (= mM); κ_f is in m·(m³/mol)^κ_ef.

    References: Xu & Lenhoff, J. Chromatogr. A 1216 (2009) 6177 ("Binary adsorption of globular
    proteins on ion-exchange media"); Oberholzer & Lenhoff, Langmuir 15 (1999) 3905; Kumar et al.,
    J. Chromatogr. A 1660 (2021) 462669; CADET documentation (Multi Component Colloidal).
    """

    name = "Colloidal"
    cadet_class = "MultiComponentColloidal"
    requires_salt = True
    n_modifiers = 1
    parameters = {
        # --- per binding component -------------------------------------------------------
        "radius": Parameter(unit="m", description="protein radius r_i", bounds=(0.0, None)),
        "ln_ke_salt_factor": Parameter(
            unit="dimensionless",
            description="k_a in ln K_e = k_a c_s^(−k_b) + k_c exp(k_d c_s)",
            bounds=(None, None),
        ),
        "ln_ke_salt_exponent": Parameter(
            unit="dimensionless",
            description="k_b, exponent of salt in the power term of ln K_e",
            bounds=(None, None),
        ),
        "ln_ke_exp_factor": Parameter(
            required=False,
            unit="dimensionless",
            description="k_c, factor of the exponential term of ln K_e (default 0)",
            bounds=(None, None),
        ),
        "ln_ke_exp_salt_coefficient": Parameter(
            required=False,
            unit="m^3/mol",
            description="k_d, salt coefficient in the exponential term of ln K_e (default 0)",
            bounds=(None, None),
        ),
        "b_pp_salt_factor": Parameter(
            unit="dimensionless",
            description="b_a in b_pp = b_a c_s^(b_b) + b_c exp(b_d c_s) (protein–protein)",
            bounds=(None, None),
        ),
        "b_pp_salt_exponent": Parameter(
            unit="dimensionless",
            description="b_b, exponent of salt in the power term of b_pp",
            bounds=(None, None),
        ),
        "b_pp_exp_factor": Parameter(
            required=False,
            unit="dimensionless",
            description="b_c, factor of the exponential term of b_pp (default 0)",
            bounds=(None, None),
        ),
        "b_pp_exp_salt_coefficient": Parameter(
            required=False,
            unit="m^3/mol",
            description="b_d, salt coefficient in the exponential term of b_pp (default 0)",
            bounds=(None, None),
        ),
        "k_kin": Parameter(
            unit="1/s",
            description="rate constant of the driving force c − c* (omit for rapid equilibrium)",
            log_scale=True,
            required=False,
        ),
        # --- whole isotherm ----------------------------------------------------------------
        "phase_ratio": Parameter(
            unit="m^2/m^3",
            description="adsorber surface area per solid volume φ",
            log_scale=True,
            per_component=False,
        ),
        "kappa_factor": Parameter(
            unit="m",
            description="κ_f in 1/κ = κ_f c_s^κ_ef + κ_c (Debye length, m; c_s in mM)",
            per_component=False,
        ),
        "kappa_exponent": Parameter(
            unit="dimensionless",
            description="κ_ef; CADET-Process 0.12 allows only κ_ef ≥ 0",
            per_component=False,
        ),
        "kappa_constant": Parameter(
            unit="m", description="κ_c, constant part of the Debye length", per_component=False
        ),
        "coordination_number": Parameter(
            unit="dimensionless",
            description="nearest neighbours n in the adsorbed lattice (6: hexagonal)",
            bounds=(1, 12),
            per_component=False,
            required=False,
        ),
        "linear_threshold": Parameter(
            unit="mol/m^3",
            description="q_tot below which the linear form is used c_ε (default 1e-8)",
            bounds=(0.0, None),
            per_component=False,
            required=False,
        ),
    }
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {
            "A": {
                "radius": 4.5e-9,
                "ln_ke_salt_factor": 60.0,
                "ln_ke_salt_exponent": 0.5,
                "b_pp_salt_factor": 10.0,
                "b_pp_salt_exponent": 0.0,
                "k_kin": 1.0,
            }
        },
        "global_parameters": {
            "phase_ratio": 4e7,
            "kappa_factor": 0.0,
            "kappa_exponent": 0.0,
            "kappa_constant": 1e-9,
        },
        "load": {"salt": "50 mM", "A": "0.05 mM"},
        "elute": {"salt": "1000 mM"},
    }

    def bound_states(self, component_names, parameters):
        """Salt has no bound state (CADET requires NBOUND = 0 for it)."""
        return [0] + [1 if name in parameters else 0 for name in component_names[1:]]

    def modifier_components(self, component_names):
        """Salt (component 0) modulates binding."""
        return [component_names[0]]

    def cadet_parameters(
        self, component_names, parameters, ionic_capacity_mol_per_m3, global_parameters=None
    ):
        g = global_parameters or {}

        def values(name):
            return per_component(component_names, parameters, lambda v: v.get(name, 0.0))

        return {
            "phase_ratio": float(g["phase_ratio"]),
            "kappa_factor": float(g["kappa_factor"]),
            "kappa_exponential": float(g["kappa_exponent"]),
            "kappa_constant": float(g["kappa_constant"]),
            "coordination_number": int(g.get("coordination_number", 6)),
            "linear_threshold": float(g.get("linear_threshold", 1e-8)),
            "use_ph": False,
            # Names verified in CADET-Core ColloidalBinding.cpp (v5.1):
            # ln K_e = c^(−POWEXP)·POWFACT + exp(c·EXPARGMULT)·EXPFACT. CADET's docs describe
            # POWEXP/POWFACT the other way round; tests/backends pins the behaviour.
            "logkeq_power_factor": values("ln_ke_salt_factor"),
            "logkeq_power_exponent": values("ln_ke_salt_exponent"),
            "logkeq_exponent_factor": values("ln_ke_exp_factor"),
            "logkeq_exponent_multiplier": values("ln_ke_exp_salt_coefficient"),
            "logkeq_ph_exponent": [0.0] * len(component_names),
            "bpp_power_factor": values("b_pp_salt_factor"),
            "bpp_power_exponent": values("b_pp_salt_exponent"),
            "bpp_exponent_factor": values("b_pp_exp_factor"),
            "bpp_exponent_multiplier": values("b_pp_exp_salt_coefficient"),
            "bpp_ph_exponent": [0.0] * len(component_names),
            "protein_radius": values("radius"),
            # Rapid equilibrium: CADET still needs a value; it is not used.
            "kinetic_rate_constant": per_component(
                component_names, parameters, lambda v: v.get("k_kin", 1.0)
            ),
        }

    # ------------------------------------------------------------------------------------
    # Reference backend (same equations as CADET-Core)
    # ------------------------------------------------------------------------------------

    def flux(self, c: np.ndarray, q: np.ndarray, p: BindingParameters) -> np.ndarray:
        """dq/dt = k_kin (c − c*), shape (n_components, n_points); see the class docstring."""
        g = p.global_
        binding = p.binding.copy()
        binding[0] = False  # salt
        c_s = np.maximum(c[0], 1e-12)  # c_{p,0}, mol/m³ (CADET evaluates c_s^(−k_b) as well)

        def per(name):
            return p[name] if name in p.per_component else np.zeros((len(binding), 1))

        ln_ke = per("ln_ke_salt_factor") * c_s ** (-per("ln_ke_salt_exponent")) + per(
            "ln_ke_exp_factor"
        ) * np.exp(per("ln_ke_exp_salt_coefficient") * c_s)
        b_pp = per("b_pp_salt_factor") * c_s ** per("b_pp_salt_exponent") + per(
            "b_pp_exp_factor"
        ) * np.exp(per("b_pp_exp_salt_coefficient") * c_s)
        r = per("radius")

        q_b = np.where(binding[:, np.newaxis], q, 0.0)
        q_tot = q_b.sum(axis=0)
        threshold = g.get("linear_threshold", 1e-8)
        lattice = q_tot > threshold
        q_safe = np.where(lattice, q_tot, 1.0)

        kappa = 1.0 / (g["kappa_factor"] * c_s ** g["kappa_exponent"] + g["kappa_constant"])
        R = np.sqrt(R_FACTOR * g["phase_ratio"] / q_safe)
        s_factor = g.get("coordination_number", 6) * 0.125 / q_safe * (3.0 / R + kappa)

        # S_i = Σ_j q_j √(b_pp,i b_pp,j) (r_i + r_j) exp(−κ (R − r_i − r_j)); CADET skips pairs
        # with b_pp,i b_pp,j ≤ 1e-16.
        pair_b = b_pp[:, np.newaxis, :] * b_pp[np.newaxis, :, :]  # (i, j, points)
        r_sum = (r[:, np.newaxis, :] + r[np.newaxis, :, :]) * np.ones_like(pair_b)
        pair = np.where(
            pair_b > 1e-16,
            np.sqrt(np.maximum(pair_b, 0.0)) * r_sum * np.exp(-kappa * (R - r_sum)),
            0.0,
        )
        S = np.einsum("ijp,jp->ip", pair, q_b)

        exponent = np.where(lattice, S * s_factor, 0.0) - ln_ke
        c_star = q * np.exp(exponent)
        rate = per("k_kin") * (c - c_star)
        return np.where(binding[:, np.newaxis], rate, 0.0)
