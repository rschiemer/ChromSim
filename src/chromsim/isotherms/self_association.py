"""Self-association steric mass action (SAS) isotherm: SMA plus protein binding to bound protein.

Mollerup (2008) added a term in which a protein adsorbs onto an already adsorbed protein of the
same kind (dimerisation on the surface, "self-association"). The adsorption rate then grows with
the liquid concentration, so the isotherm turns *convex* (anti-Langmuirian) at moderate
loading before the steric SMA term makes it saturate. Peaks then move to later retention with
increasing load before they turn back (anti-Langmuir-to-Langmuir, AL-L, behaviour), as reported
for bispecific antibodies and Fc-fusion proteins (see docs/theory/isotherms/sas.md).

CADET-Core implements the model as ``SELF_ASSOCIATION``. The flux below repeats its equation
(``SelfAssociationBinding.cpp``, v5.1), so the model also runs on the reference backend.
"""

import numpy as np

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import (
    K_KIN,
    CadetNativeIsotherm,
    per_component,
    rate_constants,
)
from chromsim.isotherms.registry import register_isotherm


@register_isotherm
class SelfAssociationSMA(CadetNativeIsotherm):
    """Self-association steric mass action (SAS) isotherm. Salt must be component 0.

    .. math:: \\frac{dq_i}{dt} = k_{a,i}\\,(1 + k_{sa,i} c_i)\\, c_i
              \\left(\\frac{\\bar q_0}{q_{ref}}\\right)^{\\nu_i}
              - k_{d,i} q_i \\left(\\frac{c_0}{c_{ref}}\\right)^{\\nu_i},
              \\qquad \\bar q_0 = \\Lambda - \\sum_{j \\ge 1} (\\nu_j + \\sigma_j) q_j

    At equilibrium (c_ref = q_ref = Λ, as for ChromSim's SMA):

    .. math:: q_i = k_{eq,i}\\,(1 + k_{sa,i} c_i)
              \\left(\\frac{\\bar q_0}{c_0}\\right)^{\\nu_i} c_i

    k_sa = k_a2/k_a1 (m³/mol) is the self-association constant: the factor (1 + k_sa c) raises the
    affinity with concentration and makes the isotherm convex. With k_sa = 0 this is exactly SMA.
    k_eq = k_a/k_d and k_kin = k_d as for SMA (omit k_kin for rapid equilibrium).

    References: Mollerup, Chem. Eng. Technol. 31 (2008) 864 ("A review of the thermodynamics of
    protein association to ligands, protein adsorption, and adsorption isotherms"); Westerberg,
    Broberg Hansen, Degerman, Budde Hansen & Nilsson, Chem. Eng. Technol. 35 (2012) 183;
    Seelinger, Wittkopp, von Hirschheydt & Frech, J. Chromatogr. A 1689 (2022) (anti-Langmuir
    elution of a bispecific antibody); CADET documentation (Self Association).
    """

    name = "SAS"
    cadet_class = "SelfAssociation"
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {
            "A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "self_association": 10.0, "k_kin": 10.0}
        },
        "load": {"salt": "50 mM", "A": "0.1 mM"},
        "elute": {"salt": "1000 mM"},
    }
    requires_salt = True
    requires_ionic_capacity = True
    n_modifiers = 1
    parameters = {
        "k_eq": Parameter(unit="dimensionless", description="equilibrium constant", log_scale=True),
        "nu": Parameter(unit="dimensionless", description="characteristic charge", bounds=(0, 50)),
        "sigma": Parameter(unit="dimensionless", description="steric factor", bounds=(0, 500)),
        "self_association": Parameter(
            unit="m^3/mol",
            description="self-association constant k_sa = k_a2/k_a (0 gives SMA)",
        ),
        "k_kin": K_KIN,
    }

    def modifier_components(self, component_names):
        """Salt (component 0) modulates binding."""
        return [component_names[0]]

    def bound_states(self, component_names, parameters):
        """Salt (component 0) always has a bound state; others only if they bind."""
        return [1] + [1 if name in parameters else 0 for name in component_names[1:]]

    def initial_bound_state(self, component_names, parameters, ionic_capacity_mol_per_m3):
        """Equilibrated, protein-free column: all sites carry salt, q_0 = Λ."""
        q = [0.0] * sum(self.bound_states(component_names, parameters))
        q[0] = float(ionic_capacity_mol_per_m3)
        return q

    def cadet_parameters(
        self, component_names, parameters, ionic_capacity_mol_per_m3, global_parameters=None
    ):
        capacity = float(ionic_capacity_mol_per_m3)
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            # k_a2 = k_sa · k_a (CADET-Core SelfAssociationBinding.cpp: (k_a1 + k_a2 c) c …)
            "adsorption_rate_dimerization": per_component(
                component_names,
                parameters,
                lambda v: v["self_association"] * rate_constants(v)[0],
            ),
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
            "characteristic_charge": per_component(component_names, parameters, lambda v: v["nu"]),
            "steric_factor": per_component(component_names, parameters, lambda v: v["sigma"]),
            "capacity": capacity,
            "reference_liquid_phase_conc": capacity,
            "reference_solid_phase_conc": capacity,
        }

    def flux(self, c, q, p):
        """SAS kinetics with c_ref = q_ref = Λ (reference backend).

        dq_i/dt = k_a (1 + k_sa c_i) c_i (q̄_0/Λ)^ν_i − k_d q_i (c_0/Λ)^ν_i for binding i, with
        q̄_0 = q_0 − Σ σ_j q_j; the bound salt follows from electroneutrality,
        dq_0/dt = −Σ ν_j dq_j/dt (as for SMA).
        """
        capacity = p.ionic_capacity
        nu, sigma = p["nu"], p["sigma"]
        k_a, k_d = p["k_eq"] * p["k_kin"], p["k_kin"]
        k_sa = p["self_association"]
        binding = p.binding.copy()
        binding[0] = False  # salt
        free_sites = q[0] - np.sum(sigma[binding] * q[binding], axis=0)
        salt = np.maximum(c[0], 0.0)
        rate = (
            k_a * (1.0 + k_sa * c) * c * (np.maximum(free_sites, 0.0) / capacity) ** nu
            - k_d * q * (salt / capacity) ** nu
        )
        dq = np.where(binding[:, np.newaxis], rate, 0.0)
        dq[0] = -np.sum(nu[binding] * dq[binding], axis=0)
        return dq
