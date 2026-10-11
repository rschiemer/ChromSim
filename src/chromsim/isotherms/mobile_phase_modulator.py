"""Mobile phase modulator (MPM) Langmuir isotherm: salt-dependent ion exchange and hydrophobic
interaction in one Langmuir model. Used for HIC and as a mixed-mode approximation (SPEC §5)."""

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import (
    K_KIN,
    CadetNativeIsotherm,
    per_component,
    rate_constants,
)
from chromsim.isotherms.registry import register_isotherm


@register_isotherm
class MobilePhaseModulator(CadetNativeIsotherm):
    """Mobile phase modulator Langmuir isotherm. Salt (component 0) is the modulator.

    .. math:: \\frac{dq_i}{dt} = k_{a,i}\\, e^{\\gamma_i c_{p,0}}\\, c_{p,i}\\, q_{max,i}
              \\left(1 - \\sum_{j \\ge 1} \\frac{q_j}{q_{max,j}}\\right)
              - k_{d,i}\\, c_{p,0}^{\\beta_i}\\, q_i

    The exponential term (γ, hydrophobicity) strengthens binding with salt as in HIC; the power
    law (β, ion-exchange characteristic) weakens it as in ion exchange. Salt is inert (no bound
    state). k_eq = k_a/k_d, k_kin = k_d (see :mod:`chromsim.isotherms.cadet_native`).
    Concentrations are in mol/m³ (= mM), so γ is in m³/mol and k_eq in m³/mol (with c_{p,0}^β in
    (mol/m³)^β).

    References: Melander, El Rassi & Horváth, J. Chromatogr. 469 (1989) 3 ("Interplay of
    hydrophobic and electrostatic interactions in biopolymer chromatography"); Karlsson,
    Jakobsson, Axelsson & Nilsson, J. Chromatogr. A 1055 (2004) 29; CADET documentation
    (Mobile Phase Modulator Langmuir).
    """

    name = "MPM"
    cadet_class = "MobilePhaseModulator"
    requires_salt = True
    n_modifiers = 1
    parameters = {
        "k_eq": Parameter(
            unit="m^3/mol", description="equilibrium constant k_a/k_d", log_scale=True
        ),
        "q_max": Parameter(unit="mol/m^3", description="saturation capacity (per solid volume)"),
        "beta": Parameter(
            unit="dimensionless", description="ion-exchange characteristic β", bounds=(0, 50)
        ),
        "gamma": Parameter(unit="m^3/mol", description="hydrophobicity γ", bounds=(None, None)),
        "k_kin": K_KIN,
    }
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {"A": {"k_eq": 1e3, "q_max": 10.0, "beta": 2.0, "gamma": 1e-3, "k_kin": 1.0}},
        "load": {"salt": "50 mM", "A": "0.1 mM"},
        "elute": {"salt": "1000 mM"},
    }

    def bound_states(self, component_names, parameters):
        """Salt has no bound state; other components bind if they have parameters."""
        return [0] + [1 if name in parameters else 0 for name in component_names[1:]]

    def modifier_components(self, component_names):
        """Salt (component 0) modulates binding."""
        return [component_names[0]]

    def cadet_parameters(
        self, component_names, parameters, ionic_capacity_mol_per_m3, global_parameters=None
    ):
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
            "capacity": per_component(component_names, parameters, lambda v: v["q_max"]),
            "ion_exchange_characteristic": per_component(
                component_names, parameters, lambda v: v["beta"]
            ),
            "hydrophobicity": per_component(component_names, parameters, lambda v: v["gamma"]),
        }
