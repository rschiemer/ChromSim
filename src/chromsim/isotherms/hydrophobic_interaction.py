"""Hydrophobic interaction chromatography (HIC) isotherms of CADET (SPEC §5).

Both models share the parameters: per binding component k_eq, q_max, ν (number of ligands per
interaction) and optionally k_kin; isotherm-wide β₀ and β₁ (``isotherm.global_parameters``) for
the salt dependence of the number of ordered water molecules, β = β₀ exp(β₁ c_{p,0}).
Salt is component 0 and has no bound state. Concentrations are in mol/m³ (= mM).
"""

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import (
    K_KIN,
    CadetNativeIsotherm,
    per_component,
    rate_constants,
)
from chromsim.isotherms.registry import register_isotherm


class _HICBase(CadetNativeIsotherm):
    requires_salt = True
    n_modifiers = 1
    parameters = {
        "k_eq": Parameter(
            unit="dimensionless", description="equilibrium constant k_a/k_d", log_scale=True
        ),
        "q_max": Parameter(unit="mol/m^3", description="saturation capacity (per solid volume)"),
        "nu": Parameter(
            unit="dimensionless",
            description="ligands per protein–ligand interaction ν",
            bounds=(0, 100),
        ),
        "k_kin": K_KIN,
        "beta_0": Parameter(
            unit="dimensionless",
            description="ordered water molecules at infinite dilution β₀",
            per_component=False,
        ),
        "beta_1": Parameter(
            unit="m^3/mol", description="salt dependence of ordered water β₁", per_component=False
        ),
    }
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {"A": {"k_eq": 1.0, "q_max": 10.0, "nu": 2.0, "k_kin": 1.0}},
        "global_parameters": {"beta_0": 0.1, "beta_1": 1e-3},
        # HIC binds at high salt and elutes at low salt.
        "load": {"salt": "1000 mM", "A": "0.1 mM"},
        "elute": {"salt": "0 mM"},
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
        global_parameters = global_parameters or {}
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
            "capacity": per_component(component_names, parameters, lambda v: v["q_max"]),
            "hic_characteristic": per_component(component_names, parameters, lambda v: v["nu"]),
            "beta_0": float(global_parameters["beta_0"]),
            "beta_1": float(global_parameters["beta_1"]),
        }


@register_isotherm
class HICConstantWaterActivity(_HICBase):
    """HIC isotherm with constant water activity.

    .. math:: \\frac{dq_i}{dt} = k_{a,i}\\, c_{p,i}
              \\left(1 - \\sum_j \\frac{q_j}{q_{max,j}}\\right)^{\\nu_i}
              - k_{d,i}\\, q_i\\, 0.1^{\\nu_i \\beta},
              \\qquad \\beta = \\beta_0 e^{\\beta_1 c_{p,0}}

    Reference: Jäpel & Buyel, J. Chromatogr. A 1679 (2022) 463408; CADET documentation (HIC
    Constant Water Activity).
    """

    name = "HICConstantWaterActivity"
    cadet_class = "HICConstantWaterActivity"


@register_isotherm
class HICWaterOnHydrophobicSurfaces(_HICBase):
    """HIC isotherm "water on hydrophobic surfaces" (slightly modified multi-component form).

    .. math:: \\frac{dq_i}{dt} = k_{a,i}\\, c_{p,i}
              \\left(1 - \\sum_j \\frac{q_j}{q_{max,j}}\\right)^{\\nu_i}
              - k_{d,i}\\, q_i \\left(\\sum_j q_j\\right)^{\\nu_i \\beta},
              \\qquad \\beta = \\beta_0 e^{\\beta_1 c_{p,0}}

    Reference: Wang, Hahn & Hubbuch, J. Chromatogr. A 1465 (2016) 71; CADET documentation (HIC
    Water on Hydrophobic Surfaces).
    """

    name = "HICWaterOnHydrophobicSurfaces"
    cadet_class = "HICWaterOnHydrophobicSurfaces"
