"""Generalized ion exchange (GIEX): SMA with salt- and pH-dependent parameters (SPEC §5)."""

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import (
    K_KIN,
    CadetNativeIsotherm,
    per_component,
    rate_constants,
)
from chromsim.isotherms.registry import register_isotherm


@register_isotherm
class GeneralizedIonExchange(CadetNativeIsotherm):
    """Generalized ion exchange isotherm. Component 0 is salt, component 1 the pH pseudo-component.

    CADET's model (binding components i ≥ 2, salt c_{p,0}, pH c_{p,1}):

    .. math:: \\frac{dq_i}{dt} = k_{a,i}(c_{p,0}, c_{p,1})\\, c_{p,i}
              \\left(\\frac{\\bar q_0}{q_{ref}}\\right)^{\\nu_i(c_{p,1})}
              - k_{d,i}(c_{p,0}, c_{p,1})\\, q_i
              \\left(\\frac{c_{p,0}}{c_{ref}}\\right)^{\\nu_i(c_{p,1})}

    .. math:: \\bar q_0 = \\Lambda - \\sum_{j \\ge 2} (\\nu_j + \\sigma_j) q_j, \\qquad
              \\nu_i(c_{p,1}) = \\nu_i + \\nu_{i,lin}\\, c_{p,1} + \\nu_{i,quad}\\, c_{p,1}^2

    with k_a and k_d each multiplied by
    exp(lin·c_{p,1} + quad·c_{p,1}² + salt·c_{p,0}/c_{ref} + prot·c_{p,i}).

    ChromSim exposes the **equilibrium** dependence only: the pH and salt terms act on k_a,
    so that

    .. math:: \\frac{k_a}{k_d} = k_{eq}\\, \\exp\\left(k_{eq,lin}\\, pH + k_{eq,quad}\\, pH^2
              + k_{eq,salt}\\, \\frac{c_{p,0}}{c_{ref}}\\right)

    and k_d = k_kin is constant (cubic and protein terms are 0). As for SMA, c_ref = q_ref = Λ.
    Concentrations are in mol/m³ (= mM); pH is the pH number (unit ``pH``).

    References: Huuk, Hahn, Doninger, Griesbach, Hepbildikler & Hubbuch, Biotechnol. J. 12 (2017)
    1600336; Mollerup, Chem. Eng. Technol. 31 (2008) 864; CADET documentation (Generalized Ion
    Exchange).
    """

    name = "GIEX"
    cadet_class = "GeneralizedIonExchange"
    requires_salt = True
    requires_ionic_capacity = True
    n_modifiers = 2
    parameters = {
        "k_eq": Parameter(
            unit="dimensionless",
            description="equilibrium constant at pH 0, no salt term",
            log_scale=True,
        ),
        "nu": Parameter(
            unit="dimensionless",
            description="characteristic charge ν (constant part)",
            bounds=(None, None),
        ),
        "nu_lin": Parameter(
            unit="1/pH",
            description="linear pH dependence of ν",
            bounds=(None, None),
            required=False,
        ),
        "nu_quad": Parameter(
            unit="1/pH^2",
            description="quadratic pH dependence of ν",
            bounds=(None, None),
            required=False,
        ),
        "sigma": Parameter(unit="dimensionless", description="steric factor σ", bounds=(0, 500)),
        "k_eq_lin": Parameter(
            unit="1/pH",
            description="linear pH term of ln(k_a/k_d)",
            bounds=(None, None),
            required=False,
        ),
        "k_eq_quad": Parameter(
            unit="1/pH^2",
            description="quadratic pH term of ln(k_a/k_d)",
            bounds=(None, None),
            required=False,
        ),
        "k_eq_salt": Parameter(
            unit="dimensionless",
            description="salt term of ln(k_a/k_d), per c_s/Λ",
            bounds=(None, None),
            required=False,
        ),
        "k_kin": K_KIN,
    }
    example = {
        "components": [{"name": "salt"}, {"name": "pH"}, {"name": "A"}],
        "parameters": {"A": {"k_eq": 1e-2, "nu": 5.0, "nu_lin": 0.1, "sigma": 20.0, "k_kin": 10.0}},
        "load": {"salt": "50 mM", "pH": "5 pH", "A": "0.1 mM"},
        "elute": {"salt": "1000 mM", "pH": "5 pH"},
    }

    def bound_states(self, component_names, parameters):
        """Salt has a bound state (free sites), pH none, other components if they bind."""
        return [1, 0] + [1 if name in parameters else 0 for name in component_names[2:]]

    def modifier_components(self, component_names):
        """Salt and pH (components 0 and 1) modulate binding."""
        return component_names[:2]

    def initial_bound_state(self, component_names, parameters, ionic_capacity_mol_per_m3):
        """Equilibrated, protein-free column: all sites carry salt, q_0 = Λ."""
        q = [0.0] * sum(self.bound_states(component_names, parameters))
        q[0] = float(ionic_capacity_mol_per_m3)
        return q

    def cadet_parameters(
        self, component_names, parameters, ionic_capacity_mol_per_m3, global_parameters=None
    ):
        capacity = float(ionic_capacity_mol_per_m3)

        def value(name):
            return per_component(component_names, parameters, lambda v: v.get(name, 0.0))

        zeros = [0.0] * len(component_names)
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            "adsorption_rate_linear": value("k_eq_lin"),
            "adsorption_rate_quadratic": value("k_eq_quad"),
            "adsorption_rate_cubic": zeros,
            "adsorption_rate_salt": value("k_eq_salt"),
            "adsorption_rate_protein": zeros,
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
            "desorption_rate_linear": zeros,
            "desorption_rate_quadratic": zeros,
            "desorption_rate_cubic": zeros,
            "desorption_rate_salt": zeros,
            "desorption_rate_protein": zeros,
            "characteristic_charge": value("nu"),  # one pH piece (n_pieces = 1): flat list
            "characteristic_charge_linear": value("nu_lin"),
            "characteristic_charge_quadratic": value("nu_quad"),
            "characteristic_charge_cubic": zeros,
            "steric_factor": value("sigma"),
            "capacity": capacity,
            "reference_liquid_phase_conc": capacity,
            "reference_solid_phase_conc": capacity,
        }
