"""Bi-Langmuir isotherm: two independent types of binding site (SPEC §5, affinity resins)."""

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.cadet_native import CadetNativeIsotherm, rate_constants
from chromsim.isotherms.registry import register_isotherm

N_SITES = 2


def _site(values: dict[str, float], site: int) -> dict[str, float]:
    """Parameters of one site in the single-site names (k_eq, q_max, k_kin)."""
    entry = {"k_eq": values[f"k_eq_{site}"]}
    if f"k_kin_{site}" in values:
        entry["k_kin"] = values[f"k_kin_{site}"]
    return entry


@register_isotherm
class BiLangmuir(CadetNativeIsotherm):
    """Multi-component Bi-Langmuir isotherm with two site types j = 1, 2.

    .. math:: \\frac{dq_{i,j}}{dt} = k_{a,i}^{(j)}\\, c_{p,i}\\, q_{max,i}^{(j)}
              \\left(1 - \\sum_k \\frac{q_{k,j}}{q_{max,k}^{(j)}}\\right) - k_{d,i}^{(j)}\\, q_{i,j}

    Components compete only for sites of the same type. Per binding component: k_eq_1, q_max_1,
    k_eq_2, q_max_2 and optionally k_kin_1, k_kin_2 (k_eq = k_a/k_d in m³/mol, q_max in mol/m³
    of solid).

    References: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
    Chromatography*, 2nd ed., Academic Press 2006; CADET documentation (Multi Component
    Bi-Langmuir).
    """

    name = "BiLangmuir"
    cadet_class = "BiLangmuir"
    parameters = {
        "k_eq_1": Parameter(
            unit="m^3/mol", description="equilibrium constant, site 1", log_scale=True
        ),
        "q_max_1": Parameter(unit="mol/m^3", description="capacity of site 1"),
        "k_eq_2": Parameter(
            unit="m^3/mol", description="equilibrium constant, site 2", log_scale=True
        ),
        "q_max_2": Parameter(unit="mol/m^3", description="capacity of site 2"),
        "k_kin_1": Parameter(
            unit="1/s", description="desorption rate, site 1", log_scale=True, required=False
        ),
        "k_kin_2": Parameter(
            unit="1/s", description="desorption rate, site 2", log_scale=True, required=False
        ),
    }
    example = {
        "components": [{"name": "A"}],
        "parameters": {
            "A": {
                "k_eq_1": 5.0,
                "q_max_1": 5.0,
                "k_eq_2": 0.5,
                "q_max_2": 10.0,
                "k_kin_1": 1.0,
                "k_kin_2": 1.0,
            }
        },
        "load": {"A": "1 mM"},
        "elute": {},
    }

    def bound_states(self, component_names, parameters):
        """One bound state per component and site, state-major: [site 1 …, site 2 …]."""
        per_site = [1 if name in parameters else 0 for name in component_names]
        return per_site * N_SITES

    def cadet_parameters(
        self, component_names, parameters, ionic_capacity_mol_per_m3, global_parameters=None
    ):
        # CADET-Process stores one entry per *bound* state, state-major (all components of
        # site 1, then of site 2).
        adsorption, desorption, capacity = [], [], []
        for site in range(1, N_SITES + 1):
            for name in component_names:
                if name not in parameters:
                    continue
                k_a, k_d = rate_constants(_site(parameters[name], site))
                adsorption.append(k_a)
                desorption.append(k_d)
                capacity.append(parameters[name][f"q_max_{site}"])
        return {"adsorption_rate": adsorption, "desorption_rate": desorption, "capacity": capacity}
