"""Isotherms that CADET-Core implements natively, described for ChromSim.

Each class maps ChromSim's user-facing parameters (``k_eq``, ``k_kin``, ``q_max`` …) to the
attributes of a CADET-Process binding model. CADET-Process names were read from the installed
version (0.12) and are pinned by ``tests/backends/test_cadet_names.py`` (CLAUDE.md rule 4).

Kinetic convention (all isotherms here)
---------------------------------------
``k_eq`` is the equilibrium constant k_a/k_d and ``k_kin`` the desorption rate constant k_d:

.. math:: k_a = k_{eq} \\, k_{kin}, \\qquad k_d = k_{kin}

Without ``k_kin`` the isotherm is solved in rapid equilibrium; CADET then uses only the ratio
k_a/k_d, so ChromSim passes k_a = k_eq and k_d = 1.

Bound-phase concentrations q are per volume of solid phase (mol/m³), as in CADET.
"""

from abc import ABC, abstractmethod
from typing import ClassVar

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.registry import register_isotherm

K_KIN = Parameter(
    unit="1/s",
    description="desorption rate constant k_d (omit for rapid equilibrium)",
    bounds=(0.0, None),
    log_scale=True,
    required=False,
)


class CadetNativeIsotherm(ABC):
    """Base class for isotherms solved by CADET-Core.

    Subclasses set ``name``, ``cadet_class`` (the CADET-Process binding class name) and
    ``parameters``, and implement :meth:`cadet_parameters`.
    """

    name: ClassVar[str]
    cadet_class: ClassVar[str]
    parameters: ClassVar[dict[str, Parameter]]
    requires_salt: ClassVar[bool] = False
    requires_ionic_capacity: ClassVar[bool] = False

    def is_kinetic(self, parameters: dict[str, dict[str, float]]) -> bool:
        """``True`` if binding is kinetic (``k_kin`` given), ``False`` for rapid equilibrium."""
        return any("k_kin" in values for values in parameters.values())

    def bound_states(
        self, component_names: list[str], parameters: dict[str, dict[str, float]]
    ) -> list[int]:
        """Number of bound states per component: 1 for components with parameters, else 0."""
        return [1 if name in parameters else 0 for name in component_names]

    def modifier_components(self, component_names: list[str]) -> list[str]:
        """Mobile-phase modifiers (e.g. salt) the isotherm needs at t = 0. Default: none."""
        return []

    def initial_bound_state(
        self,
        component_names: list[str],
        parameters: dict[str, dict[str, float]],
        ionic_capacity_mol_per_m3: float | None,
    ) -> list[float]:
        """Bound-phase concentrations of an equilibrated, protein-free column, mol/m³.

        One value per bound state (in component order). Default: nothing bound.
        """
        return [0.0] * sum(self.bound_states(component_names, parameters))

    @abstractmethod
    def cadet_parameters(
        self,
        component_names: list[str],
        parameters: dict[str, dict[str, float]],
        ionic_capacity_mol_per_m3: float | None,
    ) -> dict:
        """Values for the CADET-Process binding model attributes.

        Parameters
        ----------
        component_names : list of str
            Project components in order.
        parameters : dict
            Component name → {parameter → value}; components without an entry don't bind.
        ionic_capacity_mol_per_m3 : float or None
            Resin ionic capacity Λ, mol/m³.

        Returns
        -------
        dict
            CADET-Process attribute name → value (lists have one entry per component).
        """


def rate_constants(values: dict[str, float]) -> tuple[float, float]:
    """(k_a, k_d) from k_eq and the optional k_kin (see the module docstring)."""
    if "k_kin" in values:
        return values["k_eq"] * values["k_kin"], values["k_kin"]
    return values["k_eq"], 1.0


def per_component(
    component_names: list[str], parameters: dict[str, dict[str, float]], value_of
) -> list[float]:
    """``value_of(values)`` for every binding component, 0.0 for non-binding components."""
    return [
        float(value_of(parameters[name])) if name in parameters else 0.0 for name in component_names
    ]


@register_isotherm
class Linear(CadetNativeIsotherm):
    """Linear isotherm.

    .. math:: \\frac{dq_i}{dt} = k_{a,i} c_i - k_{d,i} q_i,
              \\qquad q_i^* = k_{eq,i} \\, c_i \\text{ at equilibrium}

    Reference: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
    Chromatography*, 2nd ed., Academic Press 2006, ch. 3.
    """

    name = "Linear"
    cadet_class = "Linear"
    parameters = {
        "k_eq": Parameter(
            unit="dimensionless",
            description="Henry constant q*/c (q per solid volume)",
            log_scale=True,
        ),
        "k_kin": K_KIN,
    }

    def cadet_parameters(self, component_names, parameters, ionic_capacity_mol_per_m3):
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
        }


@register_isotherm
class Langmuir(CadetNativeIsotherm):
    """Multi-component Langmuir isotherm.

    .. math:: \\frac{dq_i}{dt} = k_{a,i} c_i q_{max,i}
              \\left(1 - \\sum_j \\frac{q_j}{q_{max,j}}\\right) - k_{d,i} q_i

    At equilibrium, single component: :math:`q^* = q_{max} k_{eq} c / (1 + k_{eq} c)`.

    Reference: Langmuir, J. Am. Chem. Soc. 40 (1918) 1361; multi-component form as implemented in
    CADET (Leweke & von Lieres, Comput. Chem. Eng. 113 (2018) 274).
    """

    name = "Langmuir"
    cadet_class = "Langmuir"
    parameters = {
        "k_eq": Parameter(
            unit="m^3/mol", description="equilibrium constant k_a/k_d (= 1/mM)", log_scale=True
        ),
        "q_max": Parameter(unit="mol/m^3", description="saturation capacity (per solid volume)"),
        "k_kin": K_KIN,
    }

    def cadet_parameters(self, component_names, parameters, ionic_capacity_mol_per_m3):
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
            ),
            "desorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[1]
            ),
            "capacity": per_component(component_names, parameters, lambda v: v["q_max"]),
        }


@register_isotherm
class StericMassAction(CadetNativeIsotherm):
    """Steric mass action (SMA) isotherm. Salt must be component 0.

    .. math:: \\frac{dq_i}{dt} = k_{a,i} c_i \\left(\\frac{\\bar q_0}{q_{ref}}\\right)^{\\nu_i}
              - k_{d,i} q_i \\left(\\frac{c_0}{c_{ref}}\\right)^{\\nu_i},
              \\qquad \\bar q_0 = \\Lambda - \\sum_{j \\ge 1} (\\nu_j + \\sigma_j) q_j

    with salt concentration c_0, ionic capacity Λ, characteristic charge ν and steric factor σ.
    ChromSim sets both reference concentrations to Λ (c_ref = q_ref = Λ). The equilibrium constant
    then keeps its classical meaning,

    .. math:: k_{eq,i} = \\frac{k_{a,i}}{k_{d,i}}
              = \\frac{q_i \\, c_0^{\\nu_i}}{c_i \\, \\bar q_0^{\\nu_i}},

    and k_kin = k_d is a rate in 1/s instead of a number scaled by (mM)^ν. The bound salt is
    q_0 = Λ − Σ ν_j q_j (electroneutrality).

    Reference: Brooks & Cramer, AIChE J. 38 (1992) 1969; reference-concentration form as in
    CADET (Leweke & von Lieres, Comput. Chem. Eng. 113 (2018) 274).
    """

    name = "SMA"
    cadet_class = "StericMassAction"
    requires_salt = True
    requires_ionic_capacity = True
    parameters = {
        "k_eq": Parameter(unit="dimensionless", description="equilibrium constant", log_scale=True),
        "nu": Parameter(unit="dimensionless", description="characteristic charge", bounds=(0, 50)),
        "sigma": Parameter(unit="dimensionless", description="steric factor", bounds=(0, 500)),
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

    def cadet_parameters(self, component_names, parameters, ionic_capacity_mol_per_m3):
        capacity = float(ionic_capacity_mol_per_m3)
        return {
            "adsorption_rate": per_component(
                component_names, parameters, lambda v: rate_constants(v)[0]
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
