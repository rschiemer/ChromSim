"""Worked example of a Python isotherm: a mixed-mode (ion exchange + hydrophobic) Langmuir model.

The equations are those of CADET's mobile phase modulator Langmuir isotherm, so the Python route
can be checked against CADET (tests/backends/test_reference_vs_cadet.py). Use it as a template:
copy this file, change ``name``, ``parameters``, ``example`` and ``flux``, and add a docs page.
"""

import numpy as np

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.python_isotherm import BindingParameters, PythonIsotherm
from chromsim.isotherms.registry import register_isotherm


@register_isotherm
class MixedModeLangmuir(PythonIsotherm):
    """Mixed-mode Langmuir isotherm (Python). Salt (component 0) modulates binding.

    .. math:: \\frac{dq_i}{dt} = k_{a,i}\\, e^{\\gamma_i c_{p,0}}\\, c_{p,i}\\, q_{max,i}
              \\left(1 - \\sum_{j \\ge 1} \\frac{q_j}{q_{max,j}}\\right)
              - k_{d,i}\\, c_{p,0}^{\\beta_i}\\, q_i

    Hydrophobic interaction strengthens binding with salt (e^{γ c_s}); electrostatic interaction
    weakens it (c_s^β). k_a = k_eq k_kin, k_d = k_kin. Concentrations in mol/m³ (= mM).

    References: Melander, El Rassi & Horváth, J. Chromatogr. 469 (1989) 3 (retention model for
    combined hydrophobic and electrostatic interactions); Karlsson, Jakobsson, Axelsson &
    Nilsson, J. Chromatogr. A 1055 (2004) 29 (kinetic Langmuir form).
    """

    name = "MixedModeLangmuir"
    requires_salt = True
    n_modifiers = 1
    parameters = {
        "k_eq": Parameter(
            unit="m^3/mol", description="equilibrium constant k_a/k_d", log_scale=True
        ),
        "q_max": Parameter(unit="mol/m^3", description="saturation capacity (per solid volume)"),
        "beta": Parameter(
            unit="dimensionless",
            description="electrostatic (ion-exchange) exponent β",
            bounds=(0, 50),
        ),
        "gamma": Parameter(
            unit="m^3/mol", description="hydrophobic coefficient γ", bounds=(None, None)
        ),
        "k_kin": Parameter(unit="1/s", description="desorption rate constant k_d", log_scale=True),
    }
    example = {
        "components": [{"name": "salt"}, {"name": "A"}],
        "parameters": {"A": {"k_eq": 1e3, "q_max": 10.0, "beta": 2.0, "gamma": 1e-3, "k_kin": 1.0}},
        "load": {"salt": "50 mM", "A": "0.1 mM"},
        "elute": {"salt": "1000 mM"},
    }

    def bound_states(self, component_names, parameters):
        """Salt has no bound state."""
        return [0] + [1 if name in parameters else 0 for name in component_names[1:]]

    def flux(self, c: np.ndarray, q: np.ndarray, p: BindingParameters) -> np.ndarray:
        salt = np.maximum(c[0], 0.0)  # c_{p,0}
        binding = p.binding.copy()
        binding[0] = False
        k_a, k_d = p["k_eq"] * p["k_kin"], p["k_kin"]
        q_max, beta, gamma = p["q_max"], p["beta"], p["gamma"]
        free = 1.0 - np.sum(q[binding] / q_max[binding], axis=0)
        rate = k_a * np.exp(gamma * salt) * c * q_max * free - k_d * salt**beta * q
        return np.where(binding[:, np.newaxis], rate, 0.0)
