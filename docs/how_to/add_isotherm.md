# How to add an isotherm

An isotherm is **one new file + a docs page**. Registration and the shared contract tests do the
rest (CLAUDE.md rule 6). Don't edit core modules.

There are two kinds:

| Kind | When | Runs on |
|---|---|---|
| **CADET-native** | CADET already implements the model (check `CADETProcess.processModel`) | CADET backend |
| **Python** | your own or a literature model CADET lacks | reference backend (LRM, LRMP) |

## A Python isotherm

Copy `src/chromsim/isotherms/custom/mixed_mode_langmuir.py` to
`src/chromsim/isotherms/custom/my_isotherm.py` and adapt it:

```python
"""One-line description of the model."""

import numpy as np

from chromsim.isotherms.base import Parameter
from chromsim.isotherms.python_isotherm import BindingParameters, PythonIsotherm
from chromsim.isotherms.registry import register_isotherm


@register_isotherm
class MyIsotherm(PythonIsotherm):
    """Short name of the model.

    .. math:: \\frac{dq_i}{dt} = \\dots

    Reference: Author, Journal Volume (Year) Page.
    """

    name = "MyIsotherm"                      # used in project files: isotherm: {model: MyIsotherm}
    parameters = {
        "k_eq": Parameter(unit="m^3/mol", description="equilibrium constant", log_scale=True),
        "q_max": Parameter(unit="mol/m^3", description="capacity"),
        "k_kin": Parameter(unit="1/s", description="desorption rate", log_scale=True),
    }
    example = {                              # used by the contract tests
        "components": [{"name": "A"}],
        "parameters": {"A": {"k_eq": 1.0, "q_max": 10.0, "k_kin": 1.0}},
        "load": {"A": "1 mM"},
        "elute": {},
    }

    def flux(self, c: np.ndarray, q: np.ndarray, p: BindingParameters) -> np.ndarray:
        # c, q: shape (n_components, n_points), mol/m³; p["k_eq"]: shape (n_components, 1)
        k_a, k_d = p["k_eq"] * p["k_kin"], p["k_kin"]
        rate = k_a * c * (p["q_max"] - q) - k_d * q
        return np.where(p.binding[:, np.newaxis], rate, 0.0)   # zero for non-binding components
```

Rules:

- **Units:** use concentrations in mol/m³ (= mM) and time in s, so values from project files pass
  through unchanged. State the unit of every parameter.
- **Shape:** `flux` returns dq/dt with the same shape as `q`, and zero rows for components
  without parameters.
- **Salt or pH first?** Set `requires_salt = True` and `n_modifiers = 1` (salt) or `2` (salt and
  pH). Return 0 bound states for them in `bound_states` unless the model tracks bound salt.
- **Global parameters** (one value for the whole isotherm): `Parameter(..., per_component=False)`.
  Users give them under `isotherm.global_parameters`.

## A CADET-native isotherm

1. Look up the class and its attribute names in the installed CADET-Process. Never write them
   from memory (CLAUDE.md rule 4):

   ```bash
   pixi run python -c "import CADETProcess.processModel as pm; print(list(pm.MyModel(pm.ComponentSystem(2), name='b').parameters))"
   ```

2. Subclass `CadetNativeIsotherm` (see `src/chromsim/isotherms/mobile_phase_modulator.py`). Set
   `cadet_class` and implement `cadet_parameters(...)`, which returns `{cadet_attribute: value}`.
   Use `rate_constants()` for the `k_eq`/`k_kin` convention and `per_component()` for lists.
3. Copy the model equation and its reference from the
   [CADET documentation](https://cadet.github.io) into the docstring.

## Register, document, test

1. Import your module in `src/chromsim/isotherms/__init__.py` (Python isotherms:
   `from chromsim.isotherms.custom import my_isotherm`).
2. Add `docs/theory/isotherms/<name in snake_case>.md` with the equation, a parameter table and
   references, and link it in `mkdocs.yml` and `docs/theory/isotherms/index.md`.
3. Run `pixi run test-all`. The contract tests (`tests/plugins/test_isotherm_contract.py`) check
   your isotherm automatically:
   - docstring with equation and reference; units on every parameter; docs page present;
   - the `example` builds a valid project;
   - CADET-native: every mapped attribute exists, and a load–elute simulation conserves mass;
   - Python: `flux` has the right shape and zero rows for non-binding components, and a
     reference-backend simulation conserves mass.
4. Add at least one **physics check**: a special case where your model reduces to a known one
   (see `tests/physics/test_isotherm_equivalences.py`), or a comparison with published data.
