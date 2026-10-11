# Mixed-mode Langmuir (Python example)

A **Python isotherm**, solved by the reference backend (`solver: {backend: reference}`, LRM or
LRMP columns). It's the worked example for writing your own isotherm; see
[How to add an isotherm](../../how_to/add_isotherm.md).

$$
\frac{dq_i}{dt} = k_{a,i}\, e^{\gamma_i c_{p,0}}\, c_{p,i}\, q_{max,i}
\left(1 - \sum_{j\ge1} \frac{q_j}{q_{max,j}}\right) - k_{d,i}\, c_{p,0}^{\beta_i}\, q_i
$$

Hydrophobic interaction strengthens binding as salt rises (e^{γ c_s}). Electrostatic
interaction weakens it (c_s^β). Salt is component 0 and has no bound state.

| Parameter | Unit | Meaning |
|---|---|---|
| `k_eq` | m³/mol | k_a/k_d |
| `q_max` | mol/m³ solid | capacity |
| `beta` | dimensionless | electrostatic exponent β |
| `gamma` | m³/mol | hydrophobic coefficient γ |
| `k_kin` | 1/s | desorption rate k_d (required: Python isotherms are kinetic) |

**Check in the test suite.** These are the equations of CADET's built-in [MPM](mpm.md). The
Python version on the reference backend matches CADET's to 0.06 % NRMSE
(`tests/backends/test_reference_vs_cadet.py`).

References: Melander, El Rassi & Horváth, *J. Chromatogr.* 469 (1989) 3; Karlsson, Jakobsson,
Axelsson & Nilsson, *J. Chromatogr. A* 1055 (2004) 29.
