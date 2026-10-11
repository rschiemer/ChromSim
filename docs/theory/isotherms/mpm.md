# Mobile phase modulator Langmuir (MPM)

Salt-dependent Langmuir isotherm for HIC, and an approximation for mixed-mode resins. **Salt is
component 0** and is inert (no bound state).

$$
\frac{dq_i}{dt} = k_{a,i}\, e^{\gamma_i c_{p,0}}\, c_{p,i}\, q_{max,i}
\left(1 - \sum_{j\ge1} \frac{q_j}{q_{max,j}}\right) - k_{d,i}\, c_{p,0}^{\beta_i}\, q_i
$$

The exponential term (γ, hydrophobicity) strengthens binding as salt increases (HIC-like). The
power law (β, ion-exchange characteristic) weakens it (IEX-like).

| Parameter | Unit | Meaning | Bounds |
|---|---|---|---|
| `k_eq` | m³/mol | k_a/k_d (salt enters as c_s^β with c_s in mol/m³) | ≥ 0 |
| `q_max` | mol/m³ solid | saturation capacity | ≥ 0 |
| `beta` | dimensionless | ion-exchange characteristic β | 0–50 |
| `gamma` | m³/mol (= 1/mM) | hydrophobicity γ | any |
| `k_kin` | 1/s | desorption rate k_d (optional) | ≥ 0 |

**Check in the test suite:** with β = γ = 0 the outlet equals the Langmuir isotherm's
(`tests/physics/test_isotherm_equivalences.py`).

References: Melander, El Rassi & Horváth, *J. Chromatogr.* 469 (1989) 3; Karlsson, Jakobsson,
Axelsson & Nilsson, *J. Chromatogr. A* 1055 (2004) 29; CADET documentation.
