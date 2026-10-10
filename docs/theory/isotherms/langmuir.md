# Langmuir isotherm

Multi-component competitive Langmuir, as implemented in CADET:

$$
\frac{dq_i}{dt} = k_{a,i}\, c_i\, q_{max,i}\left(1 - \sum_j \frac{q_j}{q_{max,j}}\right) - k_{d,i}\, q_i
$$

At equilibrium, for a single component:

$$
q^* = \frac{q_{max}\, k_{eq}\, c}{1 + k_{eq}\, c}
$$

| Parameter | Unit | Meaning | Bounds |
|---|---|---|---|
| `k_eq` | m³/mol (= 1/mM) | equilibrium constant k_a/k_d | ≥ 0 |
| `q_max` | mol/m³ (= mM) solid | saturation capacity | ≥ 0 |
| `k_kin` | 1/s | desorption rate k_d (optional) | ≥ 0 |

```yaml
isotherm:
  model: Langmuir
  parameters:
    protein: {k_eq: 1.0, q_max: 10.0}
```

**Check in the test suite:** a breakthrough curve reaches the feed concentration, and its
stoichiometric time is (V_c/Q)·[ε_t + φ·q*/c_feed], where φ is the solid fraction of the column.

References: Langmuir, J. Am. Chem. Soc. 40 (1918) 1361. Leweke & von Lieres, Comput. Chem. Eng.
113 (2018) 274 (CADET).
