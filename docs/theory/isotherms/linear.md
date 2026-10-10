# Linear isotherm

$$
\frac{dq_i}{dt} = k_{a,i}\, c_i - k_{d,i}\, q_i ,
\qquad q_i^* = k_{eq,i}\, c_i \ \text{at equilibrium}
$$

| Parameter | Unit | Meaning | Bounds |
|---|---|---|---|
| `k_eq` | dimensionless | Henry constant q*/c (q per solid volume) | ≥ 0 |
| `k_kin` | 1/s | desorption rate k_d (optional) | ≥ 0 |

```yaml
isotherm:
  model: Linear
  parameters:
    protein: {k_eq: 2.0}
```

**Check in the test suite:** the first moment of a pulse on a GRM column equals
t_pulse/2 + (ε_c V_c/Q)·[1 + F(ε_p + (1 − ε_p) k_eq)], with F = (1 − ε_c)/ε_c
(`tests/physics/test_column_physics.py`).

Reference: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
Chromatography*, 2nd ed., Academic Press 2006.
