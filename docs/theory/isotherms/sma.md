# Steric mass action (SMA)

For ion exchange. **Salt must be component 0**; it has no parameters of its own.

$$
\frac{dq_i}{dt} = k_{a,i}\, c_i \left(\frac{\bar q_0}{q_{ref}}\right)^{\nu_i}
 - k_{d,i}\, q_i \left(\frac{c_0}{c_{ref}}\right)^{\nu_i},
\qquad
\bar q_0 = \Lambda - \sum_{j \ge 1} (\nu_j + \sigma_j)\, q_j
$$

Here c_0 is the salt concentration, Λ the ionic capacity of the resin (`resin.ionic_capacity`),
ν the characteristic charge and σ the steric factor. Bound salt follows from electroneutrality,
q_0 = Λ − Σ ν_j q_j.

**Reference concentrations.** ChromSim sets c_ref = q_ref = Λ. The equilibrium constant then keeps
its classical, dimensionless meaning,

$$
k_{eq,i} = \frac{k_{a,i}}{k_{d,i}} = \frac{q_i\, c_0^{\nu_i}}{c_i\, \bar q_0^{\nu_i}},
$$

and `k_kin` = k_d is a rate in 1/s. Without the reference concentrations, k_d would have units
of (mM)^−ν·s^−1 and span many orders of magnitude.

| Parameter | Unit | Meaning | Bounds |
|---|---|---|---|
| `k_eq` | dimensionless | equilibrium constant | ≥ 0 |
| `nu` | dimensionless | characteristic charge ν | 0–50 |
| `sigma` | dimensionless | steric factor σ | 0–500 |
| `k_kin` | 1/s | desorption rate k_d (optional) | ≥ 0 |

```yaml
resin: {..., ionic_capacity: 1200 mM}
isotherm:
  model: SMA
  parameters:
    mAb: {k_eq: 1.0e-2, nu: 6.0, sigma: 30.0, k_kin: 1.0e+2}
```

**Initial state.** The column starts equilibrated with the first step's salt, with every site
carrying a salt counter-ion (q_0 = Λ) and no protein bound.

References: Brooks & Cramer, AIChE J. 38 (1992) 1969. Leweke & von Lieres, Comput. Chem. Eng.
113 (2018) 274 (CADET).
