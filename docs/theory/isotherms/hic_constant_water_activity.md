# HIC, constant water activity

**Salt is component 0** (no bound state). Binding increases with salt through
β = β₀ exp(β₁ c_{p,0}):

$$
\frac{dq_i}{dt} = k_{a,i}\, c_{p,i}\left(1 - \sum_j \frac{q_j}{q_{max,j}}\right)^{\nu_i}
- k_{d,i}\, q_i\, 0.1^{\nu_i \beta}
$$

| Parameter | Where | Unit | Meaning |
|---|---|---|---|
| `k_eq` | per component | dimensionless | k_a/k_d |
| `q_max` | per component | mol/m³ solid | capacity |
| `nu` | per component | dimensionless | ligands per protein–ligand interaction |
| `k_kin` | per component | 1/s | desorption rate (optional) |
| `beta_0` | `isotherm.global_parameters` | dimensionless | ordered water at infinite dilution |
| `beta_1` | `isotherm.global_parameters` | m³/mol | salt dependence of ordered water |

```yaml
isotherm:
  model: HICConstantWaterActivity
  parameters:
    mAb: {k_eq: 20, q_max: 10, nu: 2}
  global_parameters: {beta_0: 0.1, beta_1: 1e-3}
```

**Check in the test suite:** with ν = 1 and β₀ = 0 it equals Langmuir with k_eq/q_max.

References: Jäpel & Buyel, *J. Chromatogr. A* 1679 (2022) 463408; CADET documentation.
