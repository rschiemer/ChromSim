# Generalized ion exchange (GIEX)

SMA with salt and **pH** dependence. **Component 0 is salt, component 1 the pH
pseudo-component.** Give its inlet values in the unit `pH`, e.g. `inlet: {salt: 50 mM, pH: 5.5 pH}`.

CADET's model, for binding components i ≥ 2:

$$
\frac{dq_i}{dt} = k_{a,i}\, c_{p,i}\left(\frac{\bar q_0}{q_{ref}}\right)^{\nu_i(pH)}
- k_{d,i}\, q_i\left(\frac{c_{p,0}}{c_{ref}}\right)^{\nu_i(pH)},
\qquad
\nu_i(pH) = \nu_i + \nu_{i,lin}\, pH + \nu_{i,quad}\, pH^2
$$

with $\bar q_0 = \Lambda - \sum_j (\nu_j + \sigma_j) q_j$ and c_ref = q_ref = Λ, as for SMA.

**What ChromSim exposes.** CADET lets k_a and k_d each depend exponentially on pH, salt and
protein. ChromSim exposes the dependence of the *equilibrium constant*, which is what
equilibrium data can identify:

$$
\frac{k_a}{k_d} = k_{eq}\exp\!\left(k_{eq,lin}\, pH + k_{eq,quad}\, pH^2 + k_{eq,salt}\,\frac{c_{p,0}}{\Lambda}\right),
\qquad k_d = k_{kin}.
$$

| Parameter | Unit | Meaning | Required |
|---|---|---|---|
| `k_eq` | dimensionless | equilibrium constant at pH 0 (no salt term) | yes |
| `nu`, `nu_lin`, `nu_quad` | –, 1/pH, 1/pH² | characteristic charge and its pH dependence | `nu` |
| `sigma` | dimensionless | steric factor | yes |
| `k_eq_lin`, `k_eq_quad` | 1/pH, 1/pH² | pH dependence of ln k_eq | no |
| `k_eq_salt` | dimensionless | salt dependence of ln k_eq (per c_s/Λ) | no |
| `k_kin` | 1/s | desorption rate (optional; rapid equilibrium otherwise) | no |

**Check in the test suite:** at constant pH, GIEX equals SMA with ν = ν + ν_lin·pH and
k_eq = k_eq·exp(k_eq_lin·pH).

References: Huuk, Hahn, Doninger, Griesbach, Hepbildikler & Hubbuch, *Biotechnol. J.* 12 (2017)
1600336; Mollerup, *Chem. Eng. Technol.* 31 (2008) 864; CADET documentation.
