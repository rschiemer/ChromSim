# HIC, water on hydrophobic surfaces

Same parameters as [HIC constant water activity](hic_constant_water_activity.md); the desorption
term depends on the bound amount:

$$
\frac{dq_i}{dt} = k_{a,i}\, c_{p,i}\left(1 - \sum_j \frac{q_j}{q_{max,j}}\right)^{\nu_i}
- k_{d,i}\, q_i\left(\sum_j q_j\right)^{\nu_i \beta},
\qquad \beta = \beta_0\, e^{\beta_1 c_{p,0}}
$$

This is CADET's slightly modified multi-component form of the original single-component model.
Salt is component 0 and has no bound state. Concentrations are in mol/m³.

**Check in the test suite:** with ν = 1 and β₀ = 0 it equals Langmuir with k_eq/q_max.

References: Wang, Hahn & Hubbuch, *J. Chromatogr. A* 1465 (2016) 71; CADET documentation.
