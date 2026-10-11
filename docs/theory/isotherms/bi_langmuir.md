# Bi-Langmuir

Two independent types of binding site (j = 1, 2). Components compete only for sites of the same
type. Typical for affinity resins with high- and low-affinity sites.

$$
\frac{dq_{i,j}}{dt} = k_{a,i}^{(j)}\, c_{p,i}\, q_{max,i}^{(j)}
\left(1 - \sum_k \frac{q_{k,j}}{q_{max,k}^{(j)}}\right) - k_{d,i}^{(j)}\, q_{i,j}
$$

| Parameter | Unit | Meaning |
|---|---|---|
| `k_eq_1`, `k_eq_2` | m³/mol | k_a/k_d of site 1 and 2 |
| `q_max_1`, `q_max_2` | mol/m³ solid | capacity of site 1 and 2 |
| `k_kin_1`, `k_kin_2` | 1/s | desorption rates (optional, both or none) |

**Check in the test suite:** two identical sites behave like one Langmuir site with double
capacity. This also checks how parameters are ordered when passed to CADET.

References: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
Chromatography*, 2nd ed., Academic Press 2006; CADET documentation.
