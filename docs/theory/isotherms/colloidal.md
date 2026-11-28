# Colloidal isotherm (Xu & Lenhoff)

`isotherm: {model: Colloidal}` · **CADET-native** (`MultiComponentColloidal`) — also runs on the
reference backend (LRM, LRMP), where it reproduces CADET to ≤ 0.3 % NRMSE.

Proteins are charged colloidal spheres. At equilibrium the adsorbed proteins form a hexagonal
lattice whose spacing R shrinks as loading increases. A protein–surface term K_e sets the
low-loading (Henry) behaviour; a screened Yukawa repulsion between neighbours (constant b_pp,
range 1/κ) makes the isotherm bend over without a fixed capacity. See
[Colloidal adsorption models](../colloidal_models.md) for the background and how it compares to
CPA.

$$
\frac{dq_i}{dt} = k_{kin,i}\left(c_{p,i} - c_{p,i}^\star\right), \qquad
c_{p,i}^\star = q_i \exp\left[\frac{n(3+\kappa R)}{4 q_{tot}}\sum_j q_j\sqrt{b_{pp,i}b_{pp,j}}
\,\frac{r_i+r_j}{2R}\,e^{-\kappa(R-r_i-r_j)} - \ln K_{e,i}\right]
$$

$$
R = \sqrt{\frac{2\phi}{\sqrt{3}\,N_A\,q_{tot}}}, \qquad
\kappa = \frac{1}{\kappa_f\,c_{p,0}^{\kappa_{ef}} + \kappa_c}
$$

$$
\ln K_{e,i} = k_{a,i}\,c_{p,0}^{-k_{b,i}} + k_{c,i}\,e^{k_{d,i}c_{p,0}}, \qquad
b_{pp,i} = b_{a,i}\,c_{p,0}^{b_{b,i}} + b_{c,i}\,e^{b_{d,i}c_{p,0}}
$$

Below `linear_threshold` (q_tot ≤ c_ε) the lattice term is dropped: c*_i = q_i/K_e,i.
Salt is component 0 and has no bound state.

## Parameters

Per binding component (`isotherm.parameters`):

| Parameter | Symbol | Unit | Meaning |
|---|---|---|---|
| `radius` | r_i | m | protein radius (mAb ≈ 4.5e-9–5.5e-9) |
| `ln_ke_salt_factor` | k_a | – | factor of the salt power term of ln K_e |
| `ln_ke_salt_exponent` | k_b | – | exponent of that term (ln K_e ∝ c_s^(−k_b)) |
| `ln_ke_exp_factor` | k_c | – | factor of the exponential term (optional, default 0) |
| `ln_ke_exp_salt_coefficient` | k_d | m³/mol | salt coefficient in it (optional, default 0) |
| `b_pp_salt_factor` | b_a | – | factor of the salt power term of b_pp |
| `b_pp_salt_exponent` | b_b | – | exponent of that term |
| `b_pp_exp_factor` | b_c | – | optional, default 0 |
| `b_pp_exp_salt_coefficient` | b_d | m³/mol | optional, default 0 |
| `k_kin` | k_kin | 1/s | rate of the driving force c − c* (omit for rapid equilibrium) |

Whole isotherm (`isotherm.global_parameters`):

| Parameter | Symbol | Unit | Meaning |
|---|---|---|---|
| `phase_ratio` | φ | m²/m³ | adsorber surface per solid volume |
| `kappa_factor` | κ_f | m·(m³/mol)^κ_ef | salt-dependent part of the Debye length |
| `kappa_exponent` | κ_ef | – | exponent; **must be ≥ 0 in CADET-Process 0.12** (below) |
| `kappa_constant` | κ_c | m | constant part of the Debye length |
| `coordination_number` | n | – | neighbours in the lattice (default 6, hexagonal) |
| `linear_threshold` | c_ε | mol/m³ | switch to the linear form (default 1e-8) |

```yaml
isotherm:
  model: Colloidal
  parameters:
    mAb: {radius: 4.5e-9, ln_ke_salt_factor: 60, ln_ke_salt_exponent: 0.5,
          b_pp_salt_factor: 10, b_pp_salt_exponent: 0, k_kin: 1.0}
  global_parameters: {phase_ratio: 4.0e+7, kappa_factor: 0, kappa_exponent: 0,
                      kappa_constant: 1.0e-9}
```

## Things to know

- **k_kin is not a desorption rate here.** It multiplies the liquid-phase driving force
  c − c*, exactly as in CADET. In the linear limit the desorption rate is k_kin/K_e. This
  differs from the `k_eq`/`k_kin` convention of the other ChromSim isotherms
  ([ADR 0011](../../decisions/0011-colloidal-isotherms.md)).
- **Screening and salt.** Physically the Debye length falls as c_s^(−1/2), which needs
  κ_ef = −0.5. CADET-Process 0.12 declares `kappa_exponential` non-negative and refuses
  negative values, so ChromSim does too. Within one salt window, use a constant κ
  (κ_f = 0, κ_c = Debye length at the elution salt) and let the salt dependence of b_pp absorb
  the rest, or use the [CPA model](cpa.md), which computes κ from the ionic strength.
- **pH** (CADET's `COL_USE_PH`) is not exposed yet; ChromSim passes `use_ph: false`.
- **Parameter names.** CADET's interface docs describe `COL_LOGKEQ_SALT_POWEXP` as a factor and
  `…_POWFACT` as an exponent; the C++ source does the opposite. ChromSim follows the source,
  and `tests/backends/test_reference_vs_cadet.py` pins it with all terms switched on.

Tutorial 06 (`examples/06_colloidal_models.py`) calibrates this model to CPA and compares both
at low and high load.

## Checks in the test suite

- Python flux = the documented equation, term by term (`tests/unit/test_colloidal_isotherms.py`).
- b_pp = 0 and k_b = 0 give exactly the Linear isotherm with k_eq = e^{k_a}
  (`tests/physics/test_isotherm_equivalences.py`).
- Reference backend vs CADET, two proteins, overloaded enough that the lattice term changes the
  outlet by 16–29 %: ≤ 0.3 % NRMSE.

## References

- Xu, X. & Lenhoff, A. M. Binary adsorption of globular proteins on ion-exchange media.
  *J. Chromatogr. A* 1216 (2009) 6177–6195. doi:10.1016/j.chroma.2009.06.082
- Oberholzer, M. R. & Lenhoff, A. M. Protein adsorption isotherms through colloidal energetics.
  *Langmuir* 15 (1999) 3905–3914. doi:10.1021/la981199k
- Kumar, V. et al. (Lenhoff group) Robust mechanistic modeling of protein ion-exchange
  chromatography. *J. Chromatogr. A* 1660 (2021) 462669.
- Bhoyar, S. et al. Predictive mechanistic modeling of loading and elution in protein A
  chromatography. *J. Chromatogr. A* (2024) 464558. doi:10.1016/j.chroma.2023.464558
- CADET-Core documentation, *Multi Component Colloidal* (binding model and interface pages),
  and `src/libcadet/model/binding/ColloidalBinding.cpp` (v5.1).
