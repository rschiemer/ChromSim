# Colloidal particle adsorption (CPA)

`isotherm: {model: CPA}` · **Python isotherm** — reference backend (`solver: {backend:
reference}`, LRM or LRMP).

The CPA model of Briskot et al. (2020, 2021) describes ion exchange without a stoichiometric
exchange of counter-ions. The protein is a charged hard sphere, the adsorber a charged plane,
and binding follows from double-layer physics. Its parameters have a physical meaning (radius,
net charge), so they can be checked against independent data, and it tends to reproduce
high-load isotherms better than SMA, which overestimates nonlinearity (Briskot et al. 2021).
See [Colloidal adsorption models](../colloidal_models.md) for the background.

## Equations

Kinetic form (Briskot et al. 2021b, eqs. 12–19):

$$
\frac{dq_i}{dt} = k_{kin,i}\left(K_{H,i}\,B_i(\Theta)\,e^{-u_{lat,i}/k_BT}\,c_{p,i} - q_i\right),
\qquad
K_{H,i} = \Delta_i\,\frac{e^{y_i}-1}{y_i}
$$

**1 · Screening.** Ionic strength I = c_{p,0} + I_0 (1:1 salt plus buffer background):

$$\kappa = \sqrt{2e^2N_AI/(\varepsilon_r\varepsilon_0k_BT)}$$

**2 · Surface potentials** from the Grahame equation, for the adsorber
(σ_A = z_L e N_A Γ_L) and the protein (σ_i = Z_i e / 4πa_i²):

$$\psi = \frac{2k_BT}{e}\operatorname{arsinh}\frac{e\,\sigma}{2\varepsilon_r\varepsilon_0\kappa k_BT}$$

**3 · Protein–adsorber well.** Sphere–plane double-layer energy at constant charge (Derjaguin
approximation, linearised Poisson–Boltzmann), x = e^{−κz}:

$$u_A(z) = \pi\varepsilon_r\varepsilon_0a_i\left[2\psi_A\psi_i\ln\frac{1+x}{1-x}
-(\psi_A^2+\psi_i^2)\ln(1-x^2)\right]$$

For opposite charges the minimum is at x* = −2ψ_Aψ_i/(ψ_A²+ψ_i²) (closed form, no
optimisation). y_i = −u_A(δ_m)/k_BT is its depth; like charges give y_i = 0, i.e. K_H = Δ_i
(partitioning only). The well-width integral is idealised as linear (Briskot 2021a, eq. 10),
which gives the factor (e^y − 1)/y.

**4 · Crowding.** With Γ_j = q_j/A_s and Θ = πN_A Σ a_j²Γ_j, the hard-disc available surface
function from scaled particle theory:

$$\ln B_i = \ln(1-\Theta) - \frac{2\pi a_i\rho_1 + \pi a_i^2\rho_0}{1-\Theta}
- \frac{(\pi a_i\rho_1)^2}{(1-\Theta)^2},\quad \rho_0 = N_A\sum_j\Gamma_j,\;
\rho_1 = N_A\sum_j a_j\Gamma_j$$

and the screened repulsion of neighbours on a hexagonal lattice of spacing
D = (2/(√3 N_A ΣΓ_j))^{1/2}:

$$u_{lat,i} = 3\sqrt3\,D\,N_A\,\frac{e^{-\kappa D}}{1-e^{-\frac{3\sqrt3}{2\pi}\kappa D}}
\sum_j\Gamma_j\beta_{ij},\qquad
\beta_{ij} = \frac{Z_{lat,i}Z_{lat,j}e^2}{4\pi\varepsilon_r\varepsilon_0}
\frac{e^{\kappa(a_i+a_j)}}{(1+\kappa a_i)(1+\kappa a_j)}$$

Each step is a public function in `chromsim.isotherms.custom.colloidal_particle_adsorption`
(`debye_kappa`, `grahame_potential`, `well_depth`, `henry_factor`,
`available_surface_function`, `lateral_energy`), and `CPA.equilibrium_factor(c_salt, q, p)`
gives q*/c, e.g. to plot K_H against salt.

## Parameters

Per binding component (`isotherm.parameters`):

| Parameter | Symbol | Unit | Meaning |
|---|---|---|---|
| `radius` | a_i | m | protein radius (mAb ≈ 5e-9) |
| `charge` | Z_i | – | net charge at the operating pH; opposite sign to the ligand to bind |
| `Delta` | Δ_i | – | A_s(δ* − δ_m): adsorption-layer width × surface; scales K_H |
| `lateral_charge` | Z_lat,i | – | effective charge for protein–protein repulsion (optional, 0 = off) |
| `k_kin` | k_kin | 1/s | desorption rate constant (required: Python isotherms are kinetic) |

Whole isotherm (`isotherm.global_parameters`):

| Parameter | Symbol | Unit | Meaning |
|---|---|---|---|
| `ligand_density` | Γ_L | mol/m² | charged ligands per surface (≈ ionic capacity / A_s; 1e-6 = 0.6 nm⁻²) |
| `ligand_charge` | z_L | – | −1 cation exchanger, +1 anion exchanger |
| `specific_surface` | A_s | m²/m³ | adsorber surface per solid volume (Poros XS ≈ 2.4e8) |
| `background_ionic_strength` | I_0 | mol/m³ | buffer ionic strength without salt (default 0) |
| `temperature` | T | K | default 298.15 |
| `relative_permittivity` | ε_r | – | default 78.4 |

```yaml
isotherm:
  model: CPA
  parameters:
    mAb: {radius: 5.0e-9, charge: 30, Delta: 0.5, lateral_charge: 5, k_kin: 1.0}
  global_parameters: {ligand_density: 2.0e-6, ligand_charge: -1, specific_surface: 2.0e+8}
solver: {backend: reference}
```

**Correlations.** Δ and Z both shift retention; Z sets the slope of log K_H against salt (it
plays the role of SMA's ν), Δ its offset. Fit them from LGE runs at low load first, then
Z_lat (and A_s if not measured) from breakthrough or high-load elution. Lorenz-Cristea et al.
(2024) describe a systematic calibration route.

## Simplifications against the publications

ChromSim's CPA is the published model with these deliberate simplifications; each can be
lifted later without changing the parameter names above.

1. **No charge regulation.** The ligand is taken as fully dissociated (true for strong
   exchangers, e.g. sulfopropyl pK ≈ 2.3, quaternary amine) and Z_i fixed. Briskot et al.
   (2020) regulate both surfaces.
2. **No pH dependence.** Briskot et al. (2021a, b) use polynomial Z_i(pH) and log Δ_i(pH).
   Simulate one pH per project for now.
3. **Constant k_kin.** Briskot et al. (2021b, eq. 18) let k_kin depend on the well depth.
4. **One A_s** for all components (the papers allow a size-dependent accessible surface).
5. **Grahame for both surfaces.** Equivalent to Debye–Hückel for weakly charged proteins.

Tutorial 06 (`examples/06_colloidal_models.py`) walks through every step on synthetic data and
compares CPA with the Colloidal isotherm at low and high load.

## Checks in the test suite

`tests/unit/test_colloidal_isotherms.py`:

- Debye length 0.96 nm at 100 mM; Grahame → Debye–Hückel at low charge.
- Closed-form well depth = numerical minimum of the energy profile (three potential pairs).
- B(Θ) = SPT (1−Θ)exp[−3Θ/(1−Θ) − Θ²/(1−Θ)²] for one component, 1 − 4Θ at low coverage,
  identical species add up, larger proteins are excluded more.
- u_lat → six nearest neighbours for κD ≫ 1.
- K_H falls monotonically with salt; like charges only partition; q*(c) monotonic and saturating.
- Plugin contract: reference-backend load–elute run conserves mass.

## References

- Briskot, T., Hahn, T., Huuk, T. & Hubbuch, J. Protein adsorption on ion exchange adsorbers:
  A comparison of a stoichiometric and non-stoichiometric modeling approach. *J. Chromatogr. A*
  1653 (2021) 462397. doi:10.1016/j.chroma.2021.462397 (2021a)
- Briskot, T., Hahn, T., Huuk, T., Wang, G., Kluters, S., Studts, J., Wittkopp, F., Winderl,
  J., Schwan, P., Hagemann, I., Kaiser, K., Trapp, A., Stamm, S. M., Koehn, J., Malmquist, G. &
  Hubbuch, J. Analysis of complex protein elution behavior in preparative ion exchange processes
  using a colloidal particle adsorption model. *J. Chromatogr. A* 1654 (2021) 462439. (2021b)
- Briskot, T., Hahn, T., Huuk, T. & Hubbuch, J. Adsorption of colloidal proteins in
  ion-exchange chromatography under consideration of charge regulation. *J. Chromatogr. A* 1611
  (2020) 460608. doi:10.1016/j.chroma.2019.460608
- Briskot, T. *Science and Risk-Based Development of Ion Exchange Chromatography.* PhD thesis,
  KIT (2022), open access: publikationen.bibliothek.kit.edu/1000144198 — chapters 4–6 contain
  the three papers with full derivations.
- Helfand, E., Frisch, H. L. & Lebowitz, J. L. Theory of the two- and one-dimensional rigid
  sphere fluids. *J. Chem. Phys.* 34 (1961) 1037. (scaled particle theory, hard discs)
- Usui, S. Interaction of electrical double layers at constant surface charge. *J. Colloid
  Interface Sci.* 44 (1973) 107; Hogg, R., Healy, T. W. & Fuerstenau, D. W. *Trans. Faraday
  Soc.* 62 (1966) 1638.
- Lorenz-Cristea, O. et al. A systematic approach for estimating colloidal particle adsorption
  model parameters. *J. Chromatogr. A* (2024) 465512. doi:10.1016/j.chroma.2024.465512
