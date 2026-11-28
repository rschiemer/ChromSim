# Isotherms

| Name in project file | Isotherm | Typical use |
|---|---|---|
| `Linear` | [Linear](linear.md) | low loading, tracers, Henry constants |
| `Langmuir` | [Langmuir](langmuir.md) | affinity resins, single- and multi-component saturation |
| `SMA` | [Steric mass action](sma.md) | ion exchange with salt gradients |
| `SAS` | [Self-association SMA](sas.md) | ion exchange with anti-Langmuirian (AL-L) elution |
| `GIEX` | [Generalized ion exchange](giex.md) | ion exchange with salt **and pH** dependence |
| `MPM` | [Mobile phase modulator Langmuir](mpm.md) | HIC; mixed-mode approximation |
| `HICConstantWaterActivity` | [HIC, constant water activity](hic_constant_water_activity.md) | HIC |
| `HICWaterOnHydrophobicSurfaces` | [HIC, water on hydrophobic surfaces](hic_water_on_hydrophobic_surfaces.md) | HIC |
| `BiLangmuir` | [Bi-Langmuir](bi_langmuir.md) | affinity resins with two site types |
| `MixedModeLangmuir` | [Mixed-mode Langmuir](mixed_mode_langmuir.md) | Python example; reference backend |
| `Colloidal` | [Colloidal (Xu & Lenhoff)](colloidal.md) | non-stoichiometric IEX, Protein A; CADET + reference |
| `CPA` | [Colloidal particle adsorption](cpa.md) | physically parameterised IEX (charge, radius); reference backend |

**Conventions shared by all isotherms**

- Every binding component has its own parameter set under `isotherm.parameters`. Components
  without one don't bind.
- Bound concentrations q are per volume of **solid phase** (mol/m³), as in CADET.
- `k_eq` is the equilibrium constant k_a/k_d. `k_kin` is the desorption rate constant k_d in
  1/s, so k_a = k_eq·k_kin. Without `k_kin` the isotherm is solved in rapid equilibrium. Give
  `k_kin` for all binding components or for none.
- Parameter values are plain numbers in the units listed on each page; they coincide with SI
  (mol/m³ = mM, m³/mol = 1/mM, 1/s).
- Isotherm-wide parameters (e.g. HIC β₀, β₁) go under `isotherm.global_parameters`.
- Modifiers come first: salt is component 0 for SMA, SAS, GIEX, MPM, HIC, Colloidal and CPA models; for GIEX the pH
  pseudo-component is component 1 (inlet unit `pH`).

Background and literature on the colloidal models: [Colloidal adsorption models](../colloidal_models.md).
Note that `Colloidal` uses CADET's meaning of `k_kin` (rate of the driving force c − c\*), not the
desorption rate ([ADR 0011](../../decisions/0011-colloidal-isotherms.md)).
