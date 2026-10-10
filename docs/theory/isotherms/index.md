# Isotherms

| Name in project file | Isotherm | Typical use |
|---|---|---|
| `Linear` | [Linear](linear.md) | low loading, tracers, Henry constants |
| `Langmuir` | [Langmuir](langmuir.md) | affinity resins, single- and multi-component saturation |
| `SMA` | [Steric mass action](sma.md) | ion exchange with salt gradients |

**Conventions shared by all isotherms**

- Every binding component has its own parameter set under `isotherm.parameters`. Components
  without one don't bind.
- Bound concentrations q are per volume of **solid phase** (mol/m³), as in CADET.
- `k_eq` is the equilibrium constant k_a/k_d. `k_kin` is the desorption rate constant k_d in
  1/s, so k_a = k_eq·k_kin. Without `k_kin` the isotherm is solved in rapid equilibrium. Give
  `k_kin` for all binding components or for none.
- Parameter values are plain numbers in the units listed on each page. All of them coincide with
  SI (mol/m³ = mM, m³/mol = 1/mM, 1/s).
