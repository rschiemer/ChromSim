# 0008 — Reference backend and the route to production for custom isotherms

- Status: accepted
- Date: 2026-10-11

## Context
SPEC §5 asks for isotherms written in Python, for models CADET doesn't implement (e.g.
literature mixed-mode models). CADET-Core is C++ and can't call Python binding functions. A second,
transparent solver also lets us cross-check CADET and is useful for teaching.

## Decision
- A **reference backend** in pure numpy/scipy for **LRM and LRMP** columns, with mixer and
  tubing, the same flow sheet, method handling and initial state as the CADET backend. It has
  one file per operator (convection, dispersion, film transfer, binding) and a solver file.
- **Numerics.** Finite volumes with a second-order van Leer (TVD) reconstruction, closed
  Danckwerts boundaries, and scipy BDF per method step with a sparse Jacobian pattern. Only
  states that influence the equations are integrated: no bound phase for non-binding
  components, and SMA uses the bound-salt state for the free sites. Dead states made scipy's
  finite-difference Jacobian overflow.
- **Kinetic binding only.** Python isotherms provide `flux(c, q, p)`. A large k_kin approximates
  rapid equilibrium.
- **Acceptance** (tests/backends/test_reference_vs_cadet.py): ≤ 1 % NRMSE against CADET on
  shared cases (measured ≤ 0.4 %), and mass balance ≤ 0.1 %. The Python example isotherm
  MixedModeLangmuir matches CADET's MPM to 0.06 %.

## Route to production (not implemented in v1)
The reference backend is for development and validation, not for large studies. A Python isotherm
that proves useful has two ways to production speed:

1. **CADET-Core 6 data-driven binding** (spline/GPR/NN surrogates), once 6.0 is final. Fit the
   surrogate to the Python model's equilibrium data, and record the switch in an ADR together with
   the CADET-Core 6 upgrade.
2. **Upstreaming a C++ binding model** to CADET-Core, with the Python version as its reference
   implementation and the cross-check test as acceptance.

## Consequences
- GRM columns, rapid-equilibrium binding and CADET-only isotherms (GIEX, HIC, Bi-Langmuir) are
  refused by the reference backend with an actionable message.
- The reference backend is 10–100× slower than CADET, but readable in an afternoon.
