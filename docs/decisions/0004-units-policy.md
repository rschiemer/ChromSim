# 0004 — Units policy: user units at the boundary, SI inside

- Status: accepted
- Date: 2026-10-10

## Context
Scientists think in mL, CV, cm/h, mM and g/L. CADET works in SI units (m, s, mol/m³). Unit bugs
are among the most common and hardest-to-find errors in simulation code, and a lossless project
file has to keep the units the user wrote.

## Decision
- All unit handling lives in `chromsim/units.py`, built on one pint `UnitRegistry` that is
  created once and never modified.
- A `Quantity` stores the magnitude and the unit text **exactly as the user wrote them**. It is
  serialised back to the same text, so project files round-trip without loss. It is a pydantic
  model, not a pint object, so it stays easy to serialise and compare.
- Domain models check dimensions with annotated types (`Length`, `Volume`, `Duration`,
  `FlowRate`, `Concentration`, …) and expose SI values as properties with the unit in the name
  (`length_m`, `volume_m3`, `flow_rate_m3_per_s`).
- **Column volumes** (`CV`) are not a pint unit, because they depend on the column. They are
  resolved explicitly in `duration_to_si` using V_c = π·(d/2)²·L.
- **Linear velocity** means superficial velocity. Q = u·A with A = π·(d/2)².
- **Mass concentrations** become molar using the component's molecular weight. 1 kDa is taken
  as 1 kg/mol; pint's `Da` is a molecular mass, so that conversion is explicit.
- `UnitError` subclasses both `ChromSimError` and `ValueError`. pydantic therefore reports it
  as an ordinary validation error inside models, and callers can still catch it as a ChromSim
  error.

## Consequences
- Science code never sees user units. Any conversion outside `units.py` is a review finding.
- Project files keep the user's units. The same project written in different units (e.g.
  `0.2 m` and `20 cm`) produces different but equivalent files.
