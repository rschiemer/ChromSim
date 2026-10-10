# 0005 — Conventions of the CADET backend

- Status: accepted
- Date: 2026-10-10

## Context
To turn a ChromSim project into a CADET simulation, the project has to be mapped onto
CADET-Process objects. Several choices in that mapping change results or what parameters mean,
so they are recorded here.

## Decisions
1. **Kinetic parameters.** `k_eq` = k_a/k_d and `k_kin` = k_d (1/s), so k_a = k_eq·k_kin. Without
   `k_kin`, binding is in rapid equilibrium and ChromSim passes k_a = k_eq and k_d = 1. This
   matches SPEC §5, where `k_kin` maps to CADET's `desorption_rate`.
2. **SMA reference concentrations.** c_ref = q_ref = Λ (the ionic capacity). k_eq keeps its
   classical, dimensionless meaning, and k_d becomes a well-scaled rate in 1/s.
3. **Initial state.** The column starts clean, equilibrated with the first step's mobile-phase
   modifiers (salt for SMA). Everything else starts at zero, and for SMA q_salt = Λ. This lets
   tracer and breakthrough methods start at t = 0, and it avoids a salt-free column (with SMA,
   binding diverges at zero salt).
4. **System.** The flow sheet is inlet → mixer (CSTR) → ½ tubing (DPFR) → column → ½ tubing
   (DPFR) → outlet. The tubing gets two new optional `System` fields with defaults:
   `tubing_diameter` (0.5 mm, typical ÄKTA) and `tubing_axial_dispersion` (1e-6 m²/s,
   placeholder). These are additive, so the project schema stays at version 1. The physics test
   checks that system volumes delay the first moment by exactly V/Q. Detector delays are applied
   when UV/conductivity are predicted (Phase 3).
5. **Output grid.** CADET-Process's default 1 s output resolution is used.
6. **Temporary files.** CADET always runs in its own temporary folder (see ADR 0001).

## Consequences
- Parameter meanings are documented on each isotherm page (`docs/theory/isotherms/`).
- Tubing dispersion defaults only matter for very small columns. Calibrate them with a system
  tracer pulse (U3).
- Projects whose first step loads protein onto a column that is *not* clean cannot be expressed
  yet. Add an explicit initial-state section if this is ever needed.
