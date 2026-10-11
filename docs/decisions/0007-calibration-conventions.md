# 0007 — Calibration conventions

- Status: accepted
- Date: 2026-10-11

## Context
Calibration results depend on choices that are easy to make silently: how residuals from
different data are weighted, how values below the LOQ are treated, how uncertainties are
computed, and which tracer quantities are trusted. These choices are recorded here.

## Decisions
1. **Volume axis.** All calibration works on pumped volume (m³), like the importers and results.
2. **Tracer moments.** Moments are used for volumes and porosities (first moments). D_ax from
   variances is flagged as unreliable when the column contributes less than 50 % of the measured
   variance; the documented remedy is an inverse fit of the pulse. The system volume measured at
   the UV detector includes `uv_delay_volume`, which `apply_to` subtracts.
3. **Yamamoto relation.** It is derived for ChromSim's model, with GH = g(1 − ε_t)V_c and
   I_R ≫ c_0 (docs/theory/yamamoto.md). It is validated against CADET to < 1 %.
4. **Residuals.** Each target returns residuals divided by √n. Signals are normalised by
   max|y_measured|; assays are weighted by the replicate SD. A below-LOQ value contributes 0, or
   (sim − LOQ)/SD when the simulation exceeds the LOQ (SD = LOQ/2 if none is given). The residual
   vector always has the same length.
5. **Search space.** Each parameter is mapped to [0, 1], linearly or in log10. The finite
   difference step is 1e-3 of the range.
6. **Uncertainty.** Linearised covariance s²(JᵀJ)⁻¹ with s² = SSE/(m − p), and Student-t
   intervals (computed in log space for log parameters). |R| > 0.95 and parameters at a bound
   are flagged. Documentation and Tutorial 03 state that the intervals exclude model error.
7. **Optimizer.** scipy `least_squares` (TRF) is called directly until the optimizer interface
   exists (Phase 5).

## Consequences
- Fits of different data types are comparable and reproducible.
- A profile-likelihood analysis (SPEC §7.3, optional) is not implemented yet; it would quantify
  the non-linear k_eq–ν valley better than the linearised intervals.
