# 0009 — Optimizer interface

- Status: accepted
- Date: 2026-10-11

## Decision
SPEC §6 is implemented as specified:

- `Variable`, `Constraint` (g ≤ 0 feasible), `OptimizationProblem` (objectives and/or residuals,
  `n_objectives`, `parallel`) and `OptimizationResult` (best point, full history, Pareto set,
  residuals and Jacobian for least squares).
- `OptimizerBase` with `name`, `capabilities` and `optimize()`. A `ScaledProblem` handles scaling
  to [0, 1] (log10 where requested), history, callbacks, a cache of evaluated points, and batch
  evaluation in a spawn-context `ProcessPoolExecutor` (works on Windows and macOS).
- Adapters: `least_squares` (TRF, parallel forward-difference Jacobian), `nelder_mead`, `lbfgsb`
  (parallel central-difference gradient), `differential_evolution` (batch-parallel population,
  constraints), `nsga2` and `unsga3` (pymoo; multi-objective, constraints, parallel).
- Finite-difference steps are 1e-3 of each variable's range, deliberately well above simulation
  noise (CADET reltol 1e-6). This limits gradient-based accuracy on very curved analytic test
  functions (L-BFGS-B reaches f ≈ 1e-5 on Rosenbrock), which doesn't matter for noisy simulation
  objectives.
- `inverse_fit` builds an `OptimizationProblem` and accepts `optimizer=` and `parallel=`.
- Results are written to HDF5 with their full history (`write_optimization_result`).

## Consequences
- One interface serves parameter estimation and process optimization (Phase 6).
- Objectives must be picklable for parallel runs: module-level functions or objects.
- Resuming a stopped run is not automatic. Start again from the best point in the stored history.
