# Calibrating a model

Calibration usually goes in four steps; Tutorial 03 (`examples/03_calibrate.py`) runs all of
them on example data.

| Step | Gives | Function |
|---|---|---|
| Tracer pulses | system volume, ε_c, ε_t, ε_p (moments); D_ax, mixer volume (fit) | `analyse_tracers`, then `inverse_fit` |
| Low-load gradient runs | starting values for SMA ν, k_eq | `yamamoto_fit` |
| Preparative runs | refined isotherm (and other) parameters | `inverse_fit` |
| Uncertainty | confidence intervals, correlations | `parameter_uncertainty` |

## 1. Tracer pulses

Give each run as a `TracerRun`: a baseline-free signal on the volume axis, the injection start
and volume, and its kind:

- `system`: no column (bypass). Gives the system volume up to the detector.
- `non_penetrating`: a pore-excluded tracer such as large dextran. Gives ε_c.
- `pore_penetrating`: acetone, or a salt pulse at high salt. Gives ε_t.

ε_p follows from ε_t = ε_c + (1 − ε_c)ε_p. `analysis.apply_to(project)` writes the results into a
project. The system volume measured at the UV detector includes `system.uv_delay_volume`, which is
subtracted before the remainder is split into mixer and tubing.

!!! warning "Variances on small columns"
    First moments (volumes, porosities) are robust. Variances are not: on small columns the
    mixer broadens the peak more than the column, and baseline noise in long tails dominates.
    `analysis.notes` warns when the column contributes less than half of the measured variance.
    Then fit D_ax (and the mixer volume) with `inverse_fit` on the dextran pulse, as in
    Tutorial 03.

## 2. Yamamoto method

For each run you need the gradient slope g (`gradient_slope(step)`) and the outlet salt at the
peak maximum, I_R (`elution_salt_at_peak(uv, salt, uv_delay, conductivity_delay)`). Then
`yamamoto_fit(runs, V_c, ε_t, Λ)` returns ν, k_eq and the regression, and `plot_yamamoto` shows
it. The derivation is on the [theory page](../theory/yamamoto.md).

Without a conductivity trace, take the outlet salt from the method:
`inlet_profile(project, "salt", delay_volume_m3=V_sys + ε_t·V_c)` shifts the programmed
inlet salt by the volume the salt front needs to reach the detector (case study 02).

## 3. Inverse fit

```python
experiment = api.FitExperiment("run01", project, targets=[
    api.SignalTarget(uv, window=(60e-6, 140e-6)),          # UV in AU, baseline-corrected
    api.AssayTarget(sec, fractions=run.fractions),          # SD-weighted
])
# A measured concentration (e.g. from another detector or a published data set):
#   api.SignalTarget(lysozyme_mM, component="lysozyme")
parameters = [
    api.FitParameter(path="isotherm.parameters.mAb.k_eq", lower=1e-3, upper=1e-1, scale="log"),
    api.FitParameter(path="isotherm.parameters.mAb.nu", lower=3, upper=10),
    api.FitParameter(path="method.1.inlet.mAb", lower=0.5, upper=3, unit="g/L", per_experiment=True),
]
fit = api.inverse_fit([experiment, ...], parameters, multi_start=2)
```

- **Paths** point into the project file (`column.film_diffusion.mAb`, `method.1.inlet.mAb`, …).
  Quantities need a `unit`.
- **Shared or per-experiment.** By default a parameter is shared by all experiments. Set
  `per_experiment=True` for things like the load concentration of each run.
- **Targets.** Each target's residuals are divided by √n, so a long trace and a few assay values
  weigh comparably; `weight` changes that. Signals use `normalized` residuals and assays use
  `sd_weighted`. Below-LOQ values only add a penalty when the simulation exceeds the LOQ.
- **Scaling.** Use `scale="log"` for rate and equilibrium constants. The search runs in [0, 1]
  for every parameter.
- **`multi_start`** adds random starting points, in case the fit gets stuck in a local minimum.

## 4. Uncertainty

`parameter_uncertainty(fit)` gives standard errors, 95 % intervals (in log space for log
parameters), the correlation matrix, and flags for |R| > 0.95 and for parameters at a bound.
`plot_correlation(report)` draws the matrix.

!!! warning "What the intervals mean"
    The intervals are linearised and cover **measurement noise only**, not model error. For
    strongly correlated parameters (SMA k_eq and ν typically have |R| > 0.99), even a tiny model
    error, such as a coarser grid than the "truth", moves the estimates far along the correlation
    valley while the fit to the data stays excellent. Tutorial 03 shows this. Report correlated
    parameters together, fit on the grid you will predict with, and add experiments that separate
    the parameters (several gradient slopes, isocratic runs).

## Optimizers and parallel fits

`inverse_fit(..., optimizer="least_squares", parallel=1)`:

| Optimizer | Kind | Notes |
|---|---|---|
| `least_squares` (default) | local, least squares | trust-region-reflective; Jacobian columns in parallel |
| `nelder_mead` | local, derivative-free | robust for small, noisy problems |
| `lbfgsb` | local, gradient | central-difference gradient in parallel |
| `differential_evolution` | global | many simulations; use with `parallel` |

With `parallel=N`, simultaneous simulations run in N processes. Fits don't depend on this
setting: the same evaluations are made, only faster. `api.list_optimizers()` lists all
optimizers, including the multi-objective `nsga2` and `unsga3` used for process optimization.
