# Fractionation and process optimization

## Pools and KPIs

A **pool** collects everything leaving the outlet between two cut points. For a pool,
`api.pool_performance(result, pool, impurities)` reports:

| KPI | Definition | Attribute |
|---|---|---|
| recovery (yield) | target in pool / target fed in the cycle | `recovery` |
| purity | target / (target + impurities); by mass if all have `mw_kDa`, else by amount | `purity` |
| productivity | target in pool / (cycle time · packed-bed volume) | `productivity_kg_per_m3_s` (= g/(L·s)) |
| buffer consumption | pumped volume of the cycle / target in pool | `buffer_consumption_m3_per_kg` (= L/g) |
| pool concentration | target in pool / pool volume | `concentration_mol_per_m3` |

Definitions follow Schmidt-Traub et al., *Preparative Chromatography*, 2nd ed. 2012, ch. 7.
Impurities default to every component except the target, salt and pH
(`api.default_impurities(project, target)`).

## Best cut points

```python
pool = api.optimal_pool(result, "mAb", ["HMW"], purity=0.99, pool_steps=("elute",))
```

This returns the contiguous pool with the **highest recovery at the required purity**, found
exactly on the simulation's output grid (no iterative search). It returns `None` if no pool
reaches the purity. Use `pool_steps` to restrict collection to the elution step: otherwise traces
of product collected during the wash can add recovery but dilute the pool.

Numerical undershoots of the solver (tiny negative concentrations) are treated as zero, so they
can't make a pool look purer than it is.

## Optimising the method

```python
result = api.optimize_process(
    project,
    [api.FitParameter(path="method.3.duration", lower=5, upper=40, unit="CV"),
     api.FitParameter(path="method.3.inlet_end.salt", lower=250, upper=1000, unit="mM")],
    target="mAb",
    purity=0.99,
    objectives=("recovery",),            # or ("recovery", "productivity") for a Pareto front
    pool_steps=("elute",),
)
result.best.project, result.best.pool, result.best.performance, result.table
```

- Decision variables are project paths, exactly as in fitting.
- For every candidate method, the cut points are chosen exactly, so **every reported point meets
  the purity**. Candidates where no pool reaches it get the worst objective values, plus a
  constraint that measures the shortfall for optimizers that handle constraints.
- Objectives: `recovery`, `productivity`, `buffer_consumption`, `concentration`.
- Default optimizers: `nelder_mead` for one objective, `nsga2` for several. Pass `optimizer=` and
  options such as `max_evaluations=` or `pop_size=`/`n_gen=`. `parallel=N` runs simulations side
  by side; in scripts, put the call under `if __name__ == "__main__":`.
- Search on a coarse grid, then **confirm the optimum on a finer grid** (Tutorial 04).

Tutorial 04 (`examples/04_optimize_process.py`) takes the calibrated mAb/HMW model from 77 % to
94 % recovery at 99 % purity, and maps the recovery–productivity trade-off.
