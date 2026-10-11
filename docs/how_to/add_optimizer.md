# How to add an optimizer

An optimizer is **one new file**. The shared contract tests check it automatically
(`tests/plugins/test_optimizer_contract.py`).

## 1. Write the adapter

Create `src/chromsim/optimizers/my_optimizer.py`:

```python
"""My optimizer (reference to the algorithm's paper)."""

from chromsim.optimizers.base import OptimizationProblem, OptimizationResult, OptimizerBase
from chromsim.optimizers.registry import register_optimizer


@register_optimizer
class MyOptimizer(OptimizerBase):
    """One line: what it is and when to use it."""

    name = "my_optimizer"
    capabilities = {"bounds"}            # see below
    default_options = {"max_evaluations": 200}

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)       # checks capabilities, scales to [0, 1]
        x = scaled.start_point()
        for _ in range(self.options["max_evaluations"]):
            f = scaled.objectives(x)                   # or scaled.evaluate_batch([x1, x2, ...])
            ...                                        # your algorithm, in x ∈ [0, 1]ⁿ
        return self.result(scaled, x, success=True, message="done")
```

What the base class does for you:

- **Scaling.** Your algorithm only sees x ∈ [0, 1]ⁿ (log10 for `log_scale` variables). Results
  come back in physical units.
- **History and callbacks.** Every evaluation is recorded, and `callback(n, f)` drives progress
  bars.
- **Parallel batches.** `scaled.evaluate_batch(points)` evaluates many points, in
  `problem.parallel` processes. Use it for populations and finite-difference columns.
- **Caching.** Repeated points are not simulated twice.

Capabilities: `"bounds"`, `"least_squares"` (needs residuals; use `scaled.residuals(x)`),
`"multi_objective"` (return `pareto_x`/`pareto_f` via `self.result(..., pareto_x=..., pareto_f=...)`),
`"constraints"` (g(x) ≤ 0 from `evaluate_batch(...)[i][1]`), `"parallel"`, `"global"`.

## 2. Register it

Import the module in `src/chromsim/optimizers/__init__.py`.

## 3. Test it

Run `pixi run test-all`. Depending on the capabilities you declare, the contract solves:

- Rosenbrock (every optimizer; f < 1e-6 for local, < 1e-2 for global methods);
- ZDT1 (multi-objective): the found front lies within 0.05 of the true front;
- constrained Rosenbrock (constraints): the result is feasible;
- a small Langmuir fit (least-squares and local methods): parameters within 1 %.

It also checks history, callbacks and refusal of problems the optimizer can't handle. Add a
docs entry to the optimizer table in `docs/user_guide/calibration.md` if users should choose it
for fitting.
