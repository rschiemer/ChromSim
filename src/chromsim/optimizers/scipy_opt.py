"""scipy optimizers: least squares (TRF), Nelder–Mead, L-BFGS-B, differential evolution.

References: Branch, Coleman & Li, SIAM J. Sci. Comput. 21 (1999) 1 (TRF); Nelder & Mead, Comput.
J. 7 (1965) 308; Byrd, Lu, Nocedal & Zhu, SIAM J. Sci. Comput. 16 (1995) 1190 (L-BFGS-B);
Storn & Price, J. Global Optim. 11 (1997) 341 (differential evolution); Virtanen et al., Nat.
Methods 17 (2020) 261 (SciPy).
"""

import numpy as np
from scipy.optimize import NonlinearConstraint, differential_evolution, least_squares, minimize

from chromsim.optimizers.base import OptimizationProblem, OptimizationResult, OptimizerBase
from chromsim.optimizers.registry import register_optimizer

# Finite-difference step in the scaled space: 0.1 % of each variable's range.
DIFF_STEP = 1e-3


def finite_difference_jacobian(
    scaled, x: np.ndarray, r0: np.ndarray, step: float = DIFF_STEP
) -> np.ndarray:
    """Forward differences dr/dx, all columns evaluated as one (parallel) batch.

    Steps backwards where x + step would leave [0, 1].
    """
    steps = np.where(x + step <= 1.0, step, -step)
    points = [x + steps[j] * np.eye(len(x))[j] for j in range(len(x))]
    outputs = scaled.evaluate_batch(points)
    return np.column_stack([(output[2] - r0) / steps[j] for j, output in enumerate(outputs)])


@register_optimizer
class LeastSquares(OptimizerBase):
    """Trust-region-reflective least squares (scipy ``least_squares``), parallel Jacobian."""

    name = "least_squares"
    capabilities = {"bounds", "least_squares", "parallel"}
    default_options = {"max_evaluations": 200, "x_tol": 1e-8, "f_tol": 1e-8}

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)
        cache: dict = {}

        def residuals(x):
            cache["x"], cache["r"] = np.array(x), scaled.residuals(x)
            return cache["r"]

        def jacobian(x):
            r0 = (
                cache["r"]
                if "x" in cache and np.array_equal(cache["x"], x)
                else scaled.residuals(x)
            )
            return finite_difference_jacobian(scaled, np.asarray(x), r0)

        solution = least_squares(
            residuals,
            scaled.start_point(),
            jac=jacobian,
            bounds=(0.0, 1.0),
            method="trf",
            max_nfev=self.options["max_evaluations"],
            xtol=self.options["x_tol"],
            ftol=self.options["f_tol"],
        )
        return self.result(
            scaled,
            solution.x,
            solution.success,
            solution.message,
            residuals=solution.fun,
            jacobian=solution.jac,
        )


def _scalar(scaled):
    """Single objective in the scaled space, for scalar minimizers."""
    return lambda x: float(scaled.objectives(x)[0])


def central_difference_gradient(scaled, x: np.ndarray, step: float = DIFF_STEP) -> np.ndarray:
    """df/dx by central differences (one-sided at the bounds), evaluated as one (parallel) batch.

    Second-order accurate, so the step can stay large compared with simulation noise.
    """
    n = len(x)
    plus = np.minimum(x + step, 1.0)
    minus = np.maximum(x - step, 0.0)
    points = [np.where(np.arange(n) == j, plus[j], x) for j in range(n)]
    points += [np.where(np.arange(n) == j, minus[j], x) for j in range(n)]
    values = [float(output[0][0]) for output in scaled.evaluate_batch(points)]
    return np.array([(values[j] - values[n + j]) / (plus[j] - minus[j]) for j in range(n)])


@register_optimizer
class NelderMead(OptimizerBase):
    """Nelder–Mead simplex with bounds (scipy ``minimize``); derivative-free, local."""

    name = "nelder_mead"
    capabilities = {"bounds"}
    # initial_step: size of the starting simplex in the scaled space (fraction of each range).
    # scipy's default (5 % of the start value) can be tiny; objectives that are flat at small
    # scales (e.g. cut points on an output grid) then stop the search at once.
    default_options = {"max_evaluations": 400, "x_tol": 1e-6, "f_tol": 1e-10, "initial_step": 0.1}

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)
        start = scaled.start_point()
        simplex = [start]
        for j in range(len(start)):
            vertex = start.copy()
            step = self.options["initial_step"]
            vertex[j] = start[j] + step if start[j] + step <= 1.0 else start[j] - step
            simplex.append(vertex)
        solution = minimize(
            _scalar(scaled),
            start,
            method="Nelder-Mead",
            bounds=[(0.0, 1.0)] * scaled.n_variables,
            options={
                "maxfev": self.options["max_evaluations"],
                "xatol": self.options["x_tol"],
                "fatol": self.options["f_tol"],
                "adaptive": True,
                "initial_simplex": np.array(simplex),
            },
        )
        return self.result(scaled, solution.x, solution.success, solution.message)


@register_optimizer
class LBFGSB(OptimizerBase):
    """L-BFGS-B quasi-Newton with bounds (scipy ``minimize``), central-difference gradient."""

    name = "lbfgsb"
    capabilities = {"bounds", "parallel"}
    default_options = {"max_evaluations": 400, "g_tol": 1e-8}

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)
        solution = minimize(
            _scalar(scaled),
            scaled.start_point(),
            jac=lambda x: central_difference_gradient(scaled, np.asarray(x)),
            method="L-BFGS-B",
            bounds=[(0.0, 1.0)] * scaled.n_variables,
            options={
                "maxfun": self.options["max_evaluations"],
                "gtol": self.options["g_tol"],
            },
        )
        return self.result(scaled, solution.x, solution.success, solution.message)


@register_optimizer
class DifferentialEvolution(OptimizerBase):
    """Differential evolution (scipy); global, with constraints and parallel populations."""

    name = "differential_evolution"
    capabilities = {"bounds", "constraints", "global", "parallel"}
    default_options = {
        "max_generations": 100,
        "pop_size": 15,
        "seed": 0,
        "tol": 1e-8,
        "polish": True,
    }

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)

        def batch_map(function, points):
            # scipy hands over batches of points; evaluate them together (in parallel).
            return [float(output[0][0]) for output in scaled.evaluate_batch(list(points))]

        constraints = []
        if problem.constraints:

            def constraint_values(x):
                return scaled.evaluate_batch([x])[0][1]

            constraints = [NonlinearConstraint(constraint_values, -np.inf, 0.0)]
        solution = differential_evolution(
            _scalar(scaled),
            bounds=[(0.0, 1.0)] * scaled.n_variables,
            maxiter=self.options["max_generations"],
            popsize=self.options["pop_size"],
            seed=self.options["seed"],
            tol=self.options["tol"],
            polish=self.options["polish"] and not problem.constraints,
            constraints=constraints,
            workers=batch_map,
            updating="deferred",
            x0=scaled.start_point()
            if any(v.initial is not None for v in problem.variables)
            else None,
        )
        return self.result(scaled, solution.x, solution.success, solution.message)
