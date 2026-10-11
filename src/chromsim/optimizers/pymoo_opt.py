"""pymoo evolutionary optimizers: NSGA-II and U-NSGA-III, for one or several objectives with
constraints. Populations are evaluated as batches (in parallel when configured).

References: Deb, Pratap, Agarwal & Meyarivan, IEEE Trans. Evol. Comput. 6 (2002) 182 (NSGA-II);
Seada & Deb, IEEE Trans. Evol. Comput. 20 (2016) 358 (U-NSGA-III); Blank & Deb, IEEE Access 8
(2020) 89497 (pymoo).
"""

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2 as PymooNSGA2
from pymoo.algorithms.moo.unsga3 import UNSGA3 as PymooUNSGA3
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions

from chromsim.optimizers.base import OptimizationProblem, OptimizationResult, OptimizerBase
from chromsim.optimizers.registry import register_optimizer


class _PymooProblem(Problem):
    """The scaled problem as a vectorised pymoo problem."""

    def __init__(self, scaled):
        problem = scaled.problem
        super().__init__(
            n_var=scaled.n_variables,
            n_obj=problem.n_objectives,
            n_ieq_constr=len(problem.constraints),
            xl=np.zeros(scaled.n_variables),
            xu=np.ones(scaled.n_variables),
        )
        self.scaled = scaled

    def _evaluate(self, X, out, *args, **kwargs):
        outputs = self.scaled.evaluate_batch(list(X))
        out["F"] = np.array([output[0] for output in outputs])
        if self.n_ieq_constr:
            out["G"] = np.array([output[1] for output in outputs])


class _PymooOptimizer(OptimizerBase):
    capabilities = {"bounds", "multi_objective", "constraints", "global", "parallel"}
    default_options = {"pop_size": 40, "n_gen": 50, "seed": 0}

    def algorithm(self, problem: OptimizationProblem):
        raise NotImplementedError

    def optimize(self, problem: OptimizationProblem, callback=None) -> OptimizationResult:
        scaled = self.prepare(problem, callback)
        solution = minimize(
            _PymooProblem(scaled),
            self.algorithm(problem),
            ("n_gen", self.options["n_gen"]),
            seed=self.options["seed"],
            verbose=False,
        )
        if solution.X is None:  # no feasible point found
            best = min(scaled.history, key=lambda e: float(np.sum(np.maximum(e.g, 0.0))))
            return self.result(
                scaled, scaled.to_scaled(best.x), False, "no feasible solution found"
            )
        pareto_x_scaled = np.atleast_2d(solution.X)
        pareto_f = np.atleast_2d(solution.F)
        first = int(np.argmin(pareto_f[:, 0]))
        return self.result(
            scaled,
            pareto_x_scaled[first],
            True,
            f"{self.options['n_gen']} generations",
            pareto_x=np.array([scaled.to_physical(x) for x in pareto_x_scaled]),
            pareto_f=pareto_f,
        )


@register_optimizer
class NSGA2(_PymooOptimizer):
    """NSGA-II genetic algorithm (pymoo)."""

    name = "nsga2"

    def algorithm(self, problem):
        return PymooNSGA2(pop_size=self.options["pop_size"])


@register_optimizer
class UNSGA3(_PymooOptimizer):
    """U-NSGA-III with Das–Dennis reference directions (pymoo)."""

    name = "unsga3"
    default_options = {"pop_size": 40, "n_gen": 50, "seed": 0, "n_partitions": 12}

    def algorithm(self, problem):
        directions = get_reference_directions(
            "das-dennis", problem.n_objectives, n_partitions=self.options["n_partitions"]
        )
        return PymooUNSGA3(ref_dirs=directions, pop_size=self.options["pop_size"])
