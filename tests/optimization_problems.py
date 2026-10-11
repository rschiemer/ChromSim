"""Benchmark problems for the optimizer contract tests (module-level, so they can be pickled for
parallel evaluation)."""

import numpy as np

from chromsim import api


def rosenbrock(x):
    """f = (1 − x₀)² + 100 (x₁ − x₀²)², minimum 0 at (1, 1) (Rosenbrock, Comput. J. 3 (1960))."""
    return np.array([(1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2])


def rosenbrock_residuals(x):
    """Residual form of Rosenbrock: ½‖r‖² = f/2."""
    return np.array([1 - x[0], 10 * (x[1] - x[0] ** 2)])


def disk_constraint(x):
    """x₀² + x₁² ≤ 2 (feasible when ≤ 0); the Rosenbrock minimum (1, 1) lies on its boundary."""
    return x[0] ** 2 + x[1] ** 2 - 2.0


def zdt1(x):
    """ZDT1 (Zitzler, Deb & Thiele, Evol. Comput. 8 (2000) 173): Pareto front f₂ = 1 − √f₁."""
    f1 = x[0]
    g = 1 + 9 * np.mean(x[1:])
    return np.array([f1, g * (1 - np.sqrt(f1 / g))])


class LangmuirFit:
    """Residuals of a small Langmuir fit (k_eq, q_max) against noise-free CADET data."""

    def __init__(self, project, measured):
        self.project = project
        self.measured = measured

    def __call__(self, x):
        project = api.set_values(
            self.project,
            [
                ("isotherm.parameters.A.k_eq", x[0], None),
                ("isotherm.parameters.A.q_max", x[1], None),
            ],
        )
        result = api.simulate(project)
        return (result.outlet("A") - self.measured) / np.max(self.measured)
