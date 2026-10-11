"""Optimizer interface and the optimizer plugin registry (SPEC §6).

Importing this package registers the built-in optimizers.
"""

from chromsim.optimizers import pymoo_opt, scipy_opt  # noqa: F401  (registration)
from chromsim.optimizers.base import (
    Constraint,
    OptimizationProblem,
    OptimizationResult,
    OptimizerBase,
    Variable,
)
from chromsim.optimizers.registry import get_optimizer, list_optimizers, register_optimizer

__all__ = [
    "Constraint",
    "OptimizationProblem",
    "OptimizationResult",
    "OptimizerBase",
    "Variable",
    "get_optimizer",
    "list_optimizers",
    "register_optimizer",
]
