"""Contract every registered optimizer must meet (SPEC §6, §9 plugins tier).

Generated for ``list_optimizers()`` according to each optimizer's capabilities: Rosenbrock
(single objective), ZDT1 (multi-objective), a constrained Rosenbrock, history and callbacks,
and a small Langmuir fit (least squares and local optimizers; slow tier).
"""

import numpy as np
import pytest
from optimization_problems import (
    LangmuirFit,
    disk_constraint,
    rosenbrock,
    rosenbrock_residuals,
    zdt1,
)
from project_examples import make_binding_project

from chromsim import api
from chromsim.errors import ChromSimError
from chromsim.optimizers import (
    Constraint,
    OptimizationProblem,
    Variable,
    get_optimizer,
    list_optimizers,
)

OPTIMIZERS = list_optimizers()
# Generous settings so that every optimizer can show it works on the benchmark problems.
GLOBAL_OPTIONS = {
    "differential_evolution": {"max_generations": 200},
    "nsga2": {"n_gen": 150},
    "unsga3": {"n_gen": 150},
}


def optimizer(name):
    return get_optimizer(name, **GLOBAL_OPTIONS.get(name, {}))


def capabilities(name):
    return get_optimizer(name).capabilities


def rosenbrock_problem(name, **extra):
    variables = [Variable("x0", -2.0, 2.0, initial=-1.2), Variable("x1", -2.0, 2.0, initial=1.0)]
    if "least_squares" in capabilities(name):
        return OptimizationProblem(variables, residuals=rosenbrock_residuals, **extra)
    return OptimizationProblem(variables, objectives=rosenbrock, **extra)


@pytest.mark.parametrize("name", OPTIMIZERS)
def test_rosenbrock(name):
    result = optimizer(name).optimize(rosenbrock_problem(name))
    # Least squares and Nelder–Mead reach f < 1e-8. L-BFGS-B uses a central-difference gradient
    # with a step of 1e-3 of the range (chosen to stay above simulation noise); its truncation
    # error stops it at f ~ 1e-5 in Rosenbrock's curved valley (measured 1.0e-5). Evolutionary
    # methods only get close in a finite budget (measured f ~ 1e-4).
    tolerance = {"lbfgsb": 1e-4}.get(name, 1e-2 if "global" in capabilities(name) else 1e-6)
    assert rosenbrock(result.x)[0] < tolerance, f"{name}: f = {rosenbrock(result.x)[0]:.3g}"
    assert np.all((result.x >= -2.0) & (result.x <= 2.0))


@pytest.mark.parametrize("name", OPTIMIZERS)
def test_history_and_callback(name):
    calls = []
    result = optimizer(name).optimize(
        rosenbrock_problem(name), callback=lambda n, f: calls.append(n)
    )
    assert result.n_evaluations == len(calls) > 0
    assert calls == list(range(1, len(calls) + 1))
    assert result.variable_names == ["x0", "x1"]
    times = [evaluation.time_s for evaluation in result.history]
    assert times == sorted(times)


@pytest.mark.parametrize("name", [n for n in OPTIMIZERS if "multi_objective" in capabilities(n)])
def test_zdt1_pareto_front(name):
    variables = [Variable(f"x{i}", 0.0, 1.0) for i in range(5)]
    result = optimizer(name).optimize(
        OptimizationProblem(variables, objectives=zdt1, n_objectives=2)
    )
    f1, f2 = result.pareto_f[:, 0], result.pareto_f[:, 1]
    # Mean distance of the found front to the true front f₂ = 1 − √f₁; measured ~0.01.
    assert np.mean(np.abs(f2 - (1 - np.sqrt(f1)))) < 0.05
    assert f1.max() - f1.min() > 0.5  # the front is spread, not one point


@pytest.mark.parametrize("name", [n for n in OPTIMIZERS if "constraints" in capabilities(n)])
def test_constrained_rosenbrock(name):
    problem = rosenbrock_problem(name, constraints=[Constraint("disk", disk_constraint)])
    result = optimizer(name).optimize(problem)
    assert disk_constraint(result.x) <= 1e-6
    assert rosenbrock(result.x)[0] < 0.05


@pytest.mark.parametrize(
    "name", [n for n in OPTIMIZERS if "multi_objective" not in capabilities(n)]
)
def test_multi_objective_problem_is_refused(name):
    variables = [Variable("x", 0.0, 1.0), Variable("y", 0.0, 1.0)]
    with pytest.raises(ChromSimError, match="one objective only"):
        get_optimizer(name).optimize(
            OptimizationProblem(variables, objectives=zdt1, n_objectives=2)
        )


def test_unknown_optimizer_and_option():
    with pytest.raises(ChromSimError, match="Unknown optimizer 'cmaes'"):
        get_optimizer("cmaes")
    with pytest.raises(ChromSimError, match="no option\\(s\\) generations"):
        get_optimizer("nsga2", generations=3)


def test_parallel_evaluation_gives_the_same_result():
    serial = get_optimizer("least_squares").optimize(rosenbrock_problem("least_squares"))
    parallel = get_optimizer("least_squares").optimize(
        rosenbrock_problem("least_squares", parallel=2)
    )
    assert parallel.x == pytest.approx(serial.x, abs=1e-10)
    assert parallel.n_evaluations == serial.n_evaluations


@pytest.mark.slow
@pytest.mark.parametrize("name", [n for n in OPTIMIZERS if "global" not in capabilities(n)])
def test_small_langmuir_fit(name):
    # Breakthrough experiment: load until the column is saturated, then wash. The breakthrough
    # front and plateau make both parameters well determined (a short pulse that binds
    # completely leaves a flat cost surface towards strong binding, where gradients vanish).
    project = make_binding_project(
        [{"name": "A"}],
        {"model": "Langmuir", "parameters": {"A": {"k_eq": 1.0, "q_max": 10.0, "k_kin": 1.0}}},
        {"A": "1 mM"},
        n_col=20,
    )
    project.method[1].duration = "30 min"
    measured = api.simulate(project).outlet("A")
    variables = [
        Variable("k_eq", 0.05, 20.0, log_scale=True, initial=0.3),
        Variable("q_max", 1.0, 50.0, initial=5.0),
    ]
    fit = LangmuirFit(project, measured)
    if "least_squares" in capabilities(name):
        problem = OptimizationProblem(variables, residuals=fit)
    else:
        problem = OptimizationProblem(
            variables, objectives=lambda x: np.array([0.5 * np.sum(fit(x) ** 2)])
        )
    result = get_optimizer(name).optimize(problem)
    # Noise-free data: the optimum is the truth; 1 % allows for the solvers' stopping criteria.
    assert result.x == pytest.approx([1.0, 10.0], rel=0.01)
