"""HDF5 round trips of simulation and optimization results."""

import h5py
import numpy as np
import pytest
from optimization_problems import zdt1
from project_examples import make_physics_project

from chromsim import api
from chromsim.errors import ChromSimError
from chromsim.io.results import (
    read_optimization_result,
    read_simulation_result,
    write_optimization_result,
    write_simulation_result,
)
from chromsim.optimizers import OptimizationProblem, Variable, get_optimizer


def test_simulation_result_round_trip(tmp_path):
    result = api.simulate(make_physics_project("LRM", n_col=20))
    path = write_simulation_result(result, tmp_path / "run.h5")
    loaded = read_simulation_result(path)
    assert loaded.component_names == result.component_names
    assert loaded.time_s == pytest.approx(result.time_s)
    assert loaded.outlet_mol_per_m3 == pytest.approx(result.outlet_mol_per_m3)
    assert loaded.steps == result.steps
    assert loaded.molar_masses_kg_per_mol == [None]
    assert loaded.provenance == result.provenance
    assert api.project_from_yaml(loaded.provenance["project_yaml"]) == make_physics_project(
        "LRM", n_col=20
    )


def test_optimization_result_round_trip(tmp_path):
    variables = [Variable(f"x{i}", 0.0, 1.0) for i in range(3)]
    result = get_optimizer("nsga2", n_gen=5, pop_size=10).optimize(
        OptimizationProblem(variables, objectives=zdt1, n_objectives=2)
    )
    loaded = read_optimization_result(write_optimization_result(result, tmp_path / "opt.h5"))
    assert loaded.optimizer == "nsga2" and loaded.variable_names == ["x0", "x1", "x2"]
    assert loaded.n_evaluations == result.n_evaluations == 50
    assert loaded.history[7].x == pytest.approx(result.history[7].x)
    assert loaded.pareto_f == pytest.approx(result.pareto_f)
    assert loaded.residuals is None
    with h5py.File(tmp_path / "opt.h5") as file:  # readable with plain h5py
        assert file["history/f"].shape == (50, 2)


def test_wrong_kind_is_rejected(tmp_path):
    result = api.simulate(make_physics_project("LRM", n_col=20))
    path = write_simulation_result(result, tmp_path / "run.h5")
    with pytest.raises(ChromSimError, match="not a ChromSim optimization result"):
        read_optimization_result(path)
    assert np.isfinite(read_simulation_result(path).outlet_mol_per_m3).all()
