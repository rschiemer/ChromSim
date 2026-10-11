"""HDF5 files for simulation and optimization results, with provenance (SPEC §2.5, §6).

One file per result. Arrays are HDF5 datasets; metadata (provenance, solver statistics, names) are
attributes, stored as JSON text where they are nested. Files can be opened with any HDF5 tool
(h5py, HDFView) as well as read back with the functions here.

Layout of a simulation result::

    /                 attrs: kind="simulation", component_names, column_volume_m3,
                             provenance (JSON), solver_stats (JSON), steps (JSON)
    /time_s           (n_time,)
    /volume_m3        (n_time,)
    /outlet_mol_per_m3 (n_time, n_components)
    /molar_masses_kg_per_mol (n_components,)  NaN where unknown

Layout of an optimization result::

    /                 attrs: kind="optimization", optimizer, variable_names, success, message,
                             provenance (JSON)
    /x, /f            best point and its objectives
    /history/x, /history/f, /history/g, /history/time_s   every evaluation
    /pareto_x, /pareto_f, /residuals, /jacobian            when present
"""

import datetime
import getpass
import json
from dataclasses import asdict
from pathlib import Path

import h5py
import numpy as np

from chromsim.backends.base import SimulationResult
from chromsim.core.method import ResolvedStep
from chromsim.diagnostics import collect_versions
from chromsim.errors import ChromSimError
from chromsim.optimizers.base import Evaluation, OptimizationResult


def _check_kind(file: h5py.File, kind: str, path: Path) -> None:
    if file.attrs.get("kind") != kind:
        raise ChromSimError(
            f"{path} is not a ChromSim {kind} result (kind={file.attrs.get('kind')!r})."
        )


def write_simulation_result(result: SimulationResult, path: str | Path) -> Path:
    """Write a simulation result to HDF5 (overwrites ``path``)."""
    path = Path(path)
    with h5py.File(path, "w") as file:
        file.attrs["kind"] = "simulation"
        file.attrs["component_names"] = json.dumps(result.component_names)
        file.attrs["column_volume_m3"] = result.column_volume_m3
        file.attrs["provenance"] = json.dumps(result.provenance)
        file.attrs["solver_stats"] = json.dumps(result.solver_stats)
        file.attrs["steps"] = json.dumps([asdict(step) for step in result.steps])
        file["time_s"] = result.time_s
        file["volume_m3"] = result.volume_m3
        file["outlet_mol_per_m3"] = result.outlet_mol_per_m3
        file["molar_masses_kg_per_mol"] = np.array(
            [np.nan if m is None else m for m in result.molar_masses_kg_per_mol], dtype=float
        )
    return path


def read_simulation_result(path: str | Path) -> SimulationResult:
    """Read a simulation result written by :func:`write_simulation_result`."""
    path = Path(path)
    with h5py.File(path, "r") as file:
        _check_kind(file, "simulation", path)
        steps = [
            ResolvedStep(
                **(
                    step
                    | {
                        "inlet_start_mol_per_m3": tuple(step["inlet_start_mol_per_m3"]),
                        "inlet_end_mol_per_m3": tuple(step["inlet_end_mol_per_m3"]),
                    }
                )
            )
            for step in json.loads(file.attrs["steps"])
        ]
        molar_masses = file["molar_masses_kg_per_mol"][()]
        return SimulationResult(
            component_names=json.loads(file.attrs["component_names"]),
            time_s=file["time_s"][()],
            outlet_mol_per_m3=file["outlet_mol_per_m3"][()],
            volume_m3=file["volume_m3"][()],
            column_volume_m3=float(file.attrs["column_volume_m3"]),
            steps=steps,
            molar_masses_kg_per_mol=[None if np.isnan(m) else float(m) for m in molar_masses],
            solver_stats=json.loads(file.attrs["solver_stats"]),
            provenance=json.loads(file.attrs["provenance"]),
        )


def write_optimization_result(result: OptimizationResult, path: str | Path) -> Path:
    """Write an optimization result with its full evaluation history to HDF5."""
    path = Path(path)
    provenance = {
        "versions": collect_versions(),
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "created_by": getpass.getuser(),
    }
    with h5py.File(path, "w") as file:
        file.attrs["kind"] = "optimization"
        file.attrs["optimizer"] = result.optimizer
        file.attrs["variable_names"] = json.dumps(result.variable_names)
        file.attrs["success"] = result.success
        file.attrs["message"] = result.message
        file.attrs["provenance"] = json.dumps(provenance)
        file["x"], file["f"] = result.x, result.f
        history = file.create_group("history")
        history["x"] = np.array([e.x for e in result.history])
        history["f"] = np.array([e.f for e in result.history])
        history["g"] = np.array([e.g for e in result.history])
        history["time_s"] = np.array([e.time_s for e in result.history])
        for name in ("pareto_x", "pareto_f", "residuals", "jacobian"):
            value = getattr(result, name)
            if value is not None:
                file[name] = value
    return path


def read_optimization_result(path: str | Path) -> OptimizationResult:
    """Read an optimization result written by :func:`write_optimization_result`."""
    path = Path(path)
    with h5py.File(path, "r") as file:
        _check_kind(file, "optimization", path)
        history = file["history"]
        evaluations = [
            Evaluation(x=x, f=f, g=g, time_s=float(t))
            for x, f, g, t in zip(
                history["x"][()],
                history["f"][()],
                history["g"][()],
                history["time_s"][()],
                strict=True,
            )
        ]
        optional = {
            name: file[name][()]
            for name in ("pareto_x", "pareto_f", "residuals", "jacobian")
            if name in file
        }
        return OptimizationResult(
            optimizer=str(file.attrs["optimizer"]),
            variable_names=json.loads(file.attrs["variable_names"]),
            x=file["x"][()],
            f=file["f"][()],
            success=bool(file.attrs["success"]),
            message=str(file.attrs["message"]),
            history=evaluations,
            **optional,
        )
