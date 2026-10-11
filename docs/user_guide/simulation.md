# Simulating a project

```python
from chromsim import api

project = api.load_project("my_project.yaml")
result = api.simulate(project)
api.plot_chromatogram(result, x="CV")
```

See Tutorial 01 (`examples/01_simulate_lge.py`) for a complete walk-through.

## What gets simulated

ChromSim builds a CADET-Process flow sheet. Units with zero volume are left out:

```text
inlet → mixer (CSTR) → tubing, first half → column → tubing, second half → outlet
```

- **Method → events.** At the start of each step, the inlet flow rate and concentrations are set.
  A gradient (`inlet_end`) is linear in time within its step.
- **Initial state.** The column starts clean and is equilibrated with the first step's mobile-phase
  modifiers: salt for SMA starts at its first-step value, and everything else starts at zero with
  nothing bound. Tracer pulses and breakthrough experiments can start at t = 0.
- **Tubing** is a dispersed plug-flow reactor. Its length follows from `tubing_volume` and
  `tubing_diameter` (default 0.5 mm). Its `tubing_axial_dispersion` (default 1e-6 m²/s) is a
  placeholder until a system tracer pulse is calibrated (Tutorial 03).
- **Mixer** is an ideal stirred tank of `mixer_volume`.

## Column models and the parameters they need

| Column model | Needs |
|---|---|
| `LRM` | `axial_dispersion`; uses the total porosity ε_t = ε_c + (1 − ε_c) ε_p |
| `LRMP` | `axial_dispersion`, `film_diffusion` for every component, ε_c, ε_p, particle radius |
| `GRM` (default) | as LRMP, plus `pore_diffusion` for every component |

If something is missing, `api.simulate` lists all of it at once, with example values, instead of
guessing.

## The result

`SimulationResult` holds:

- `time_s`, `volume_m3` and `column_volumes`: the axes, on a 1 s grid;
- `outlet_mol_per_m3` (time × component), plus `result.outlet("mAb", "g/L")` for single
  components in any unit;
- `steps`: the method on the time axis;
- `provenance`: the project file, package versions, solver settings, time and user.

## Is the grid fine enough?

```python
check = api.check_discretisation(project)
print(check.nrmse)   # per component, project grid vs doubled grid
```

If the NRMSE is well below 1 %, the grid is fine enough. Otherwise increase `solver.n_col` (and
`n_par` for the GRM).

## Backends

`solver.backend` chooses the simulator:

| Backend | Columns | Isotherms | Use |
|---|---|---|---|
| `cadet` (default) | LRM, LRMP, GRM | all CADET-native ones | everything |
| `reference` | LRM, LRMP | Python isotherms, and Linear, Langmuir, SMA (kinetic) | custom isotherms, cross-checks, teaching |

The reference backend is a readable numpy/scipy implementation (`src/chromsim/backends/reference/`,
one file per operator). It agrees with CADET to ≤ 1 % NRMSE (ADR 0008). It needs `k_kin` for
every binding component.

## Saving results

```python
api.write_simulation_result(result, "run01.h5")      # HDF5 with provenance
result = api.read_simulation_result("run01.h5")
```
