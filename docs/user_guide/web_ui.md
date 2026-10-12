# Web UI, workspaces and background jobs

The web UI is an optional extra (SPEC §8). Everything it does is also available from Python:
each page has a **Show Python code** button, and Tutorials 00–05 do the same as the pages.

## Start it

```bash
pixi run ui                                   # browser at http://localhost:8080
pixi run chromsim ui --workspace my-project   # another workspace folder
pixi run chromsim ui --native                 # desktop window instead of the browser
```

The UI has no login, so it only accepts connections from your own computer. `--host 0.0.0.0`
makes it reachable from the network; then everyone who can reach the port can read and change
the workspace.

Without pixi, install the extra with `pip install 'chromsim[ui]'`. If NiceGUI is missing,
`chromsim ui` says how to install it and exits.

## Pages

The pages follow the workflow. The project selected on one page stays selected on the others.

| Page | What you do | Python equivalent |
|---|---|---|
| Project | edit components, column, resin, system, isotherm (parameter table per component), method and solver; **Check and save** validates and lists every problem with its location; upload a project YAML | Tutorial 00 |
| Data | upload ÄKTA/CSV exports and assay sheets (stored read-only); view signals, fractions and injection marks; read an assay with a mapping profile and compare it with a simulation | Tutorial 02 |
| Simulate | simulate the project as a background job; plot the chromatogram | Tutorials 01, 05 |
| Calibrate | **Tracer**: porosities, system volume and dispersion from pulses; **Yamamoto**: ν and k_eq from LGE runs; **Inverse fit**: tick parameters, set bounds, fit to UV and/or assays, see the uncertainty and correlations. Results are saved as new projects | Tutorial 03 |
| Optimize | tick method settings (step lengths, gradient end concentrations), choose product, purity, objectives and the steps to pool from; see the best method, its pool and the Pareto front | Tutorial 04 |
| Results | stored simulations (plot, provenance, HDF5 download), fits, optimizations and all jobs | Tutorial 05 |

Every page also has **Download project YAML**.

## Workspaces

A workspace is a folder:

```text
my-project/
├── projects/        project YAML files
├── data/            uploaded data files (read-only, never overwritten)
├── results/         simulation results (HDF5)
├── optimizations/   optimization histories (HDF5)
└── records/         small JSON records: fits, processes, jobs
```

Each stored item records who created it and when. From Python:

```python
from chromsim import api

workspace = api.Workspace("my-project")
storage = workspace.storage
project = storage.load_project("CEX_mAb_LGE")
result = storage.load_simulation_result("run-1")
```

Data files are written once: storing a different file under an existing name is an error, so the
raw data behind a result cannot change later. The workspace interface (`api.FileStorage`) is a
protocol, so a database can replace the folder later without changes to the pages (SPEC §8.1).

## Background jobs

Simulations, fits and optimizations run in separate processes, so the UI stays responsive:

```python
job_id = workspace.submit_simulation("CEX_mAb_LGE", "run-1")
workspace.job(job_id)        # state, progress (0–1), message
workspace.cancel(job_id)     # stops at the next progress report
workspace.result(job_id)     # name of the stored result; api.JobError if not done
```

`submit_fit(api.FitSetup(...))` and `submit_process_optimization(api.ProcessSetup(...))` take
setups that name projects and files in the workspace. The UI builds the same setups, and
**Show Python code** prints them. A job's status is stored in the workspace, so a browser
refresh or a new Python session still sees it. Results are stored as soon as a job finishes.

!!! note "Scripts and worker processes"
    Worker processes start by importing your script again. In a script (not a notebook), put
    the code that submits jobs or uses `parallel=` under `if __name__ == "__main__":`.
    Generated code already does this.

## For developers

- Pages contain layout and calls to `chromsim.api` only; `tests/ui/test_ui_imports.py` enforces it.
- No module-level mutable state anywhere: `tests/unit/test_no_module_state.py`.
- Each page has a smoke test (`tests/ui/test_pages.py`, NiceGUI's simulated user). A slow test
  walks through U1–U6 in the UI (`tests/ui/test_walkthrough.py`).
- `tests/ui/test_page_parity.py` checks that each page's API calls appear in its tutorial.
- See ADR 0003 for the choice of NiceGUI.
