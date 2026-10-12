# Getting started

## Install

ChromSim uses [pixi](https://pixi.sh) to install Python, CADET-Core and all other dependencies
in one step:

```bash
pixi install
```

If you can't use pixi, `environment.yml` (exported from the pixi lockfile) works with conda or
mamba (`conda env create -f environment.yml`); it also installs ChromSim itself in editable mode.

## Check the installation

```bash
pixi run chromsim doctor
```

This prints the versions of Python, ChromSim, CADET-Process, CADET-Python and CADET-Core, then
runs a short CADET test simulation (a tracer pulse through a non-binding column, checked for
mass balance). The same check is available from Python:

```python
from chromsim import api

report = api.doctor()
print(report.ok, report.versions)
```

## Everyday commands

| Command | What it does |
|---|---|
| `pixi run test` | fast tests (unit + physics) |
| `pixi run test-all` | all tests, including slow workflow/backend tests |
| `pixi run lint` | ruff lint and format check |
| `pixi run docs` | serve this documentation locally |
| `pixi run ui` | start the optional web UI (see "Web UI, workspaces and jobs") |
