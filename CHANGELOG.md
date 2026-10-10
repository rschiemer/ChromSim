# Changelog

All notable changes to ChromSim are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Phase 1 domain model:
  - `chromsim.units`: the single unit-conversion module (pint). `Quantity` keeps the user's unit
    text for lossless files; conversions to SI for lengths, volumes, durations (time, volume or
    CV), flows (cm/h ↔ mL/min), concentrations (mM, g/L with molecular weight), diffusivities;
    column cross-section and volume.
  - `chromsim.core`: `Component`, `Column`/`ColumnModel`, `Resin`, `System`, `Step`,
    `IsothermConfig`, `SolverSettings`, `Project`, with cross-checks (unique component names,
    known components in inlets/isotherm/transport, molecular weight for g/L), and
    `resolve_method()` turning a method into an SI timeline.
  - `chromsim.io.project_file`: YAML ⇄ `Project` (`load_project`, `save_project`,
    `project_from_yaml`, `project_to_yaml`) with `schema_version` migrations and readable,
    located error messages (`ProjectFileError`).
  - New errors `UnitError` and `ProjectFileError`; all of the above exported from `chromsim.api`.
- Example project `examples/projects/cex_mab_lge.yaml` (SPEC §4.3) and notebook
  `examples/00_project_files_and_units.py`.
- Docs: user guide "Project files and units"; ADR 0004 (units policy).

### Fixed
- `chromsim doctor` no longer creates a `tmp/` folder in the current directory or changes
  Python's global temp directory (CADET-Process default); CADET now runs in its own temporary
  folder.
- Phase 0 scaffold: package layout from SPEC §4.1, pixi environment (`default` and `noui`) with
  lockfile and exported `environment.yml`, pre-commit hooks, GitHub Actions CI (lint, tests on
  Linux/macOS/Windows, no-GUI tier, docs build), MkDocs site skeleton.
- `chromsim doctor` / `chromsim.api.doctor()`: reports the Python, ChromSim, CADET-Process,
  CADET-Python and CADET-Core versions and runs a sub-second CADET tracer-pulse simulation checked
  for mass balance.
- `ChromSimError` / `EnvironmentCheckError` exception hierarchy.
- Tests pinning the CADET-Process names used by `doctor`, and a static check that `chromsim.ui`
  imports only `chromsim.api`.
- ADR 0001 (CADET-Process as process layer) and ADR 0002 (pixi).
- Example `examples/00_environment_check.py`.
