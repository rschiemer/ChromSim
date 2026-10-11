# Changelog

All notable changes to ChromSim are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Phase 4 calibration:
  - Peak analysis (`analysis/peaks.py`): moments with injection correction, plate number, HETP,
    asymmetry factor, sub-grid peak maximum.
  - Tracer workflow: system volume, ε_c, ε_t, ε_p from first moments; D_ax via the closed–closed
    dispersion model, flagged when system-dominated; `apply_to(project)`.
  - Yamamoto method derived for ChromSim's model (validated against CADET to < 1 %); helpers for
    gradient slope and delay-corrected elution salt; `plot_yamamoto`.
  - Inverse fit: `FitParameter` (project paths, lin/log scaling, shared or per-experiment),
    `SignalTarget`/`AssayTarget` with a residual registry (normalized, absolute, sd_weighted with
    LOQ censoring), multi-experiment TRF least squares with multi-start and progress callback.
  - Uncertainty: linearised covariance, 95 % intervals, correlation matrix, flags for high
    correlation and bounds; `plot_correlation`.
  - Synthetic-data recovery tests (tracer, Yamamoto, multi-experiment fit with noise).
  - Calibration fixtures (tracer pulses, four LGE runs) and Tutorial 03.
- Docs: "Calibrating a model", Yamamoto derivation, ADR 0007 (calibration conventions).

### Fixed
- Text assigned into nested dicts of a project (e.g. `step.inlet["mAb"] = "1 g/L"`) no longer
  crashes unit conversion.
- `compare_with_assay(..., warn=False)` to avoid repeated mismatch warnings during fits.
- Phase 3 experimental data:
  - `Experiment`, `Signal`, `Fraction`: imported runs in SI units, each signal on its own axis,
    fraction marks, injections, run log and provenance (incl. SHA-256 of the raw file).
  - Importer registry with `akta_csv` (UNICORN text: UTF-16/UTF-8, tab/comma/semicolon,
    decimal commas), `akta_xlsx` and `generic_csv` (column mapping); `read_chromatogram` picks
    the importer automatically. Shared importer contract tests.
  - `FractionAssay` + YAML mapping profiles (wide/long, CSV/XLSX, LOQ flags, readouts as sums of
    components); `read_assay`.
  - Preprocessing: baseline subtraction (constant/linear), UV → concentration (Lambert–Beer),
    `ConductivityCalibration` and conductivity → salt, resampling, volume → time.
  - `predict_uv`, `predict_conductivity` (with detector delay volumes); `fraction_averages`,
    `compare_with_assay` (mismatched fraction ids reported, never guessed); `plot_overlay`,
    `plot_fraction_comparison`.
  - Synthetic, realistically formatted fixtures and their generator (`tests/data/`); end-to-end
    test that the generating parameters reproduce all data within noise.
- Example project `cex_mab_hmw_lge.yaml` (mAb + HMW), NaCl conductivity calibration, two assay
  mapping profiles, Tutorial 02 (`examples/02_compare_with_data.py`).
- Docs: "Experimental data", how-to "Add a chromatogram importer", ADR 0006 (openpyxl).
- Dependency: openpyxl.
- Phase 2 CADET backend:
  - Isotherm descriptors `Linear`, `Langmuir`, `SMA` with a registry (`register_isotherm`,
    `get_isotherm`, `list_isotherms`), parameter checks (unknown/missing/out of bounds, salt as
    component 0, ionic capacity for SMA) at project load, and theory pages.
  - `CadetBackend`: Project → CADET-Process flow sheet (mixer, tubing, LRM/LRMP/GRM column) →
    CADET-Core; method steps as events with linear gradients; clear errors listing missing
    transport parameters.
  - `SimulationResult` with time, volume and CV axes, outlet in any concentration unit, and
    provenance (project YAML, versions, solver settings, time, user).
  - `api.simulate`, `api.check_discretisation`, `api.plot_chromatogram`.
  - `System.tubing_diameter` and `System.tubing_axial_dispersion` (optional, with defaults).
  - Physics tests: tracer first moments (LRM/LRMP/GRM), system dead-volume delay, gradient delay,
    linear-isotherm first moment, Langmuir breakthrough plateau and stoichiometric time, mass
    balances including SMA gradient elution.
- Tutorial 01 (`examples/01_simulate_lge.py`) simulates the SPEC §4.3 project.
- Docs: "Simulating a project", isotherm theory pages, ADR 0005 (CADET backend conventions).
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
