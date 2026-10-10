# Changelog

All notable changes to ChromSim are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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
