# 0001 — Use CADET-Process as the process layer

- Status: accepted
- Date: 2026-10-10

## Context
ChromSim has to build CADET simulations of whole chromatography runs: a column with system dead
volumes, multi-step methods with gradients, fractionation, and optimization. CADET-Core
(the C++ solver) takes a large, deeply nested HDF5 configuration. Writing that configuration by hand
is error-prone, and the CADET-Core maintainers already ship
[CADET-Process](https://github.com/fau-advanced-separations/CADET-Process), a Python layer that
handles process assembly, events and sections, fractionation and optimizer interfaces.

## Decision
ChromSim converts its own declarative `Project` into CADET-Process objects and lets CADET-Process
drive CADET-Core. ChromSim does not write CADET-Core configurations directly.

- CADET-Core 5.x from conda-forge (`cadet`), version `>=5.1,<6`.
- CADET-Process `>=0.11` from PyPI (it is not packaged on conda-forge). Version 0.11 renamed
  `run()` to `simulate()`/`optimize()`; ChromSim uses only the new names.
- Every CADET-Process name ChromSim relies on is pinned by a test (CLAUDE.md rule 4). The first
  such tests are in `tests/backends/test_doctor_cadet.py`.

## Consequences
- ChromSim stays a thin layer and inherits CADET-Process features such as multi-column
  processes (v2).
- CADET-Process API changes can break ChromSim. The pinning tests make that visible straight
  away.
- Licensing: CADET-Core and CADET-Process are GPL-3.0. ChromSim is internal-only, so no
  distribution obligations arise. Revisit this if ChromSim is ever distributed.
- Moving to CADET-Core 6 needs its own ADR once 6.0 is final.
- CADET-Process 0.12 creates `./tmp` in the current working directory for its scratch files and
  redirects Python's global `tempfile.tempdir` there. ChromSim therefore always creates
  `Cadet(temp_dir=Path(...))` with its own temporary folder; it must be a `Path`, because a
  `str` crashes. `tests/backends/test_doctor_cadet.py` guards this for `doctor`, and the
  Phase 2 backend must do the same.
