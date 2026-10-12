# ChromSim

A mechanistic chromatography modeling toolbox built on [CADET](https://cadet.github.io) and
maintained by scientists. You can set up, calibrate and use models of preparative ion-exchange,
HIC, mixed-mode and affinity chromatography from notebooks, the command line or an optional
local web UI.

- **Specification:** [`SPEC.md`](SPEC.md) (source of truth)
- **Working on the code:** [`CLAUDE.md`](CLAUDE.md)
- **Changes:** [`CHANGELOG.md`](CHANGELOG.md)

## Quick start

```bash
pixi install              # Python, CADET-Core, CADET-Process and all dependencies
pixi run chromsim doctor  # check the installation
pixi run test             # run the fast tests
pixi run ui               # optional web UI in the browser
```

No pixi? Use `conda env create -f environment.yml` (installs ChromSim in editable mode).

Status: Phase 7 (web UI) done. Tutorials are in `examples/`; see SPEC §11 for the roadmap.
