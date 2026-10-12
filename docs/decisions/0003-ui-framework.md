# 0003 — NiceGUI for the optional web UI

- Status: accepted (SPEC §3 and §8; recorded in Phase 7)
- Date: 2026-10-12

## Context
The UI is optional (SPEC §8). It is maintained by scientists who write Python, not JavaScript.
It needs forms generated from pydantic models, Plotly plots, tables, file uploads, and progress
of long background jobs. Later it may run as a lab server for several users (SPEC §8.1).

## Decision
- **NiceGUI** (3.x), installed with the `ui` extra (`chromsim[ui]`, pixi feature `ui`).
  Streamlit (re-runs the whole script on each interaction, awkward for a stateful multi-page
  app with jobs), Dash (callback wiring gets complex for generated forms) and Django (a web
  framework with ORM, templates and a separate JS frontend) were considered and rejected
  (SPEC §8).
- One root function with `ui.sub_pages` for the six pages. `create_root(workspace_dir)` builds
  the shared `FileStorage` and `LocalProcessRunner`; every browser session gets its own
  `Workspace`. Nothing is created at import time.
- Only NiceGUI's built-in elements (`ui.plotly`, `ui.table`, tabs, dialogs, uploads); no custom
  JavaScript.
- Jobs are submitted to the `JobRunner` (separate processes) and polled with `ui.timer`, so the
  event loop never blocks. Short steps that take well under a second (reading a file, a
  Yamamoto regression) run directly.
- Tests use NiceGUI's simulated user (`nicegui.testing.user_simulation`): no browser is needed,
  and the tests run in the normal pytest suite (marker `ui`; skipped without NiceGUI).

## Consequences
- The UI is about 1 300 lines of plain Python in `src/chromsim/ui/`.
- A NiceGUI major version can change element APIs; the smoke and walkthrough tests catch this.
- Moving to a server: a queue-backed `JobRunner` and a database `Storage` plug in behind the same
  `Workspace`. Authentication would be added in front of `create_root`.
