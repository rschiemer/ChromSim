# CLAUDE.md — standing instructions for working on ChromSim

Read `SPEC.md` before any task. It is the source of truth; if a request conflicts with it, say so
and ask before deviating. Significant deviations get an ADR in `docs/decisions/`.

## Who maintains this code
Scientists with good Python, not software engineers. Optimise for readability:
- Plain functions and pydantic/dataclass models; at most one level of inheritance below an ABC.
- No metaclasses, monkey-patching, dynamic imports beyond the plugin registries, or "clever" one-liners.
- Equation code uses the symbols from its docstring; everything else uses descriptive names.
- numpy-style docstrings with **units for every physical quantity** and a literature reference for every model equation.

## Hard rules
1. **No science and no state in the UI.** `src/chromsim/ui/` imports nothing from chromsim except
   `chromsim.api`; per-user state lives in a `Workspace`, long jobs go through `JobRunner`, persistence
   through `Storage` (SPEC §8.1). No module-level mutable state anywhere.
2. **SI inside, units at the boundary.** All conversions go through `chromsim/units.py`.
3. **Never weaken or delete a test to make code pass.** If a test seems wrong, stop and explain.
4. **Never hard-code CADET parameter names from memory** — check the installed CADET-Process
   (`inspect`, docs) and add a test pinning the mapping.
5. **No new dependencies** without asking; if approved, add to `pixi.toml` and write an ADR.
6. New isotherm / optimizer / importer = new file + registry entry + contract tests + docs page.
   Do not modify core modules to add a plugin.
7. Raise `ChromSimError` subclasses with actionable messages; use `logging`, never `print`, in the package.

## Commands
```bash
pixi install              # set up environment
pixi run test             # fast tests (unit + physics)
pixi run test-all         # includes slow workflow/backend tests
pixi run lint             # ruff check + ruff format --check
pixi run docs             # mkdocs serve
pixi run ui               # web UI
pixi run chromsim doctor  # environment + CADET sanity check
```

## Workflow per task
1. Restate the task and the SPEC phase it belongs to; list the files you will touch.
2. Write or update tests first for physics/workflow changes.
3. Implement; run `pixi run lint` and `pixi run test`.
4. Update docs (how-to/theory page, example notebook if behaviour changed) and `CHANGELOG.md`.
5. Summarise: what changed, how it was verified, any open questions. Commit with a descriptive message.

## Domain conventions
- Salt is component 0 for SMA-type models; pH, when modelled, is a non-binding pseudo-component
  placed directly after salt.
- Every binding component (target and each impurity) has its own isotherm parameter set.
- Raw data files (ÄKTA exports, assay sheets) are read-only; never modify them.
- `mol/m³ == mM`; concentrations in project files may be given in mM, g/L or mg/mL.
- Column volume (CV) = π·(d/2)²·L of the packed bed; linear velocity in cm/h is converted to
  volumetric flow using the column cross-section.
- Default column model GRM; reference backend supports LRM/LRMP only.
