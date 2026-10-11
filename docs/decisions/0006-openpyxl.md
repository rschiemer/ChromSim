# 0006 — openpyxl for Excel import

- Status: accepted (approved by the maintainer, 2026-10-11)
- Date: 2026-10-11

## Context
SPEC §7.2 requires importing ÄKTA UNICORN Excel exports. Assay results (SPEC §7.2b) are also
commonly delivered as Excel files. pandas reads `.xlsx` only through an engine package, and none
was installed.

## Decision
Add `openpyxl` (conda-forge, pure Python, MIT licence) as a core dependency. It is pandas'
standard `.xlsx` engine and is used directly for UNICORN sheets (read-only mode).

## Consequences
- `.xlsx` chromatograms and assays work out of the box.
- The legacy `.xls` format is not supported. Export as `.xlsx` or CSV instead.
- The fixture generator (`tests/data/generate_fixtures.py`) also uses openpyxl to write the
  synthetic Excel fixtures.
