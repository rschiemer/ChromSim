# 0010 — Fractionation implemented in ChromSim, checked against CADET-Process

- Status: accepted (approved by the maintainer, 2026-10-11)
- Date: 2026-10-11
- Deviates from SPEC §4.1/§7.4 ("wraps CADET-Process fractionation")

## Context
CADET-Process's `Fractionator` and `FractionationOptimizer` only accept CADET-Process
`SimulationResults`. ChromSim converts every backend's output into its own backend-independent
`SimulationResult` (SPEC §4.2). Wrapping CADET-Process would make fractionation CADET-only (the
reference backend's results could not be pooled) and would keep CADET objects inside results.

## Decision
- Pooling and KPIs are implemented in numpy on `SimulationResult` (`analysis/fractionation.py`,
  about 150 lines).
- Cut points are chosen **exactly** on the output grid: the maximum recovery subject to purity is
  a prefix-sum problem over all start and end pairs, so no iterative fractionation optimizer is
  needed.
- **Parity:** `tests/backends/test_fractionation_parity.py` runs one CADET simulation and
  compares pooled amount, recovery, purity, productivity and eluent consumption with
  CADET-Process's `Fractionator` for the same cut points. Agreement is about 3e-7 relative. The
  conventions differ only by documented factors: productivity per packed-bed rather than solid
  volume; buffer consumption as volume per mass rather than its inverse; purity by mass when
  molecular weights are known.
- Process optimization nests the exact cut points inside an outer optimizer over method variables
  (`analysis/process_optimization.py`), using the optimizer interface (ADR 0009).

## Consequences
- Fractionation works for every backend and is easy to read and test.
- One contiguous pool per target. Multiple pools (e.g. side fractions for re-processing) would
  need an extension.
- If CADET-Process's definitions change, the parity test shows it.
