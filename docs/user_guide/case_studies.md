# Case studies on published data

The tutorials in `examples/` use synthetic data with known answers. The case studies in
`examples/case_studies/` run ChromSim on **real data published by other groups**. The data are
not part of the repository (see [ADR 0011](../decisions/0011-third-party-data.md)); download
them first, each from a pinned commit and checked against its SHA-256:

```bash
pixi run python examples/case_studies/fetch_data.py          # all data sets
pixi run python examples/case_studies/cs01_unicorn_imports.py
pixi run python examples/case_studies/cs02_lysozyme_cex_calibration.py   # ~8 min
pixi run python examples/case_studies/cs03_yamamoto_cross_check.py
```

Respect each data set's license (listed in `fetch_data.py` and in `SOURCE.txt` next to the
data); several are for non-commercial use only.

| Case study | Data | What it shows |
|---|---|---|
| 01 UNICORN imports | 4 native UNICORN 5 `.res` files, 4 UNICORN 7 CSV exports, 3 UNICORN 5 `.asc`/CSV exports, 3 UNICORN 5 `.xls` exports (ÄKTA prime, purifier, pure) | every file imports with `read_chromatogram`; units, fractions, injections, logbook |
| 02 Lysozyme on a cation exchanger | Fraunhofer IME (Menzel & Buyel, via CADET-MCMC): tracer pulses with/without column, lysozyme gradients 6–124 CV and a step | full calibration from measured data: tracer moments → dextran-pulse fit → Yamamoto → inverse fit → prediction of runs not used in the fit |
| 03 Yamamoto cross-check | CADET-Workshop: 3 proteins × 4 gradients | ChromSim's Yamamoto result agrees with CADET-Process's to < 0.3 % (ν) and ~1 % (k_eq); how residuals reveal an inconsistent run |

## Results of case study 02

| Step | Result |
|---|---|
| Tracer moments (triplicates) | ε_c = 0.40 (dextran), ε_t = 0.925 (acetone), ε_p = 0.875 |
| Dextran-pulse fit | D_ax = 1.4·10⁻⁷ m²/s, mixer 0.17 mL, tubing 0.06 mL |
| Yamamoto (5 gradients) | ν = 3.22, k_eq = 0.058, R² = 0.999 |
| Inverse fit (10, 30, 60 CV) | ν = 3.12, k_eq = 0.073 (correlation −0.998) |
| Prediction of 6 and 124 CV gradients | peak positions within 0.3 % and 2.1 % of the elution volume |

Lessons from the real data, all discussed in the script:

- **Each tracer needs its own bypass run.** Dextran's system peak is 0.1 mL later than
  acetone's (Taylor dispersion of the large molecule in the tubing).
- **Baselines decide moments.** A 1 % offset over a long tail moves first moments by tens of µL;
  use two-window baselines and a moment window that ends at baseline.
- **Effective pore diffusion.** The measured peaks are narrower than the general rate model can
  produce with physical pore diffusion, so D_p runs to its bound: the data call for an extra
  intraparticle transport path (e.g. surface diffusion).
- **An open discrepancy.** The step elution is predicted 1.5 mL early; its peak is sharp but
  delayed, which points to a timing difference in the set-up rather than a model error.
