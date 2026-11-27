# Experimental data

ChromSim reads chromatograms and offline fraction assays, converts them to SI units, and
compares them with simulations on a common **volume** axis. Raw files are only read, never
modified: every import records the file's SHA-256 checksum.

```python
from chromsim import api

experiment = api.read_chromatogram("run01.csv")             # importer chosen automatically
assay = api.read_assay("sec.xlsx", "profiles/sec_hplc.yaml")
comparison = api.compare_with_assay(result, assay, fractions=experiment.fractions)
```

Tutorial 02 (`examples/02_compare_with_data.py`) walks through a complete example.

## Chromatogram importers

| Importer | Files | Notes |
|---|---|---|
| `akta_csv` | ÄKTA UNICORN text export (`.csv`, `.asc`, `.txt`) | UTF-16 or UTF-8; tab, comma or semicolon; decimal commas |
| `akta_xlsx` | ÄKTA UNICORN Excel export (`.xlsx`) | first worksheet unless `sheet=` |
| `unicorn_res` | native UNICORN 3–5 result file (`.res`) | ÄKTA prime/purifier/explorer; binary, read directly |
| `generic_csv` | any delimited text file, with or without header | needs `mapping=` (see below) |
| `generic_xlsx` | any Excel table | needs `mapping=`; `sheet=` optional |

List them with `api.list_importers()`.

### ÄKTA UNICORN exports

UNICORN writes two columns per curve (x and value). Above the data are a row of curve names
(`Chrom.1_UV1_280nm`, …) and a row of units (`ml`, `mAU`, …). ChromSim:

- strips the `Chrom.N_` prefix (UNICORN 6/7) or the `<result name>:<n>_` prefix (UNICORN 5)
  and classifies curves as `uv`, `conductivity`, `ph`, `concentration_b`, `pressure`,
  `concentration` (mM, g/L) or `other`;
- converts mAU → AU, mS/cm → S/m, % → fraction, ml → m³, min → s;
- turns the `Fractions` curve into fractions: each mark starts a fraction, and `Waste` marks
  end one;
- keeps `Injection` marks and the `Run Log`.

!!! note "Tested on real exports"
    Besides the synthetic files in `tests/data/akta/`, the UNICORN importers have been run on
    eleven real files from ÄKTA prime, purifier and pure systems published by other groups
    (case study 01, `examples/case_studies/`). Still check once with an export from your own
    system: UNICORN's export options (curve selection, x axis relative to injection) vary.

### Native UNICORN 3–5 result files (`.res`)

UNICORN 5 and older store each run in one binary `.res` file. `unicorn_res` reads it without
UNICORN or any extra package: every curve with its own unit and scale factor (read from the
file, not assumed), fractions, injection marks and the logbook. The volume axis is the one the
system recorded (not shifted to the injection); injection marks are in
`experiment.injections_m3`. UNICORN 6 and later do not write `.res` files: export those results
as CSV or XLSX. The layout is described in the module docstring of
`chromsim/data/importers/unicorn_res.py`.

Legacy Excel 97–2003 exports (`.xls`) cannot be read (that needs another dependency); open them
in Excel and save as `.xlsx`, or export from UNICORN as CSV.

### Generic CSV

For other instruments, describe the columns once:

```yaml
x_column: Volume (mL)
x_unit: mL
delimiter: ","
signals:
  - {column: "UV 280 (mAU)", name: UV280, unit: mAU}
  - {column: "Conductivity (mS/cm)", name: Cond, unit: mS/cm}
```

```python
experiment = api.read_chromatogram("export.csv", importer="generic_csv", mapping="mapping.yaml")
```

Files without a header row (e.g. MATLAB `save -ascii` output: time and concentration in two
aligned columns) use `header: false`; columns are then numbered from 0:

```yaml
header: false
delimiter: whitespace
x_column: 0
x_unit: s
signals:
  - {column: 1, name: lysozyme, unit: mM}
```

Concentration signals (`mM`, `g/L`, …) are converted to mol/m³ or kg/m³ and get the kind
`concentration`, or `salt` if their name starts with "salt". Excel tables use the same mapping
with `importer="generic_xlsx"` (and `sheet=` if not the first sheet).

### Data recorded against time

Workflows work on the pumped volume. For a run at constant flow, convert with
`api.time_to_volume(signal, flow_rate_m3_per_s)`.

## Preprocessing

| Function | Does |
|---|---|
| `subtract_baseline(signal, window)` | subtracts the mean of a protein-free window |
| `subtract_baseline(signal, window1, window2)` | subtracts the straight line through two windows (drift) |
| `uv_to_concentration(signal, ε, path_length)` | Lambert–Beer: c = A/(ε·l), result in kg/m³ (= g/L) |
| `conductivity_to_salt(signal, calibration)` | uses a measured salt ↔ conductivity calibration |
| `resample(signal, x)` | interpolates onto new x values |
| `volume_to_time(volume, steps)` | pumped volume → time, using the method's flow rates |
| `time_to_volume(signal, flow_rate)` | time axis → pumped-volume axis, constant flow |
| `inlet_profile(project, "salt", delay)` | the method's inlet salt on the volume axis, shifted by a delay volume (outlet salt when no conductivity was recorded) |

Windows and x values are in SI units (m³ or s), e.g. `(2e-6, 15e-6)` for 2–15 mL.

A conductivity calibration is a YAML table, interpolated linearly between points:

```yaml
points:
  - [0 mM, 1.2 mS/cm]
  - [100 mM, 11.3 mS/cm]
  - [500 mM, 45.5 mS/cm]
```

## Predicted detector signals

`predict_uv(result, project)` and `predict_conductivity(result, project, calibration)` turn a
simulation into what the detectors would record. Each signal is shifted by its delay volume
(`system.uv_delay_volume`, `system.conductivity_delay_volume`). `plot_overlay(experiment,
simulated=[...])` draws measured curves solid and simulated curves dashed, with fractions shaded.

## Fraction assays and mapping profiles

A mapping profile describes a lab's file layout once. Two layouts are supported:

**Wide**: one row per fraction, one column per analyte (`examples/assay_profiles/total_protein_a280.yaml`):

```yaml
name: Total protein A280
layout: wide
skip_rows: 2            # title lines above the header
delimiter: ";"
decimal: ","
fraction_column: Fraction
start_volume_column: Start volume [ml]   # optional; otherwise volumes come from the ÄKTA marks
end_volume_column: End volume [ml]
analytes:
  - {column: "Total protein [mg/ml]", analyte: total_protein, unit: mg/mL, sd_column: "SD [mg/ml]"}
loq: {total_protein: 0.05 mg/mL}
readouts:
  total_protein: [mAb, HMW]   # an analyte may be the sum of several model components
```

**Long**: one row per fraction and analyte (`examples/assay_profiles/sec_hplc_long.yaml`), with
`analyte_column`, `value_column`, `unit_column` (or a fixed `unit`), `sd_column` and
`analyte_names` to rename peaks.

Values like `<LOQ`, `< 0.05`, `BLQ` or `n.d.` are kept as missing values with `below_loq = True`;
they are never replaced by zero. Readouts must be concentrations (mass or molar). Percentages
need a total concentration first.

## Comparing with simulations

`compare_with_assay(result, assay, fractions=experiment.fractions)` averages the simulated
outlet over each fraction,

$$
\bar c_i = \frac{1}{V_b - V_a}\int_{V_a}^{V_b} c_i(V - V_{delay})\, dV ,
$$

then sums the components of each readout and converts to the assay's unit. Fraction ids are
matched exactly. Ids that only appear in the assay (e.g. a `Load` sample) or only in the
chromatogram are listed in `comparison.assay_fractions_not_collected` and
`comparison.collected_fractions_without_assay`, and logged. They are never matched by guessing.
`plot_fraction_comparison(comparison, analyte, result, components)` draws measured bars with SD
error bars, simulated fraction averages, and the simulated outlet curve.
