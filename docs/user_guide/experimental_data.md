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
| `generic_csv` | any delimited text file | needs `mapping=` (see below) |

List them with `api.list_importers()`.

### ÄKTA UNICORN exports

UNICORN writes two columns per curve (x and value). Above the data are a row of curve names
(`Chrom.1_UV1_280nm`, …) and a row of units (`ml`, `mAU`, …). ChromSim:

- strips the `Chrom.N_` prefix and classifies curves as `uv`, `conductivity`, `ph`,
  `concentration_b`, `pressure` or `other`;
- converts mAU → AU, mS/cm → S/m, % → fraction, ml → m³, min → s;
- turns the `Fractions` curve into fractions: each mark starts a fraction, and `Waste` marks
  end one;
- keeps `Injection` marks and the `Run Log`.

!!! warning "Validate with your own export"
    The UNICORN parser was developed against the documented layout and synthetic files in
    `tests/data/akta/`. Before relying on it, check it once with an anonymised export from your
    own system, and add that file as a test fixture.

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

## Preprocessing

| Function | Does |
|---|---|
| `subtract_baseline(signal, window)` | subtracts the mean of a protein-free window |
| `subtract_baseline(signal, window1, window2)` | subtracts the straight line through two windows (drift) |
| `uv_to_concentration(signal, ε, path_length)` | Lambert–Beer: c = A/(ε·l), result in kg/m³ (= g/L) |
| `conductivity_to_salt(signal, calibration)` | uses a measured salt ↔ conductivity calibration |
| `resample(signal, x)` | interpolates onto new x values |
| `volume_to_time(volume, steps)` | pumped volume → time, using the method's flow rates |

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
