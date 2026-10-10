# Project files and units

A **project** describes one chromatography run completely: components, column, resin, system,
isotherm, method and solver settings. It's stored as a YAML file that you can edit by hand,
keep in version control and load from Python, the CLI or the UI.

## Example

```yaml
--8<-- "examples/projects/cex_mab_lge.yaml"
```

Load it and save it from Python:

```python
from chromsim import api

project = api.load_project("examples/projects/cex_mab_lge.yaml")
project.method[3].duration = "30 CV"          # edit; the new value is validated
api.save_project(project, "cex_mab_lge_30cv.yaml")
```

A load–save cycle loses nothing. Quantities are written back exactly as you wrote them (`20 cm`
stays `20 cm`). Optional fields that are not set are left out. Defaults, such as the solver
settings, are written explicitly so the file documents the whole run.

## Sections

| Section | Content |
|---|---|
| `schema_version` | File format version (currently 1). Older files are migrated on load. |
| `name` | Project name. |
| `components` | List of components in order. Each has a `name` and optionally `mw_kDa` (needed for g/L concentrations), `extinction_coeff_L_per_g_cm` (UV) and `charge_hint`. **Salt is component 0 for SMA-type models**; pH, if modelled, comes directly after salt. |
| `column` | `model` (`LRM`, `LRMP` or `GRM`, default GRM), `length`, `diameter`, optional `axial_dispersion` and per-component `film_diffusion`, `pore_diffusion`. |
| `resin` | `name`, `particle_radius`, optional `ionic_capacity` (IEX), `col_porosity` ε_c and `par_porosity` ε_p (both between 0 and 1). |
| `system` | `tubing_volume`, `mixer_volume`, `uv_delay_volume`, `conductivity_delay_volume` (all default 0 mL); `tubing_diameter` (default 0.5 mm) and `tubing_axial_dispersion` (default 1e-6 m^2/s, placeholder until calibrated). |
| `isotherm` | `model` (`Linear`, `Langmuir`, `SMA`) and `parameters`: one parameter set per binding component; see [Isotherms](../theory/isotherms/index.md). |
| `method` | Ordered steps, each with `name`, `duration`, `flow`, `inlet` and, for a linear gradient, `inlet_end`. |
| `solver` | `backend` (`cadet`), `n_col`, `n_par`, `abstol`. |

## Units

Write every physical quantity as a number and a unit: `20 cm`, `150 cm/h`, `50 mM`. Internally,
ChromSim converts everything to SI, the units CADET uses.

| Quantity | Examples | Stored internally as |
|---|---|---|
| length | `20 cm`, `45 um`, `5 mm` | m |
| volume | `0.25 mL`, `1 L` | m³ |
| step duration | `10 min`, `2 mL`, `5 CV` | s |
| flow | `150 cm/h` (linear velocity), `0.5 mL/min` | m³/s |
| concentration | `50 mM`, `1 M`, `1 g/L`, `1 mg/mL` | mol/m³ (identical to mM) |
| diffusivity | `1e-7 m^2/s` | m²/s |
| film transfer coefficient | `1e-5 m/s` | m/s |

How conversions work:

- **Column volume:** 1 CV = π·(d/2)²·L of the packed bed. For the example column
  (20 cm × 0.5 cm), 1 CV = 3.93 mL.
- **Linear velocity:** this is the superficial velocity, so the volumetric flow is
  Q = u · π·(d/2)². 150 cm/h on a 0.5 cm column is 0.491 mL/min.
- **Durations** in CV or mL become times using the step's flow rate. At 150 cm/h on a 20 cm
  column, 1 CV takes L/u = 480 s, so `5 CV` lasts 2400 s.
- **Mass concentrations** (g/L) become mol/m³ using the component's `mw_kDa`
  (1 g/L of a 148 kDa mAb is 6.76·10⁻³ mol/m³). Without `mw_kDa`, use mM.

## When something is wrong

Loading reports every problem at once, with where it is and how to fix it:

```text
bad.yaml has 2 problem(s):
  - column.length: '20 mL' is not a length; use e.g. cm or mm.
  - resin.par_porosity: Input should be less than 1
```

Misspelled keys are errors too (`lenght` is not silently ignored).

!!! note "YAML and scientific notation"
    YAML reads `1.0e2` (no sign in the exponent) as text, not a number. ChromSim converts such
    values for you, but `1.0e+2` is the safer way to write them.
