# ChromSim — Specification

**A mechanistic chromatography modeling toolbox built on CADET, maintained by scientists.**

Status: v1 specification · Owner: R. Schiemer · Implementation: Claude Code, phase by phase (see §11)

---

## 1. Purpose and audience

ChromSim lets bioprocess scientists set up, calibrate and use mechanistic models of preparative
chromatography of single- and multi-component solutions (typically a target protein plus its
impurities, each component with its own parameter set) on **ion-exchange, HIC, mixed-mode and
affinity** resins without hand-writing CADET configurations. It is used in three ways, which must always stay equivalent:

1. **Notebooks / scripts** – the primary interface for model development.
2. **A CLI** – for batch runs, CI and reproducibility.
3. **An optional local web UI (NiceGUI)** – for routine use by people who do not want to write code.

**The GUI is strictly optional.** Every capability (U1–U7) must be fully usable from Python code
alone. The core package installs and runs without any GUI dependency: `pip install chromsim` gives
the code-only toolbox; `chromsim[ui]` (pixi feature `ui`) adds NiceGUI. Nothing in the GUI may be the
only way to do something.

**Maintainers are scientists with solid Python (numpy/scipy/pandas, classes, pytest) but no
professional software-engineering background.** Every design decision in this document is judged
first by: *"Can a chromatography scientist read, debug and extend this in an afternoon?"*

### Primary use cases (v1)
| # | Use case | Example |
|---|----------|---------|
| U1 | Simulate a single-column run from a method description | Load 30 g/L resin, wash 5 CV, 20 CV linear salt gradient |
| U2 | Import experimental chromatograms and overlay with simulation | ÄKTA UNICORN CSV/XLSX export |
| U2b | Import offline fraction assays and compare with simulated fraction pools | SEC/CEX/HIC-HPLC, HCP-ELISA — any tabular format |
| U3 | Determine column/system parameters from tracer pulses | Porosities, axial dispersion, system dead volume |
| U4 | Estimate isotherm parameters | Yamamoto method from LGE runs → SMA ν, k_eq; then inverse fit of breakthrough/elution data |
| U5 | Quantify parameter uncertainty and sensitivity | Confidence intervals, correlation matrix, local sensitivities |
| U6 | Optimize a process | Gradient slope, loading and cut points vs. yield/purity/productivity |
| U7 | Add a new isotherm or optimizer without touching core code | Drop-in plugin file + test |

### Non-goals for v1 (keep the architecture open for them)
- Multi-column continuous processes (SMB, MCSGP, CaptureSMB) – CADET-Process supports them; add in v2.
- 2D/radial column models, GPU solvers, model-predictive control.
- Multi-user server deployment, user management, database backends — **but v1 must be designed so
  these can be added later without touching science code** (§8.1).
- Writing new C++ code in CADET-Core.

---

## 2. Guiding principles

1. **Thin layer, not a rewrite.** Use CADET-Process for what it already does well (process model
   assembly, simulation, fractionation, optimizer interfaces). ChromSim adds: declarative project
   files, unit handling, experiment data handling, calibration workflows, a plugin system, and the UI.
2. **One source of truth for science.** All physics, fitting and analysis live in the Python
   package. The UI and the CLI only call public package functions. **No science in the UI.**
3. **Boring, explicit Python.** Plain functions, dataclasses/pydantic models, shallow inheritance
   (max. one level below an abstract base), no metaclasses, no dynamic code generation, no clever
   decorators beyond a simple registry.
4. **Units at the boundary, SI inside.** Users write `mL`, `CV`, `cm/h`, `mM`, `g/L`; everything
   inside is SI as CADET expects (note: `mol/m³ ≡ mM`). Conversion happens in exactly one module.
5. **Reproducibility by default.** Every result stores the project file, package versions
   (ChromSim, CADET-Core, CADET-Process), solver settings and a timestamp.
6. **Tests encode the physics.** Mass balance, analytical limits and backend cross-checks are
   tests, not comments.
7. **Extend by adding files, not editing them.** New isotherms, optimizers, importers and objective
   functions are registered plugins.

---

## 3. Technology stack

| Concern | Choice | Rationale |
|---|---|---|
| Language | Python ≥ 3.11 (target 3.12) | CADET-Process requires ≥ 3.10 |
| Simulation engine | **CADET-Core 5.x** (conda-forge `cadet`) | Stable line. 6.0 is in pre-release (alpha as of mid-2026); plan an upgrade ADR when it is final |
| Process layer | **CADET-Process ≥ 0.11** | Process assembly, events, fractionation, optimizer adapters. Note 0.11 renamed `run()` → `simulate()` / `optimize()` |
| Environment | **pixi** (conda-forge + PyPI, lockfile) with an exported `environment.yml` fallback | CADET-Core binaries come from conda-forge; one command setup, reproducible lockfile |
| Data models / validation | pydantic v2 | Readable validation errors for project files and UI forms |
| Units | pint | Single conversion module (`chromsim/units.py`) |
| Numerics | numpy, scipy, pandas | Familiar to maintainers |
| Plotting | plotly (UI + notebooks), matplotlib optional for publication figures | Interactive overlays |
| Optimization | scipy (least_squares, minimize, differential_evolution), pymoo; optional: optuna / Ax | Behind one interface (§6) |
| Results storage | HDF5 via h5py (one file per run/study) + YAML project file | CADET already uses HDF5 |
| UI (optional extra) | **NiceGUI** — Streamlit, Dash and Django considered and rejected (see §8) | Pure-Python, event-driven, built on FastAPI → server-capable later |
| Docs | MkDocs Material + mkdocstrings (numpy docstrings) | Markdown, easy for scientists |
| Quality | ruff (lint+format), pytest, pytest-cov, pre-commit, mypy in non-strict mode on `chromsim/` only | Low ceremony |
| CI | GitHub Actions: lint, tests (Linux + macOS + Windows), docs build | CADET binaries differ by OS |
| Notebooks | Jupytext "percent" `.py` files are the source of truth for examples; `.ipynb` generated | Diffable in git |

Dependencies beyond this list require an ADR (§10).

---

## 4. Architecture

```
            ┌───────────────┐   ┌────────────┐   ┌─────────────┐
            │    Web UI     │   │    CLI     │   │  Notebooks  │
            └──────┬────────┘   └─────┬──────┘   └──────┬──────┘
                   └──────────── public API (chromsim.api) ─────┘
                                        │
   ┌──────────────┬───────────────┬─────┴────────┬──────────────┬─────────────┐
   │ core/        │ data/         │ estimation/  │ analysis/    │ io/         │
   │ domain model │ import &      │ calibration  │ peaks,       │ project     │
   │ (pydantic)   │ preprocessing │ workflows    │ fractions,   │ files,      │
   │ + units      │               │ + uncertainty│ KPIs, sens.  │ results     │
   └──────┬───────┴───────────────┴─────┬────────┴──────────────┴─────────────┘
          │                             │
   ┌──────┴─────────┐          ┌────────┴─────────┐
   │ isotherms/     │          │ optimizers/      │      ← plugin registries
   │ registry       │          │ registry         │
   └──────┬─────────┘          └──────────────────┘
          │
   ┌──────┴──────────────────── backends/ (SimulatorBackend protocol) ────────┐
   │ cadet   → CADET-Process → CADET-Core (production)                         │
   │ reference → pure numpy/scipy LRM/LRMP solver (custom isotherms, teaching, │
   │             cross-validation)                                             │
   └───────────────────────────────────────────────────────────────────────────┘
```

### 4.1 Repository layout

```
ChromSim/
├── CLAUDE.md                  # standing instructions for Claude Code
├── SPEC.md                    # this file
├── README.md
├── pixi.toml / pyproject.toml
├── src/chromsim/
│   ├── __init__.py
│   ├── api.py                 # the ONLY module UI/CLI import from
│   ├── units.py               # pint registry, CV/cm·h⁻¹ conversions
│   ├── core/
│   │   ├── components.py      # Component (name, MW, extinction coeff., charge hint)
│   │   ├── column.py          # Column geometry, ColumnModel (LRM/LRMP/GRM), transport params
│   │   ├── resin.py           # Resin (ionic capacity, particle radius, porosities)
│   │   ├── system.py          # System dead volume: tubing (DPFR) + mixer (CSTR), UV/cond. delays
│   │   ├── method.py          # Method = ordered Steps (Equilibrate, Load, Wash, Gradient, Step, Strip)
│   │   └── project.py         # Project = components + column + system + isotherm + method + settings
│   ├── isotherms/
│   │   ├── base.py            # IsothermBase, Parameter
│   │   ├── registry.py        # register_isotherm, get_isotherm, list_isotherms
│   │   ├── cadet_native.py    # descriptors mapping to CADET-Process binding models
│   │   └── custom/            # Python-defined isotherms (reference backend; examples)
│   ├── backends/
│   │   ├── base.py            # SimulatorBackend Protocol, SimulationResult
│   │   ├── cadet_backend.py   # Project → CADETProcess Process → simulate → SimulationResult
│   │   └── reference/         # FV discretisation + scipy solve_ivp(BDF, sparse Jacobian)
│   ├── data/
│   │   ├── experiment.py      # Experiment: signals (UV, cond, pH), metadata, fractions, offline assays
│   │   ├── importers/         # registry; akta_csv.py, akta_xlsx.py, generic_csv.py (pycorn optional)
│   │   ├── assays.py          # FractionAssay: fraction table + analyte results, mapping profiles
│   │   └── preprocessing.py   # baseline, UV→conc, conductivity→salt calibration, delay alignment, resampling
│   ├── estimation/
│   │   ├── parameters.py      # FitParameter (path, bounds, scale: lin/log, initial)
│   │   ├── objectives.py      # residual functions (SSE, normalized, shape/time-shift); registry
│   │   ├── workflows/
│   │   │   ├── tracer.py      # porosities, Dax, system dead volume from pulses
│   │   │   ├── yamamoto.py    # LGE → SMA ν and k_eq
│   │   │   └── inverse_fit.py # general multi-experiment fit
│   │   └── uncertainty.py     # Jacobian-based CI, correlation matrix, profile likelihood (opt.)
│   ├── optimizers/
│   │   ├── base.py            # OptimizationProblem, OptimizerBase, OptimizationResult
│   │   ├── registry.py
│   │   ├── scipy_opt.py       # least_squares, minimize, differential_evolution
│   │   └── pymoo_opt.py       # NSGA-II/U-NSGA-III
│   ├── analysis/
│   │   ├── peaks.py           # moments, retention time, HETP/asymmetry
│   │   ├── fractionation.py   # wraps CADET-Process fractionation; cut points; yield/purity/productivity
│   │   └── sensitivity.py     # local finite-difference sensitivities
│   ├── io/
│   │   ├── project_file.py    # YAML ⇄ Project (schema-versioned)
│   │   └── results.py         # HDF5 writer/reader with provenance
│   ├── plotting.py            # plotly helpers shared by notebooks and UI
│   ├── cli.py                 # `chromsim doctor|simulate|fit|optimize|ui`
│   ├── ui/
│   │   ├── app.py
│   │   └── pages/             # project, data, simulate, calibrate, optimize, results
│   ├── jobs/                  # JobRunner protocol; LocalProcessRunner (v1); queue runner (future)
│   └── storage/               # Storage protocol; FileStorage (v1); database storage (future)
├── tests/
│   ├── unit/  physics/  backends/  workflows/  ui/
│   └── data/                  # small synthetic + anonymised example chromatograms
├── examples/                  # jupytext .py notebooks, numbered tutorials
└── docs/
    ├── index.md  user_guide/  how_to/  theory/  api/
    └── decisions/             # ADRs (0001-use-cadet-process.md, …)
```

`jobs/` and `storage/` are package-level services used by the UI, CLI and API, not part of `ui/`
(clarified 2026-10-10, Phase 0).

### 4.2 Key data objects (simplified)

```python
class Step(BaseModel):
    name: str
    duration: Quantity              # in CV, mL or min — converted by units.py
    flow: Quantity                  # cm/h or mL/min
    inlet: dict[str, Quantity]      # component → concentration at start
    inlet_end: dict[str, Quantity] | None = None   # for linear gradients

class Project(BaseModel):
    schema_version: int
    name: str
    components: list[Component]     # convention: salt is component 0 for SMA-type models
    column: Column
    resin: Resin
    system: System
    isotherm: IsothermConfig        # {"model": "SMA", "parameters": {...}}
    method: list[Step]
    solver: SolverSettings          # discretisation, tolerances, backend name
```

`SimulationResult` holds time (s and CV/mL), outlet concentrations per component, optional UV/
conductivity predictions, solver stats, and provenance. It is backend-agnostic.

### 4.3 Example project file (YAML)

```yaml
schema_version: 1
name: CEX_mAb_LGE
components:
  - {name: salt}
  - {name: mAb, mw_kDa: 148, extinction_coeff_L_per_g_cm: 1.4}
column:  {model: GRM, length: 20 cm, diameter: 0.5 cm}
resin:   {name: Generic CEX, particle_radius: 45 um, ionic_capacity: 1200 mM,
          col_porosity: 0.36, par_porosity: 0.85}
system:  {tubing_volume: 0.25 mL, mixer_volume: 0.2 mL}
isotherm:
  model: SMA
  parameters:
    mAb: {k_eq: 1.0e-2, nu: 6.0, sigma: 30.0, k_kin: 1.0e2}
method:
  - {name: equilibrate, duration: 5 CV,  flow: 150 cm/h, inlet: {salt: 50 mM}}
  - {name: load,        duration: 2 mL,  flow: 150 cm/h, inlet: {salt: 50 mM, mAb: 1 g/L}}
  - {name: wash,        duration: 5 CV,  flow: 150 cm/h, inlet: {salt: 50 mM}}
  - {name: elute,       duration: 20 CV, flow: 150 cm/h, inlet: {salt: 50 mM}, inlet_end: {salt: 500 mM}}
solver: {backend: cadet, n_col: 100, n_par: 25, abstol: 1.0e-8}
```

---

## 5. Isotherm extensibility (the core extension point)

Two kinds of isotherm, one registry:

**(a) CADET-native** — a descriptor that maps a ChromSim name and user-facing units to a
CADET-Process binding model class. Adding one is ~20–40 lines:

```python
@register_isotherm
class StericMassAction(CadetNativeIsotherm):
    """Steric mass action (Brooks & Cramer 1992). Salt must be component 0."""
    name = "SMA"
    cadet_class = "StericMassAction"
    parameters = {
        "k_eq":  Parameter(unit="dimensionless", per_component=True, bounds=(0, None), log_scale=True),
        "nu":    Parameter(unit="dimensionless", per_component=True, bounds=(0, 50)),
        "sigma": Parameter(unit="dimensionless", per_component=True, bounds=(0, 500)),
        "k_kin": Parameter(unit="1/s", per_component=True, bounds=(0, None), log_scale=True,
                           cadet_name="desorption_rate", transform="k_kin→(k_a,k_d)"),
    }
    requires_salt = True
```

v1 must ship descriptors for, by modality:

| Modality | CADET-native isotherms in v1 | Notes |
|---|---|---|
| IEX | Linear, Langmuir, SMA, Generalized Ion Exchange (GIEX; salt **and pH** dependence) | SMA first; GIEX for pH gradients |
| HIC | Mobile Phase Modulator Langmuir, CADET's HIC models (constant water activity; water on hydrophobic surfaces) | Salt type/concentration as modifier |
| Affinity | Langmuir, Bi-Langmuir (+ kinetic forms) | Elution via pH pseudo-component; MPM or GIEX-type pH dependence |
| Mixed-mode | MPM, GIEX as approximations | **First custom Python isotherms**: literature MMC models (e.g. Nfor et al. 2010 type) — candidates for later CADET upstreaming | Claude Code must read the parameter names from the installed
CADET-Process version (not from memory) and add a test that fails if CADET-Process renames them.

**(b) Custom Python** — for models CADET does not have. The scientist writes the physics as
vectorised numpy:

```python
@register_isotherm
class MyColloidalVariant(PythonIsotherm):
    name = "colloidal_variant"
    parameters = {...}
    kinetic = True
    def flux(self, c, q, p):        # c, q: arrays (n_comp, n_points); returns dq/dt
        ...
    # optional: def equilibrium(self, c, p) -> q*  (rapid-equilibrium form)
```

Python isotherms run on the **reference backend** (LRM/LRMP, finite volume, scipy BDF with sparse
Jacobian sparsity pattern). Requirements for the reference backend:
- Correct to ≤ 1 % NRMSE vs. CADET for Linear, Langmuir and SMA on shared benchmark cases (test).
- Mass balance closure ≤ 0.1 % (test).
- Readable over fast: one file per operator (convection, dispersion, film transfer, binding).
- Documented route to production: (i) CADET-Core 6 data-driven binding (spline/GPR/NN) once stable;
  (ii) upstreaming a C++ model to CADET-Core. Neither is implemented in v1 — record in an ADR.

Every isotherm (both kinds) must provide: a docstring with equations in LaTeX and a literature
reference, parameter units and bounds, a `docs/theory/isotherms/<name>.md` page, and tests
(§9). The UI builds its parameter forms automatically from `parameters`.

---

## 6. Optimizer extensibility

```python
@dataclass
class Variable:  name: str; lower: float; upper: float; log_scale: bool = False; initial: float | None = None

@dataclass
class OptimizationProblem:
    variables: list[Variable]
    objectives: Callable[[np.ndarray], np.ndarray]     # pure, picklable
    n_objectives: int = 1
    constraints: list[Constraint] = ()
    residuals: Callable | None = None                  # for least-squares optimizers
    parallel: int = 1

class OptimizerBase(ABC):
    name: str
    capabilities: set[str]   # {"bounds","least_squares","multi_objective","constraints","parallel","global"}
    @abstractmethod
    def optimize(self, problem, callback=None) -> OptimizationResult: ...
```

- The same interface serves **parameter estimation** (residuals → least squares) and **process
  optimization** (objectives/constraints, possibly multi-objective).
- Problems work in a scaled space (log/linear transform to [0, 1]) handled by the base class, not
  by each optimizer.
- `OptimizationResult` stores the full evaluation history (x, f, time) to HDF5 so runs can be
  inspected and resumed; callbacks drive UI progress bars.
- v1 adapters: `scipy.least_squares` (TRF), `scipy.minimize` (Nelder-Mead, L-BFGS-B),
  `scipy.differential_evolution`, `pymoo` NSGA-II/U-NSGA-III. Optional adapter wrapping
  CADET-Process's optimizers. Parallel evaluation via `concurrent.futures.ProcessPoolExecutor`
  (Windows/macOS-safe spawn context).
- Adding an optimizer = one file + registry entry + the shared optimizer test suite
  (Rosenbrock, a 2-objective ZDT problem, and a small Langmuir fit) passing.

---

## 7. Scientific workflows (acceptance-level detail)

### 7.1 Simulation
- Method steps are converted to CADET-Process events/sections; gradients as piecewise-linear inlets.
- System dead volume modeled as tubing (DPFR/TubularReactor) + mixer (CSTR) before and after the
  column; detector delay applied to UV/conductivity separately.
- UV prediction: Lambert-Beer with per-component extinction coefficients and path length;
  conductivity prediction from a user-supplied salt calibration curve.
- Defaults: GRM, n_col = 100, n_par = 25 (tunable); a `check_discretisation()` helper doubles the
  grid and reports the change in outlet profile.

### 7.2 Experimental data
- Importers return an `Experiment` with time-indexed signals in SI plus original units, method
  metadata where available, and fraction marks. v1: ÄKTA UNICORN CSV and XLSX exports + generic CSV
  with a column-mapping dialog. Raw files are never modified.
- Preprocessing: baseline subtraction, UV→concentration, conductivity→salt, volume↔time↔CV axes,
  alignment with the simulated system delay, resampling onto simulation time grid.

### 7.2b Offline fraction assays
- `FractionAssay` = fraction table (id, start/end volume or time) + analyte results (any number of
  analytes, units, LOQ, replicate SD). Importer works on **any tabular format** (CSV/XLSX, wide or
  long) via a saved **mapping profile** (YAML) so a lab's assay layout is configured once.
- Fraction ids are matched to the ÄKTA fraction marks; mismatches are reported, never silently fixed.
- Simulated outlet is integrated over each fraction → fraction-averaged concentrations comparable to
  assay values; plotted as bars over the chromatogram.
- Assay readouts are mapped to model components in the mapping profile (normally one-to-one;
  a readout may also be the sum of several components).

### 7.3 Calibration
1. **Tracer** (non-pore-penetrating, e.g. dextran; pore-penetrating, e.g. acetone or salt pulse;
   with and without column for system dead volume): first-moment analysis → ε_c, ε_t, ε_p;
   second moment / fit → D_ax. Report HETP and asymmetry.
2. **Yamamoto** (linear-gradient elution at several gradient slopes, low load): regression of
   log(GH) vs. log(I_R) → ν and k_eq for SMA; plot with fit and residuals.
3. **Inverse fit** on UV/conductivity traces **and/or fraction assay data** (fraction-averaged
   residuals, weighted by assay SD), across multiple experiments simultaneously (shared vs. per-experiment parameters,
   e.g. k_kin shared, load concentration per run), weights per signal, log-scaled parameters,
   multi-start option.
4. **Uncertainty**: Jacobian-based covariance → 95 % CIs and correlation matrix; flag
   |correlation| > 0.95 and parameters at bounds; optional profile likelihood for selected
   parameters.

### 7.4 Process analysis & optimization
- Fractionation via CADET-Process (cut points, pooling), KPIs: yield, purity, productivity,
  buffer consumption, concentration in pool.
- Example optimization: gradient length + end salt + cut points to maximise yield subject to
  purity ≥ X %; 2-objective yield vs. productivity Pareto front.

---

## 8. User interface

**Framework: NiceGUI**, installed only as an optional extra (`chromsim[ui]`). Record in
`docs/decisions/0003-ui-framework.md`. Rejected: Plotly Dash (callback graphs become hard to follow
in large apps); Streamlit (rerun-the-script model scales poorly to a
stateful multi-page app); Django (built for database-backed multi-user sites; would require a separate
JS frontend for interactivity and adds ORM/migrations/templates the maintainers would have to learn).
Django or FastAPI can later be layered on top of `chromsim.api` if a lab-wide platform is wanted.

Pages mirror the workflow: **Project → Data (chromatograms + assays) → Simulate → Calibrate →
Optimize → Results.**

- Forms are generated from pydantic models and isotherm `parameters`; units shown and selectable.
- Every page offers "Download project YAML" and "Show equivalent Python code" so UI users can
  move to notebooks.
- Long jobs (fits, optimizations) run via the `JobRunner` with progress callback and a cancel
  button; results are persisted immediately, so a browser refresh loses nothing.
- Plotly overlays: experiment vs. simulation, per-component traces, residuals, fraction marks,
  assay bars.
- Runs locally via `chromsim ui` (or `pixi run ui`), in the browser or as a native window
  (NiceGUI `native=True`). No authentication in v1. If NiceGUI is not installed, `chromsim ui`
  prints how to install the extra instead of failing with an ImportError.
- Use NiceGUI's built-in elements (`ui.plotly`, `ui.table`/`ui.aggrid`, tabs, dialogs); no custom
  JavaScript/Vue components without an ADR.
- Long jobs: call `JobRunner` via `run.cpu_bound`/background tasks so the event loop never blocks;
  progress via `ui.timer` polling of job status.
- **Code-only parity test:** every UI page has a matching example script in `examples/` that does
  the same thing with `chromsim.api` alone, and runs in CI.
- UI files contain only layout and calls to `chromsim.api`. A test checks that `ui/` imports nothing
  from `chromsim` except `chromsim.api`.

### 8.1 Server-readiness (design now, implement later)
- **No session state in module globals.** All per-user state lives in a `Workspace` object passed
  explicitly; the UI keeps one per browser session.
- **`JobRunner` protocol** (`submit`, `status`, `cancel`, `result`, progress events). v1:
  `LocalProcessRunner` (ProcessPoolExecutor). Later: a queue-backed runner (e.g. Celery/RQ/arq)
  with no change to callers.
- **`Storage` protocol** for projects, experiments, results. v1: `FileStorage` (YAML + HDF5 in a
  workspace folder). Later: database + object store.
- `chromsim.api` functions take and return serialisable objects (pydantic/dataclasses, file
  paths), so they can be exposed as a REST API later without redesign.
- Every stored object records `created_by` (OS user in v1) for future multi-user provenance.

## 9. Testing & quality

| Tier | Content | Runs |
|---|---|---|
| unit | units conversion, YAML round-trip, method→events, registries, importers on fixture files | every commit, < 30 s |
| physics | mass balance; linear-isotherm first moment vs. analytical; no-binding tracer = residence time; Langmuir breakthrough plateau = feed conc. | every commit |
| backends | reference vs. CADET on benchmark cases (≤ 1 % NRMSE); discretisation convergence | CI |
| workflows | recover known parameters from synthetic noisy data (tracer, Yamamoto, inverse fit) within stated tolerance | CI (marked `slow`) |
| plugins | shared contract tests auto-applied to every registered isotherm and optimizer | CI |
| regression | golden HDF5 outputs for example projects (tolerance-based) | CI |
| ui | smoke tests per page with NiceGUI's `User` test fixture | CI |
| no-gui | full test suite (except `ui/`) runs in an environment **without NiceGUI installed**; `import chromsim` must not import nicegui | CI |

Rules: tests are never weakened to make code pass; tolerances are justified in a comment;
coverage target ≥ 80 % for `core`, `isotherms`, `estimation`, `optimizers`.

---

## 10. Maintainability conventions (for humans and Claude Code)

- **Docstrings:** numpy style, with units for every physical argument and a reference for every
  model equation.
- **Naming:** physical meaning over brevity (`col_porosity`, not `eps_c`), except within equation
  implementations, where symbols must match the docstring.
- **Errors:** raise `ChromSimError` subclasses with an actionable message ("Salt must be component
  0 for SMA; found 'mAb' at index 0").
- **Logging:** standard `logging`, no prints in the package.
- **ADRs:** every significant choice gets a short Markdown record in `docs/decisions/`
  (context – decision – consequences). Initial ADRs: CADET-Process as process layer; pixi;
  reference backend; units policy; UI framework; licensing note (CADET-Core and CADET-Process are
  GPL-3.0; ChromSim is internal-only, so no distribution obligations arise — revisit if that changes).
- **Versioning:** SemVer; `CHANGELOG.md` (Keep a Changelog); project-file `schema_version` with
  migration functions.
- **Contributor docs:** `docs/how_to/add_isotherm.md`, `add_optimizer.md`, `add_importer.md` — each
  a copy-paste template walkthrough that a scientist can follow without reading core code.

---

## 11. Implementation plan for Claude Code

Build in phases. **Each phase is one or more Claude Code sessions, ends with all tests green, a
short demo notebook, and a commit.** Do not start a phase until the previous one is accepted.

| Phase | Deliverable | Acceptance criteria |
|---|---|---|
| 0 Scaffold | repo layout, pixi env, pre-commit, CI, `chromsim doctor` | `pixi run test` passes on a fresh clone; `doctor` prints Python/CADET-Core/CADET-Process versions and runs a 1-second CADET test simulation |
| 1 Domain model | `units`, `core/*`, YAML I/O | Example YAML (§4.3) round-trips losslessly; CV/mL/min and cm/h↔mL/min conversions tested |
| 2 CADET backend | Project → CADET-Process → `SimulationResult`; Linear/Langmuir/SMA; system dead volume | Tutorial 01 notebook simulates the §4.3 project; physics tests pass |
| 3 Data | `Experiment`, ÄKTA + generic importers, preprocessing, `FractionAssay` + mapping profiles, overlay plots | Fixture chromatogram + two differently formatted assay files import; fraction-averaged comparison in Tutorial 02 |
| 4 Calibration | tracer, Yamamoto, inverse fit, uncertainty | Synthetic-data recovery tests pass; Tutorial 03 |
| 5 Plugins | isotherm + optimizer registries, remaining CADET isotherms, reference backend, one custom isotherm example, optimizer adapters | Contract tests auto-run on all plugins; reference vs. CADET ≤ 1 % NRMSE; how-to guides written |
| 6 Process | fractionation, KPIs, process optimization incl. Pareto | Tutorial 04 optimizes a gradient under a purity constraint |
| 7 UI | NiceGUI pages (optional extra), `JobRunner`/`Storage` integration, code export | Smoke tests per page pass; full U1–U6 walkthrough possible in UI; no module-level state |
| 8 Docs & release | MkDocs site, ADRs, CHANGELOG, v0.1.0 tag | Docs build in CI; a colleague completes Tutorial 01 unaided |

### Decisions recorded (2026-10-10)
1. Scope: single- or multi-component solutions (target + impurities, per-component parameters) on
   IEX (primary), HIC, mixed-mode and affinity resins.
2. Data: ÄKTA/UNICORN exports plus offline fraction assays in arbitrary tabular formats.
3. Distribution: internal only.
4. UI: NiceGUI, as an optional extra. Local single-user in v1; shared-server deployment possible
   later → §8.1 applies.
5. The toolbox must be fully usable from code without the GUI (no GUI dependency in the core).
