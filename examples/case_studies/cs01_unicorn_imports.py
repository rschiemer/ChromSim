# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Case study 01 — Importing real ÄKTA / UNICORN files
#
# ChromSim's importer tests use synthetic files. Here we import real files written by ÄKTA
# systems and published on GitHub by other groups (licenses in `fetch_data.py`):
#
#     pixi run python examples/case_studies/fetch_data.py unicorn5_res unicorn7_csv unicorn5_asc
#
# - `unicorn5_res`: native UNICORN 5 result files (`.res`, binary): an ÄKTA prime sample and
#   three SEC runs on a Superdex 200 16/60 (ÄKTA purifier).
# - `unicorn7_csv`: UNICORN 7 CSV exports (UTF-16, tab, `Chrom.1_` names): one SEC run and three
#   HisTrap (IMAC) + Superdex 75 runs.
# - `unicorn5_asc`: UNICORN 5 `.asc`/CSV exports (`<result>:<n>_` curve names): three runs.
#
# One call, `api.read_chromatogram(path)`, picks the right importer for each.

# %%
from pathlib import Path

import pandas as pd

from chromsim import api

HERE = Path(__file__).parent if "__file__" in globals() else Path.cwd()
DATA = HERE / "data"
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)
FOLDERS = ("unicorn5_res", "unicorn7_csv", "unicorn5_asc")
files = sorted(
    path
    for folder in FOLDERS
    for path in (DATA / folder).glob("*")
    if path.suffix.lower() in (".res", ".csv", ".asc")
)
if not files:
    raise SystemExit(f"No data in {DATA}; run examples/case_studies/fetch_data.py first.")

# %%
rows, experiments = [], {}
for path in files:
    experiment = api.read_chromatogram(path)
    experiments[path.name] = experiment
    uv = experiment.signal("uv")
    rows.append(
        {
            "file": path.name[:40],
            "importer": experiment.metadata["importer"],
            "curves": ", ".join(experiment.signals),
            "UV points": uv.x.size,
            "end / mL": round(experiment.end_m3 * 1e6, 1),
            "fractions": len(experiment.fractions),
            "injections": len(experiment.injections_m3),
            "log entries": len(experiment.log),
        }
    )
summary = pd.DataFrame(rows)
print(summary.to_string(index=False, max_colwidth=60))
assert len(experiments) == 11 and set(summary["importer"]) == {"unicorn_res", "akta_csv"}

# %% [markdown]
# What the importers take care of, file by file:
#
# - **`.res`** (binary): units and scale factors are read from each curve's header. The ÄKTA
#   prime file stores pressure in 0.001 MPa steps and the purifier files in 0.01 MPa steps,
#   so a fixed divisor would be wrong for one of them. The prime system had no pH probe: its
#   pH curve is a constant 14.64 — imported as recorded, so check before using it.
# - **UNICORN 7 CSV**: UTF-16 with tabs despite the `.csv` extension, a different sampling
#   grid for every curve, and text curves (Fraction, Injection, Run Log) that become fraction
#   marks, injection marks and log entries. The helper curve `UV_CUT_TEMP@100,BASEM` is kept as
#   a separate UV-kind curve; `experiment.signal("uv")` returns the real `UV` curve, which
#   comes first.
# - **UNICORN 5 `.asc`**: curve names carry the result name (`20220817 DAZAP1 …:10_UV`),
#   which is stripped to `UV`.

# %%
sec = experiments["KOKO038032 1660 S200.res"]
print({name: signal.original_unit for name, signal in sec.signals.items()})
print("first fractions:", [(f.label, round(f.start_m3 * 1e6, 2)) for f in sec.fractions[:4]])
print("injection at", [round(v * 1e6, 2) for v in sec.injections_m3], "mL")
figure = api.plot_overlay(sec, kinds=("uv",))
figure.write_html(OUTPUT / "cs01_sec_res.html")

# %%
imac = next(e for name, e in experiments.items() if "HisTrap" in name)
figure = api.plot_overlay(imac, kinds=("uv", "conductivity", "concentration_b"))
figure.write_html(OUTPUT / "cs01_imac_unicorn7.html")
figure
