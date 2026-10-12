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
# # Tutorial 05 — Workspaces, background jobs and stored results
#
# The web UI (`chromsim ui`) keeps everything in a **workspace folder**: projects (YAML), uploaded
# data (read-only), simulation and optimization results (HDF5) and small records (JSON) of fits,
# optimizations and jobs. This tutorial does in Python what the UI's Simulate and Results pages
# do, so you can open a UI workspace in a notebook, or prepare one for colleagues who use the UI.
#
# You will:
#
# 1. create a workspace and store a project and a data file,
# 2. run a simulation as a background job, with progress and a result,
# 3. look at what is stored, including provenance, and open the HDF5 file.

# %%
import tempfile
import time
from pathlib import Path

from chromsim import api

ROOT = Path.cwd() if (Path.cwd() / "examples").is_dir() else Path.cwd().parent
EXAMPLES, DATA = ROOT / "examples", ROOT / "tests" / "data"

# %% [markdown]
# Jobs run in separate worker processes. These start by importing this script again, so in a
# script (not in a notebook) every cell that does work sits under `if __name__ == "__main__":`,
# which is false in the workers.

# %%
if __name__ == "__main__":
    # A temporary folder keeps this tutorial self-contained. Use a folder of your own, e.g.
    # api.Workspace("my-workspace"), to keep the results, and `chromsim ui --workspace my-workspace`
    # to open the same folder in the web UI.
    folder = tempfile.TemporaryDirectory()
    workspace = api.Workspace(folder.name, max_workers=1)
    storage = workspace.storage

# %% [markdown]
# ## 1. Projects and data
#
# Names are plain text (letters, digits, `-`, `_`, `.`). Uploaded data files are stored
# read-only and are never overwritten: storing a different file under an existing name is an error.

# %%
if __name__ == "__main__":
    project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    project.solver.n_col, project.solver.n_par = 30, 6  # coarse grid: quick for a tutorial
    storage.save_project("hmw", project)

    run01 = DATA / "akta" / "cex_mab_hmw_run01.csv"
    storage.save_raw_file(run01.name, run01.read_bytes())
    print("projects:", storage.list_projects())
    print("data files:", storage.list_raw_files())

# %% [markdown]
# ## 2. A simulation as a background job
#
# `submit_simulation` returns at once with a job id. The job runs in a separate process; its
# status (state, progress, message) is stored in the workspace, so another Python session, or
# the UI after a browser refresh, sees it too. `workspace.cancel(job_id)` asks it to stop.

# %%
if __name__ == "__main__":
    job_id = workspace.submit_simulation("hmw", "hmw-run")
    status = workspace.job(job_id)
    while not status.finished:
        time.sleep(0.5)
        status = workspace.job(job_id)
    print(f"{job_id}: {status.state} — {status.message}")
    result_name = workspace.result(job_id)  # raises api.JobError if the job failed

# %% [markdown]
# Fits and process optimizations work the same way with `submit_fit(api.FitSetup(...))` and
# `submit_process_optimization(api.ProcessSetup(...))`. Each UI page's "Show Python code" button
# prints the setup it would submit (Tutorials 03 and 04 explain the settings).

# %% [markdown]
# ## 3. What is stored
#
# Every result records who created it and when, the package versions and the backend.

# %%
if __name__ == "__main__":
    result = storage.load_simulation_result(result_name)
    print({key: result.provenance.get(key) for key in ("backend", "created_by", "created_at")})
    print(result.provenance["versions"])
    print("HDF5 file:", storage.result_file(result_name))

    for job in workspace.jobs():
        print(job.job_id, job.state, f"{job.progress:.0%}")

    api.plot_chromatogram(result, x="CV")

# %% [markdown]
# Finally, stop the worker processes (the UI does this when it shuts down).

# %%
if __name__ == "__main__":
    workspace.close()
    folder.cleanup()
