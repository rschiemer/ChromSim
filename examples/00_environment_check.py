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
# # 00 — Check your ChromSim installation
#
# Before modeling anything, make sure Python, CADET-Process and the CADET-Core solver are
# installed and working. `api.doctor()` collects the versions and runs a short CADET simulation:
# a tracer pulse through a non-binding column, checked for mass balance.
#
# The command line does the same thing with `pixi run chromsim doctor`.

# %%
from chromsim import api

report = api.doctor()

# %% [markdown]
# ## Installed versions

# %%
for name, version in report.versions.items():
    print(f"{name:<14} {version}")
print(f"{'Platform':<14} {report.platform}")

# %% [markdown]
# ## CADET test simulation
#
# Mass recovery is the eluted amount divided by the injected amount and should be 1.

# %%
if report.ok:
    print(f"Simulation took {report.smoke_test.wall_time_s:.2f} s")
    print(f"Mass recovery: {report.smoke_test.mass_recovery:.6f}")
else:
    for problem in report.problems:
        print("PROBLEM:", problem)

assert report.ok, "The environment check failed; see the problems listed above."
