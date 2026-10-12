"""Equivalent Python code for what the web UI does (SPEC §8: "Show equivalent Python code").

Every function returns a complete, runnable script that uses only ``chromsim.api``, so a UI user
can move to notebooks. The scripts call the same functions the UI's jobs call.
"""

import pprint

_HEADER = """from chromsim import api

storage = api.FileStorage({root!r})  # the workspace folder the UI uses
"""


def _setup_literal(setup) -> str:
    return pprint.pformat(setup.model_dump(mode="json"), width=90, sort_dicts=False)


def code_for_project(root: str, project_name: str) -> str:
    """Load, check and save a project of the workspace."""
    return (
        _HEADER.format(root=root)
        + f"""
project = storage.load_project({project_name!r})
print(api.project_to_yaml(project))
# edit, e.g. project.method[3].duration = "30 CV", then save:
storage.save_project({project_name!r}, project)
"""
    )


def code_for_simulation(root: str, project_name: str, result_name: str) -> str:
    """Simulate a stored project and store and plot the result."""
    return (
        _HEADER.format(root=root)
        + f"""
project = storage.load_project({project_name!r})
result = api.simulate(project)
storage.save_simulation_result({result_name!r}, result)
api.plot_chromatogram(result, x="CV")
"""
    )


def code_for_data(root: str, chromatogram: str | None, assays: list[tuple[str, str]]) -> str:
    """Import an uploaded chromatogram and assay files."""
    lines = [_HEADER.format(root=root)]
    if chromatogram:
        lines.append(f"""experiment = api.read_chromatogram(storage.raw_file_path({chromatogram!r}))
for signal in experiment.signals.values():
    print(signal.name, signal.kind, signal.original_unit)
api.plot_overlay(experiment)
""")
    for file, profile in assays:
        lines.append(f"""assay = api.read_assay(
    storage.raw_file_path({file!r}), storage.raw_file_path({profile!r})
)
print(assay.results)
""")
    return "\n".join(lines)


def _main_guarded(body: str) -> str:
    """Put ``body`` under ``if __name__ == "__main__":`` (parallel workers re-import scripts)."""
    indented = "\n".join("    " + line if line else line for line in body.strip().splitlines())
    return (
        "\n# Parallel worker processes import this script again; the guard keeps them from\n"
        "# re-running the work (Tutorial 05).\n"
        f'if __name__ == "__main__":\n{indented}\n'
    )


def code_for_fit(root: str, setup) -> str:
    """Run an inverse fit exactly as the UI does."""
    return _HEADER.format(root=root) + _main_guarded(
        f"""
setup = api.FitSetup.model_validate(
{_setup_literal(setup)}
)
api.run_fit(storage, setup, callback=lambda n, cost: print(n, cost))
record = storage.load_record("fits", setup.result_name)
print(record["values"])
"""
    )


def code_for_process(root: str, setup) -> str:
    """Run a process optimization exactly as the UI does."""
    return _HEADER.format(root=root) + _main_guarded(
        f"""
setup = api.ProcessSetup.model_validate(
{_setup_literal(setup)}
)
api.run_process_optimization(storage, setup)
record = storage.load_record("processes", setup.result_name)
print(record["best_variables"], record["best_kpis"])
"""
    )


def code_for_results(root: str, result_name: str) -> str:
    """Open a stored simulation result."""
    return (
        _HEADER.format(root=root)
        + f"""
result = storage.load_simulation_result({result_name!r})
print(result.provenance["versions"])
api.plot_chromatogram(result, x="CV")
"""
    )
