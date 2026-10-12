"""Results page: stored simulations, fits, optimizations and jobs of the workspace."""

from nicegui import ui

from chromsim import api
from chromsim.ui import components


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Results",
        "Everything is stored in the workspace as soon as a job finishes, so nothing is lost on a "
        "browser refresh.",
    )
    view = ui.column().classes("w-full")
    state = {"result": None}

    def show(name):
        state["result"] = name
        view.clear()
        if not name:
            return
        with view:
            result = workspace.storage.load_simulation_result(name)
            provenance = result.provenance
            ui.label(
                f"{name}: backend {provenance.get('backend')}, "
                f"created {provenance.get('created_at')} by {provenance.get('created_by')}"
            )
            components.records_table([provenance.get("versions", {})])
            ui.plotly(api.plot_chromatogram(result, x="CV")).classes("w-full")
            ui.button(
                "Download HDF5",
                icon="download",
                on_click=lambda: ui.download.content(
                    workspace.storage.result_file(name).read_bytes(), f"{name}.h5"
                ),
            ).props("outline")

    with toolbar:
        ui.select(
            [None, *workspace.storage.list_simulation_results()],
            label="Simulation result",
            on_change=lambda e: show(e.value),
        ).classes("w-96").mark("result-select")
        components.page_actions(
            workspace,
            lambda: api.codegen.code_for_results(
                str(workspace.storage.root), state["result"] or "result"
            ),
        )
    with ui.expansion("Fits", icon="tune").classes("w-full"):
        components.records_table(
            [
                {"name": n} | workspace.storage.load_record("fits", n)["values"]
                for n in workspace.storage.list_records("fits")
            ]
        )
    with ui.expansion("Process optimizations", icon="rocket_launch").classes("w-full"):
        rows = []
        for name in workspace.storage.list_records("processes"):
            record = workspace.storage.load_record("processes", name)
            rows.append({"name": name} | (record["best_kpis"] or {}) | record["best_variables"])
        components.records_table(rows)
    with ui.expansion("Jobs", icon="work_history", value=True).classes("w-full"):
        components.records_table(
            [
                {
                    "job": j.job_id,
                    "state": j.state,
                    "progress": f"{j.progress:.0%}",
                    "message": j.message,
                    "created": j.created_at,
                }
                for j in workspace.jobs()
            ]
        )
