"""Simulate page: run a project as a background job and show the chromatogram (U1)."""

from nicegui import ui

from chromsim import api
from chromsim.ui import components


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Simulate",
        "Simulate the selected project. The simulation runs in the background; its result is "
        "stored in the workspace as soon as it finishes.",
    )
    jobs = ui.column().classes("w-full")
    plot = ui.column().classes("w-full")

    def show_result(name: str) -> None:
        plot.clear()
        with plot:
            result = workspace.storage.load_simulation_result(name)
            stats = result.solver_stats
            ui.label(
                f"{name}: {result.column_volumes[-1]:.1f} CV simulated in "
                f"{stats.get('wall_time_s', 0):.1f} s"
            )
            ui.plotly(api.plot_chromatogram(result, x="CV")).classes("w-full").mark("chromatogram")

    def done(status):
        if status.state == "done":
            show_result(status.result)
        else:
            ui.notify(
                f"Simulation {status.state}: {status.message}", type="negative", multi_line=True
            )

    def run():
        if workspace.current_project is None:
            ui.notify("Select a project first.", type="warning")
            return
        job_id = workspace.submit_simulation(workspace.current_project, result_name.value)
        with jobs:
            components.job_panel(workspace, job_id, done)

    with toolbar:
        components.project_select(
            workspace, on_change=lambda name: setattr(result_name, "value", name or "")
        )
        result_name = (
            ui.input("Result name", value=workspace.current_project or "")
            .classes("w-72")
            .mark("result-name")
        )
        ui.button("Simulate", icon="play_arrow", on_click=run).mark("run-simulation")
        components.page_actions(
            workspace,
            lambda: api.codegen.code_for_simulation(
                str(workspace.storage.root),
                workspace.current_project or "my project",
                result_name.value or "result",
            ),
        )
    if result_name.value in workspace.storage.list_simulation_results():
        show_result(result_name.value)
