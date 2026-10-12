"""Optimize page: method settings for the best pool at a required purity, Pareto fronts (U6)."""

import pandas as pd
from nicegui import ui

from chromsim import api
from chromsim.ui import components

OBJECTIVES = ("recovery", "productivity", "buffer_consumption", "concentration")


def form(workspace: api.Workspace, project_name: str):
    project = workspace.storage.load_project(project_name)
    names = project.component_names
    with ui.row().classes("gap-2 items-end"):
        result_name = ui.input("Result name", value=f"{project_name}-opt").classes("w-60")
        target = ui.select(
            names, value=names[1] if len(names) > 1 else names[0], label="Product"
        ).classes("w-40")
        purity = ui.number("Purity [%]", value=99.0).classes("w-32")
        objectives = ui.select(
            list(OBJECTIVES), value=["recovery"], multiple=True, label="Objectives"
        ).classes("w-72")
        steps = [step.name for step in project.method]
        pool_steps = ui.select(
            steps,
            value=[s for s in steps if "elut" in s.lower()] or None,
            multiple=True,
            label="Pool during",
        ).classes("w-60")
    ui.label("Method settings to optimise").classes("font-medium")
    choices = []
    with ui.grid(columns=4).classes("items-center gap-x-4 gap-y-0"):
        for candidate in api.process_variable_candidates(project):
            choices.append(
                (
                    candidate,
                    ui.checkbox(candidate.name or candidate.path).mark(f"vary-{candidate.path}"),
                    ui.number("lower", value=candidate.lower).classes("w-28"),
                    ui.number("upper", value=candidate.upper).classes("w-28"),
                )
            )
    with ui.row().classes("gap-2 items-end"):
        max_evaluations = (
            ui.number("Max. evaluations (one objective)", value=30, format="%.0f")
            .classes("w-56")
            .mark("process-max-evaluations")
        )
        population = ui.number("Population (several objectives)", value=12, format="%.0f").classes(
            "w-56"
        )
        generations = ui.number("Generations", value=5, format="%.0f").classes("w-32")
        parallel = ui.number("Parallel processes", value=1, format="%.0f").classes("w-40")

    def build():
        several = len(objectives.value) > 1
        options = (
            {"pop_size": int(population.value), "n_gen": int(generations.value)}
            if several
            else {"max_evaluations": int(max_evaluations.value)}
        )
        return api.ProcessSetup(
            project=project_name,
            result_name=result_name.value,
            target=target.value,
            purity=purity.value / 100.0,
            variables=[
                candidate.model_copy(update={"lower": lower.value, "upper": upper.value})
                for candidate, box, lower, upper in choices
                if box.value
            ],
            objectives=list(objectives.value),
            pool_steps=list(pool_steps.value) if pool_steps.value else None,
            options=options,
            parallel=int(parallel.value),
        )

    return build


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Optimize",
        "Choose method settings (gradient length, end salt …) for the best product pool at a "
        "required purity. Cut points are chosen exactly for every candidate, so every result "
        "meets the purity.",
    )
    content = ui.column().classes("w-full")
    jobs = ui.column().classes("w-full")
    output = ui.column().classes("w-full")
    state = {}

    def show_form(_=None):
        content.clear()
        with content:
            if workspace.current_project is None:
                ui.label("Select a project first.")
                return
            state["build"] = form(workspace, workspace.current_project)
            ui.button("Start optimization", icon="rocket_launch", on_click=submit).mark(
                "start-optimization"
            )

    def show(status):
        output.clear()
        if status.state != "done":
            ui.notify(
                f"Optimization {status.state}: {status.message}", type="negative", multi_line=True
            )
            return
        record = workspace.storage.load_record("processes", status.result)
        setup = record["setup"]
        with output:
            if record["best_kpis"] is None:
                ui.label("No setting reached the required purity.").classes("text-red-700")
                return
            ui.label(
                "Best settings: "
                + ", ".join(f"{k} = {v:.4g}" for k, v in record["best_variables"].items())
            ).classes("font-medium")
            components.records_table([record["best_kpis"]])
            table = pd.DataFrame(record["table"])
            if len(setup["objectives"]) == 2:
                ui.plotly(
                    api.plot_pareto(table, x=setup["objectives"][0], y=setup["objectives"][1])
                ).classes("w-full")
            components.records_table(record["table"][:50])
            result = workspace.storage.load_simulation_result(f"{status.result}-best")
            pool = api.Pool(
                setup["target"], record["best_pool"]["start_s"], record["best_pool"]["end_s"]
            )
            ui.plotly(api.plot_chromatogram(result, x="CV", pools=[pool])).classes("w-full")
            ui.label(f"Best method saved as project '{status.result}-best'.").mark("process-done")

    def submit():
        try:
            setup = state["build"]()
            if not setup.variables:
                ui.notify("Tick at least one setting to optimise.", type="warning")
                return
            job_id = workspace.submit_process_optimization(setup)
        except (api.ChromSimError, ValueError) as error:
            components.show_error(error)
            return
        with jobs:
            components.job_panel(workspace, job_id, show)

    with toolbar:
        components.project_select(workspace, on_change=show_form)
        components.page_actions(
            workspace,
            lambda: api.codegen.code_for_process(str(workspace.storage.root), state["build"]()),
        )
    show_form()
