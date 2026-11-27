"""Calibrate page: tracer analysis (U3), Yamamoto (U4), inverse fit with uncertainty (U4, U5)."""

from nicegui import ui

from chromsim import api
from chromsim.ui import components

CHROMATOGRAMS = (".csv", ".asc", ".txt", ".xlsx")
YAML = (".yaml", ".yml")


def tracer_tab(workspace: api.Workspace) -> None:
    ui.markdown(
        "Pulses of a system run (no column), a pore-excluded tracer (e.g. dextran) and a "
        "pore-penetrating tracer (e.g. acetone). Baselines and injection marks are taken from "
        "the files."
    )
    with ui.row().classes("gap-2"):
        system = components.file_select(workspace, "System run (no column)", CHROMATOGRAMS).mark(
            "tracer-system"
        )
        excluded = components.file_select(workspace, "Pore-excluded tracer", CHROMATOGRAMS).mark(
            "tracer-excluded"
        )
        penetrating = components.file_select(
            workspace, "Pore-penetrating tracer", CHROMATOGRAMS
        ).mark("tracer-penetrating")
        injection = ui.number("Injection volume [mL]", value=0.1).classes("w-40")
    output = ui.column().classes("w-full")
    state = {}

    def analyse():
        output.clear()
        if workspace.current_project is None:
            ui.notify("Select a project (column geometry and flow) first.", type="warning")
            return
        chosen = [
            (system.value, "system"),
            (excluded.value, "non_penetrating"),
            (penetrating.value, "pore_penetrating"),
        ]
        try:
            project = workspace.storage.load_project(workspace.current_project)
            runs = [
                api.tracer_run_from_file(
                    workspace.storage.raw_file_path(name), kind, injection.value
                )
                for name, kind in chosen
                if name
            ]
            flow = project.resolved_method()[0].flow_rate_m3_per_s
            state["analysis"] = api.analyse_tracers(runs, project.column, flow)
        except api.ChromSimError as error:
            components.show_error(error)
            return
        with output:
            components.records_table(api.tracer_summary(state["analysis"])).mark("tracer-summary")
            for note in state["analysis"].notes:
                ui.label("Note: " + note).classes("text-amber-700")
            with ui.row().classes("items-end gap-2"):
                save_as = ui.input(
                    "Save project as", value=f"{workspace.current_project}-tracer"
                ).classes("w-72")
                ui.button("Apply and save", on_click=lambda: apply(save_as.value)).props("outline")

    def apply(name):
        try:
            project = workspace.storage.load_project(workspace.current_project)
            workspace.storage.save_project(name, state["analysis"].apply_to(project))
        except api.ChromSimError as error:
            components.show_error(error)
            return
        ui.notify(f"Saved project '{name}'", type="positive")

    ui.button("Analyse tracers", icon="insights", on_click=analyse).mark("analyse-tracers")


def yamamoto_tab(workspace: api.Workspace) -> None:
    ui.markdown(
        "Low-load linear-gradient runs of the project's method with different gradient lengths. "
        "Gives the SMA ν and k_eq of one component as starting values for an inverse fit."
    )
    rows = []
    with ui.row().classes("gap-2"):
        calibration = components.file_select(workspace, "Conductivity calibration", YAML).mark(
            "lge-calibration"
        )
        project = (
            workspace.storage.load_project(workspace.current_project)
            if workspace.current_project
            else None
        )
        names = list(project.isotherm.parameters) if project else []
        component = ui.select(names, value=names[0] if names else None, label="Component").classes(
            "w-48"
        )
    for index in range(4):
        with ui.row().classes("gap-2"):
            rows.append(
                (
                    components.file_select(workspace, "LGE run", CHROMATOGRAMS).mark(
                        f"lge-file-{index}"
                    ),
                    ui.number("Gradient [CV]").classes("w-32").mark(f"lge-cv-{index}"),
                )
            )
    output = ui.column().classes("w-full")
    state = {}

    def fit():
        output.clear()
        try:
            project = workspace.storage.load_project(workspace.current_project)
            calibration_data = api.load_conductivity_calibration(
                workspace.storage.raw_file_path(calibration.value)
            )
            runs = [
                api.lge_run_from_file(
                    workspace.storage.raw_file_path(file.value), project, cv.value, calibration_data
                )
                for file, cv in rows
                if file.value and cv.value
            ]
            result = api.yamamoto_fit(
                runs,
                project.column.volume_m3,
                project.resin.total_porosity,
                project.resin.ionic_capacity_mol_per_m3,
            )
        except (api.ChromSimError, TypeError) as error:
            components.show_error(error)
            return
        state["result"] = result
        with output:
            ui.label(
                f"ν = {result.nu:.3f}, k_eq = {result.k_eq:.4g}, R² = {result.r_squared:.5f}"
            ).classes("font-medium").mark("yamamoto-result")
            ui.plotly(api.plot_yamamoto(result)).classes("w-full")
            with ui.row().classes("items-end gap-2"):
                save_as = ui.input(
                    "Save project as", value=f"{workspace.current_project}-yamamoto"
                ).classes("w-72")
                ui.button("Apply ν and k_eq and save", on_click=lambda: apply(save_as.value)).props(
                    "outline"
                )

    def apply(name):
        result = state["result"]
        try:
            project = workspace.storage.load_project(workspace.current_project)
            updated = api.set_values(
                project,
                [
                    (f"isotherm.parameters.{component.value}.nu", result.nu, None),
                    (f"isotherm.parameters.{component.value}.k_eq", result.k_eq, None),
                ],
            )
            workspace.storage.save_project(name, updated)
        except api.ChromSimError as error:
            components.show_error(error)
            return
        ui.notify(f"Saved project '{name}'", type="positive")

    ui.button("Yamamoto fit", icon="show_chart", on_click=fit).mark("yamamoto-fit")


def fit_tab(workspace: api.Workspace) -> dict:
    ui.markdown(
        "Fit parameters of the selected project to a chromatogram (UV) and/or fraction assays. "
        "Runs in the background; the fitted project is saved as `<result name>-fitted`."
    )
    if workspace.current_project is None:
        ui.label("Select a project first.")
        return {}
    project = workspace.storage.load_project(workspace.current_project)
    with ui.row().classes("gap-2 items-end"):
        result_name = ui.input("Result name", value=f"{workspace.current_project}-fit").classes(
            "w-60"
        )
        chromatogram = components.file_select(workspace, "Chromatogram", CHROMATOGRAMS).mark(
            "fit-chromatogram"
        )
        fit_uv = ui.switch("Fit UV", value=True)
        uv_from = ui.number("UV window from [mL]").classes("w-36").mark("fit-uv-from")
        uv_to = ui.number("to [mL]").classes("w-28").mark("fit-uv-to")
        baseline_to = ui.number("Baseline up to [mL]", value=15.0).classes("w-40")
    assays = []
    with ui.row().classes("gap-2"):
        for _ in range(2):
            assays.append(
                (
                    components.file_select(workspace, "Assay file"),
                    components.file_select(workspace, "Profile", YAML),
                )
            )
    ui.label("Parameters to fit").classes("font-medium")
    choices = []
    with ui.grid(columns=5).classes("items-center gap-x-4 gap-y-0"):
        for candidate in api.fit_parameter_candidates(project):
            choices.append(
                (
                    candidate,
                    ui.checkbox(candidate.name or candidate.path).mark(f"fit-{candidate.path}"),
                    ui.number("lower", value=candidate.lower).classes("w-28"),
                    ui.number("upper", value=candidate.upper).classes("w-28"),
                )
            )
            ui.label(candidate.unit or "").classes("text-xs")
            ui.label(candidate.scale).classes("text-xs")
    with ui.row().classes("gap-2 items-end"):
        optimizer = ui.select(
            ["least_squares", "nelder_mead", "lbfgsb", "differential_evolution"],
            value="least_squares",
            label="Optimizer",
        ).classes("w-56")
        max_evaluations = (
            ui.number("Max. evaluations", value=100, format="%.0f")
            .classes("w-40")
            .mark("fit-max-evaluations")
        )
        parallel = ui.number("Parallel processes", value=1, format="%.0f").classes("w-40")
    jobs = ui.column().classes("w-full")
    output = ui.column().classes("w-full")
    state = {}

    def build_setup():
        parameters = [
            candidate.model_copy(update={"lower": lower.value, "upper": upper.value})
            for candidate, box, lower, upper in choices
            if box.value
        ]
        windows = [(0.0, baseline_to.value)] if baseline_to.value else []
        uv_window = (
            (uv_from.value, uv_to.value)
            if uv_from.value is not None and uv_to.value is not None
            else None
        )
        return api.FitSetup(
            project=workspace.current_project,
            result_name=result_name.value,
            chromatogram=chromatogram.value,
            fit_uv=fit_uv.value,
            uv_baseline_ml=windows,
            uv_window_ml=uv_window,
            assays=[
                api.AssayInput(file=a.value, profile=p.value)
                for a, p in assays
                if a.value and p.value
            ],
            parameters=parameters,
            optimizer=optimizer.value,
            max_evaluations=int(max_evaluations.value),
            parallel=int(parallel.value),
        )

    def show(status):
        output.clear()
        if status.state != "done":
            ui.notify(f"Fit {status.state}: {status.message}", type="negative", multi_line=True)
            return
        record = workspace.storage.load_record("fits", status.result)
        with output:
            ui.label(
                f"Cost {record['cost']:.4g} after {record['n_evaluations']} evaluations — "
                f"{record['message']}"
            )
            components.records_table(record["uncertainty"])
            for a, b, r in record["high_correlations"]:
                ui.label(
                    f"{a} and {b} are strongly correlated (R = {r:.3f}): the data determine "
                    "them only together."
                ).classes("text-amber-700")
            correlation = record["correlation"]
            ui.plotly(
                api.plot_correlation_matrix(correlation["labels"], correlation["matrix"])
            ).classes("w-96")
            ui.label(f"Fitted project saved as '{status.result}-fitted'.").mark("fit-done")

    def submit():
        try:
            setup = build_setup()
            if not setup.parameters:
                ui.notify("Tick at least one parameter.", type="warning")
                return
            state["setup"] = setup
            job_id = workspace.submit_fit(setup)
        except (api.ChromSimError, ValueError) as error:
            components.show_error(error)
            return
        with jobs:
            components.job_panel(workspace, job_id, show)

    ui.button("Start fit", icon="play_arrow", on_click=submit).mark("start-fit")
    state["code"] = lambda: api.codegen.code_for_fit(str(workspace.storage.root), build_setup())
    return state


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Calibrate", "Determine column, system and isotherm parameters from experiments."
    )
    content = ui.column().classes("w-full")
    code = {}

    def show(_=None):
        content.clear()
        with content:
            with ui.tabs() as tabs:
                tracer = ui.tab("Tracer")
                yamamoto = ui.tab("Yamamoto")
                fit = ui.tab("Inverse fit")
            with ui.tab_panels(tabs, value=fit).classes("w-full"):
                with ui.tab_panel(tracer):
                    tracer_tab(workspace)
                with ui.tab_panel(yamamoto):
                    yamamoto_tab(workspace)
                with ui.tab_panel(fit):
                    code.update(fit_tab(workspace) or {})

    with toolbar:
        components.project_select(workspace, on_change=show)
        components.page_actions(
            workspace,
            lambda: (
                code["code"]()
                if "code" in code
                else api.codegen.code_for_project(
                    str(workspace.storage.root), workspace.current_project or "my project"
                )
            ),
        )
    show()
