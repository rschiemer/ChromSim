"""Project page: create, edit and save project files (U1)."""

from nicegui import ui

from chromsim import api
from chromsim.ui import components, forms


def isotherm_editor(isotherm: dict, component_names: list[str]) -> None:
    """Isotherm model and its parameters: one row per binding component, plus global ones."""

    @ui.refreshable
    def render():
        table = api.isotherm_parameter_table(isotherm["model"])
        per_component = [p for p in table if p["per_component"]]
        global_ones = [p for p in table if not p["per_component"]]
        parameters = isotherm.setdefault("parameters", {})
        ui.label(
            "Parameters per binding component (switch off for non-binding components)"
        ).classes("text-sm")
        for component in component_names:
            binds = component in parameters
            with ui.row().classes("items-end gap-2"):
                ui.switch(
                    component,
                    value=binds,
                    on_change=lambda e, c=component: (
                        parameters.setdefault(c, {}) if e.value else parameters.pop(c, None),
                        render.refresh(),
                    ),
                ).classes("w-40")
                if binds:
                    values = parameters[component]
                    for p in per_component:
                        ui.number(
                            f"{p['name']} [{p['unit']}]",
                            value=values.get(p["name"]),
                            on_change=lambda e, v=values, n=p["name"]: (
                                v.__setitem__(n, e.value) if e.value is not None else v.pop(n, None)
                            ),
                        ).classes("w-36").tooltip(p["description"])
        if global_ones:
            values = isotherm.get("global_parameters") or {}
            isotherm["global_parameters"] = values
            with ui.row().classes("gap-2"):
                for p in global_ones:
                    ui.number(
                        f"{p['name']} [{p['unit']}]",
                        value=values.get(p["name"]),
                        on_change=lambda e, n=p["name"]: values.__setitem__(n, e.value),
                    ).classes("w-40").tooltip(p["description"])
        else:
            isotherm.pop("global_parameters", None)

    def change_model(event):
        isotherm["model"] = event.value
        render.refresh()

    ui.select(
        api.list_isotherms(), label="Isotherm", value=isotherm["model"], on_change=change_model
    ).classes("w-72")
    render()


def editor(workspace: api.Workspace, data: dict, name: str) -> None:
    component_names = [c["name"] for c in data["components"]]
    ui.input(
        "Name", value=data["name"], on_change=lambda e: data.__setitem__("name", e.value)
    ).classes("w-96")
    with ui.expansion("Components", icon="science").classes("w-full"):
        forms.list_form(data["components"], api.Component, component_names=[])
    with ui.expansion("Column", icon="view_column").classes("w-full"):
        forms.model_form(data["column"], api.Column, component_names)
    with ui.expansion("Resin", icon="grain").classes("w-full"):
        forms.model_form(data["resin"], api.Resin, component_names)
    with ui.expansion("System", icon="settings_input_component").classes("w-full"):
        forms.model_form(data.setdefault("system", {}), api.System, component_names)
    with ui.expansion("Isotherm", icon="functions").classes("w-full"):
        isotherm_editor(data["isotherm"], component_names)
    with ui.expansion("Method", icon="list").classes("w-full"):
        forms.list_form(data["method"], api.Step, component_names)
    with ui.expansion("Solver", icon="tune").classes("w-full"):
        forms.model_form(data.setdefault("solver", {}), api.SolverSettings, component_names)

    def save():
        try:
            project = api.project_from_dict(data, source="project form")
            workspace.storage.save_project(save_as.value, project)
        except api.ChromSimError as error:
            components.show_error(error)
            return
        workspace.current_project = save_as.value
        ui.notify(f"Saved project '{save_as.value}'", type="positive")

    with ui.row().classes("items-end gap-2"):
        save_as = ui.input("Save as", value=name).classes("w-72").mark("save-as")
        ui.button("Check and save", icon="save", on_click=save).mark("save-project")


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Project",
        "Describe column, resin, system, isotherm and method. Values carry units; ChromSim checks "
        "everything when you save. Start by uploading a project file (YAML).",
    )
    content = ui.column().classes("w-full")

    def show(name):
        content.clear()
        with content:
            if name is None:
                ui.label("No project in this workspace yet: upload a project file.").classes(
                    "text-gray-500"
                )
                return
            editor(workspace, workspace.storage.load_project(name).model_dump(mode="json"), name)

    async def upload_project(event):
        text = (await event.file.read()).decode("utf-8")
        try:
            project = api.project_from_yaml(text, source=event.file.name)
        except api.ChromSimError as error:
            components.show_error(error)
            return
        name = event.file.name.rsplit(".", 1)[0]
        workspace.storage.save_project(name, project)
        workspace.current_project = name
        ui.notify(f"Loaded project '{name}'", type="positive")
        ui.navigate.reload()

    with toolbar:
        components.project_select(workspace, on_change=show)
        ui.upload(label="Upload project YAML", on_upload=upload_project, auto_upload=True).props(
            "accept=.yaml,.yml"
        )
        components.page_actions(
            workspace,
            lambda: api.codegen.code_for_project(
                str(workspace.storage.root), workspace.current_project or "my project"
            ),
        )
    show(workspace.current_project)
