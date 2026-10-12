"""Building blocks shared by all pages: navigation, page actions, job progress, selectors.

Layout only: every action calls ``chromsim.api`` (SPEC §8, CLAUDE.md rule 1).
"""

from collections.abc import Callable

from nicegui import ui

from chromsim import api

PAGES = (
    ("/", "Project"),
    ("/data", "Data"),
    ("/simulate", "Simulate"),
    ("/calibrate", "Calibrate"),
    ("/optimize", "Optimize"),
    ("/results", "Results"),
)


def header(workspace: api.Workspace) -> None:
    """Title, navigation through the workflow, and the workspace folder."""
    with ui.header().classes("items-center gap-6 bg-slate-800"):
        ui.label("ChromSim").classes("text-lg font-bold")
        for path, label in PAGES:
            ui.link(label, path).classes("text-white no-underline")
        ui.space()
        ui.label(str(workspace.storage.root)).classes("text-xs opacity-70")


def page_title(title: str, text: str):
    """Title and explanation; returns the toolbar row below them (project, page actions)."""
    ui.label(title).classes("text-2xl font-semibold")
    ui.markdown(text).classes("text-gray-600")
    return ui.row().classes("items-end gap-4 w-full")


def page_actions(workspace: api.Workspace, code: Callable[[], str]) -> None:
    """ "Download project YAML" and "Show equivalent Python code" (every page, SPEC §8)."""

    def download_yaml():
        name = workspace.current_project
        if name is None:
            ui.notify("Select a project first.", type="warning")
            return
        yaml_text = api.project_to_yaml(workspace.storage.load_project(name))
        ui.download.content(yaml_text, f"{name}.yaml")

    def show_code():
        try:
            text = code()
        except api.ChromSimError as error:
            ui.notify(str(error), type="warning")
            return
        with ui.dialog() as dialog, ui.card().classes("w-[60rem] max-w-full"):
            ui.label("Equivalent Python code").classes("text-lg")
            ui.code(text, language="python").classes("w-full")
            ui.button("Close", on_click=dialog.close)
        dialog.open()

    with ui.row().classes("gap-2"):
        ui.button("Download project YAML", icon="download", on_click=download_yaml).props(
            "outline"
        ).mark("download-yaml")
        ui.button("Show Python code", icon="code", on_click=show_code).props("outline").mark(
            "show-code"
        )


def project_select(
    workspace: api.Workspace, on_change: Callable[[str], None] | None = None, label: str = "Project"
):
    """Select one of the workspace's projects; the choice is the session's current project."""
    names = workspace.storage.list_projects()
    if workspace.current_project not in names:
        workspace.current_project = names[0] if names else None

    def changed(event):
        workspace.current_project = event.value
        if on_change is not None:
            on_change(event.value)

    return (
        ui.select(names, value=workspace.current_project, label=label, on_change=changed)
        .classes("w-72")
        .mark("project-select")
    )


def file_select(workspace: api.Workspace, label: str, suffixes: tuple[str, ...] = (), value=None):
    """Select an uploaded data file (optionally filtered by suffix)."""
    files = [
        f
        for f in workspace.storage.list_raw_files()
        if not suffixes or f.lower().endswith(suffixes)
    ]
    return ui.select([None, *files], value=value, label=label).classes("w-72")


def upload(
    workspace: api.Workspace, label: str, on_stored: Callable[[str], None] | None = None
) -> None:
    """Upload files into the workspace (stored read-only, never overwritten)."""

    async def handle(event):
        try:
            workspace.storage.save_raw_file(event.file.name, await event.file.read())
        except api.ChromSimError as error:
            ui.notify(str(error), type="negative")
            return
        ui.notify(f"Stored {event.file.name}", type="positive")
        if on_stored is not None:
            on_stored(event.file.name)

    ui.upload(label=label, on_upload=handle, auto_upload=True, multiple=True).classes("w-96")


def job_panel(
    workspace: api.Workspace, job_id: str, on_done: Callable[[api.JobStatus], None]
) -> None:
    """Progress bar, message and cancel button for a background job (polls every second)."""
    with ui.card().classes("w-full"):
        title = ui.label(f"Job {job_id}").classes("font-medium")
        progress = ui.linear_progress(value=0.0, show_value=False).classes("w-full")
        message = ui.label("queued").mark("job-message")
        cancel = ui.button("Cancel", icon="stop", on_click=lambda: workspace.cancel(job_id)).props(
            "flat"
        )

    def poll():
        status = workspace.job(job_id)
        progress.value = status.progress
        message.text = f"{status.state}: {status.message}"
        if status.finished:
            timer.deactivate()
            cancel.disable()
            title.text = f"Job {job_id} — {status.state}"
            on_done(status)

    timer = ui.timer(1.0, poll)


def show_error(error: Exception) -> None:
    """Show a ChromSim error (multi-line messages in a dialog)."""
    text = str(error)
    if "\n" not in text:
        ui.notify(text, type="negative", multi_line=True)
        return
    with ui.dialog() as dialog, ui.card():
        ui.label("Please fix:").classes("font-medium")
        ui.markdown("```\n" + text + "\n```")
        ui.button("Close", on_click=dialog.close)
    dialog.open()


def records_table(rows: list[dict], columns: list[str] | None = None):
    """A simple table of dict rows (numbers rounded for display)."""
    if not rows:
        return ui.label("(none)").classes("text-gray-500")
    columns = columns or list(rows[0])

    def fmt(value):
        if isinstance(value, float):
            return f"{value:.4g}"
        return value

    return ui.table(
        columns=[{"name": c, "label": c, "field": c, "align": "left"} for c in columns],
        rows=[{c: fmt(row.get(c)) for c in columns} for row in rows],
    ).classes("w-full")
