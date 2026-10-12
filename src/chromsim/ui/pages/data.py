"""Data page: import chromatograms and fraction assays, compare with a simulation (U2, U2b)."""

from nicegui import ui

from chromsim import api
from chromsim.ui import components

CHROMATOGRAM_SUFFIXES = (".csv", ".asc", ".txt", ".xlsx")
PROFILE_SUFFIXES = (".yaml", ".yml")


def page(workspace: api.Workspace) -> None:
    toolbar = components.page_title(
        "Data",
        "Upload ÄKTA exports, assay results, assay mapping profiles and calibrations. Files are "
        "stored read-only in the workspace and never changed.",
    )
    state = {"chromatogram": None, "assay": None, "profile": None}
    components.upload(
        workspace, "Upload data files", on_stored=lambda name: ui.navigate.to("/data")
    )

    with ui.tabs() as tabs:
        chromatogram_tab = ui.tab("Chromatogram")
        assay_tab = ui.tab("Fraction assays")
    with ui.tab_panels(tabs, value=chromatogram_tab).classes("w-full"):
        with ui.tab_panel(chromatogram_tab):
            chromatogram_view = ui.column().classes("w-full")

            def show_chromatogram(name):
                state["chromatogram"] = name
                chromatogram_view.clear()
                if not name:
                    return
                with chromatogram_view:
                    try:
                        experiment = api.read_chromatogram(workspace.storage.raw_file_path(name))
                    except api.ChromSimError as error:
                        components.show_error(error)
                        return
                    components.records_table(
                        [
                            {
                                "signal": s.name,
                                "kind": s.kind,
                                "unit": s.original_unit,
                                "points": len(s.x),
                            }
                            for s in experiment.signals.values()
                        ]
                    )
                    ui.label(
                        f"{len(experiment.fractions)} fractions, "
                        f"{len(experiment.injections_m3)} injection mark(s); "
                        f"imported with '{experiment.metadata['importer']}'"
                    )
                    ui.plotly(api.plot_overlay(experiment)).classes("w-full")

            components.file_select(
                workspace, "Chromatogram", CHROMATOGRAM_SUFFIXES
            ).on_value_change(lambda e: show_chromatogram(e.value)).mark("chromatogram-select")

        with ui.tab_panel(assay_tab):
            assay_view = ui.column().classes("w-full")
            with ui.row().classes("gap-2"):
                assay_select = components.file_select(workspace, "Assay file")
                profile_select = components.file_select(
                    workspace, "Mapping profile", PROFILE_SUFFIXES
                )
                fractions_select = components.file_select(
                    workspace, "Chromatogram (fractions)", CHROMATOGRAM_SUFFIXES
                )
                result_select = ui.select(
                    [None, *workspace.storage.list_simulation_results()],
                    label="Compare with simulation",
                ).classes("w-72")

            def show_assay():
                assay_view.clear()
                state["assay"], state["profile"] = assay_select.value, profile_select.value
                if not assay_select.value or not profile_select.value:
                    ui.notify("Choose an assay file and a mapping profile.", type="warning")
                    return
                with assay_view:
                    try:
                        storage = workspace.storage
                        assay = api.read_assay(
                            storage.raw_file_path(assay_select.value),
                            storage.raw_file_path(profile_select.value),
                        )
                        components.records_table(assay.results.to_dict(orient="records"))
                        if result_select.value:
                            fractions = None
                            if fractions_select.value:
                                fractions = api.read_chromatogram(
                                    storage.raw_file_path(fractions_select.value)
                                ).fractions
                            result = storage.load_simulation_result(result_select.value)
                            comparison = api.compare_with_assay(result, assay, fractions=fractions)
                            if comparison.assay_fractions_not_collected:
                                ui.label(
                                    "Not matched to a collected fraction: "
                                    + ", ".join(comparison.assay_fractions_not_collected)
                                )
                            for analyte in comparison.table["analyte"].unique():
                                ui.plotly(
                                    api.plot_fraction_comparison(
                                        comparison,
                                        analyte,
                                        result=result,
                                        components=assay.readouts[analyte],
                                    )
                                ).classes("w-full")
                    except api.ChromSimError as error:
                        components.show_error(error)

            ui.button("Read assay", icon="table_view", on_click=show_assay).mark("read-assay")

    with toolbar:
        components.page_actions(
            workspace,
            lambda: api.codegen.code_for_data(
                str(workspace.storage.root),
                state["chromatogram"],
                [(state["assay"], state["profile"])] if state["assay"] and state["profile"] else [],
            ),
        )
