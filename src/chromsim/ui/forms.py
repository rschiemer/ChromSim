"""Forms generated from pydantic models (SPEC §8): one editor for any ChromSim model.

The form edits a plain dict (``model.model_dump(mode="json")``); quantities are shown as a number
plus a unit selector, with the units ChromSim offers for that field's dimension. Validation stays
in the package: the page passes the dict to ``api.project_from_dict``.
"""

import copy
import enum
import types
import typing

from nicegui import ui
from pydantic import BaseModel

from chromsim import api


def _unwrap_optional(annotation):
    """X | None → X (keeps other unions)."""
    if typing.get_origin(annotation) in (typing.Union, types.UnionType):
        arguments = [a for a in typing.get_args(annotation) if a is not type(None)]
        if len(arguments) == 1:
            return arguments[0]
    return annotation


def _split_quantity(text) -> tuple[float | None, str | None]:
    if text in (None, ""):
        return None, None
    quantity = api.Quantity.parse(text)
    return quantity.magnitude, quantity.unit


def quantity_input(
    label: str, data: dict, key, dimensions: tuple[str, ...], on_empty=None, on_write=None
) -> None:
    """Number + unit selector writing ``"<number> <unit>"`` into ``data[key]``.

    An empty number writes None, or calls ``on_empty()`` instead if given (e.g. to remove a
    component from a per-component table).
    """
    magnitude, unit = _split_quantity(data.get(key))
    choices = api.unit_choices(dimensions)
    if unit is not None and unit not in choices:
        choices = [unit, *choices]
    unit = unit or choices[0]

    def write():
        if number.value is None:
            if on_empty is not None:
                on_empty()
            else:
                data[key] = None
        else:
            data[key] = str(api.Quantity(magnitude=float(number.value), unit=selector.value))
        if on_write is not None:
            on_write()

    with ui.row().classes("items-end gap-1"):
        number = ui.number(label, value=magnitude, on_change=write).classes("w-40")
        selector = ui.select(choices, value=unit, on_change=write).classes("w-28")


def _quantity_dict(
    label, data: dict, key: str, dimensions, component_names, optional: bool
) -> None:
    """Per-component quantities (inlet concentrations, film diffusion …); empty ones left out."""
    entries = dict(data.get(key) or {})

    def store():
        data[key] = entries if entries else (None if optional else {})

    ui.label(label).classes("text-sm text-gray-600")
    with ui.column().classes("pl-4"):
        for name in component_names:
            quantity_input(
                name,
                entries,
                name,
                dimensions,
                on_empty=lambda n=name: entries.pop(n, None),
                on_write=store,
            )


def model_form(
    data: dict, model_class: type[BaseModel], component_names: list[str], skip=()
) -> None:
    """Editors for every field of ``model_class``, reading and writing ``data``."""
    for name, info in model_class.model_fields.items():
        if name in skip:
            continue
        label = name.replace("_", " ")
        optional = _unwrap_optional(info.annotation) is not info.annotation
        annotation = _unwrap_optional(info.annotation)
        origin = typing.get_origin(annotation)
        dimensions = api.field_dimensions(info)
        if dimensions is not None and origin is dict:
            _quantity_dict(label, data, name, dimensions, component_names, optional)
        elif dimensions is not None:
            quantity_input(label, data, name, dimensions)
        elif isinstance(annotation, type) and issubclass(annotation, BaseModel):
            nested = data.get(name) or {}
            data[name] = nested
            with ui.expansion(label).classes("w-full"):
                model_form(nested, annotation, component_names)
        elif isinstance(annotation, type) and issubclass(annotation, enum.Enum):
            options = [member.value for member in annotation]
            ui.select(
                options,
                label=label,
                value=data.get(name),
                on_change=lambda e, n=name: data.__setitem__(n, e.value),
            )
        elif annotation is bool:
            ui.switch(
                label,
                value=bool(data.get(name)),
                on_change=lambda e, n=name: data.__setitem__(n, e.value),
            )
        elif annotation in (int, float):
            ui.number(
                label,
                value=data.get(name),
                on_change=lambda e, n=name, t=annotation: data.__setitem__(
                    n, None if e.value is None else t(e.value)
                ),
            ).classes("w-40")
        elif annotation is str:
            ui.input(
                label,
                value=data.get(name) or "",
                on_change=lambda e, n=name: data.__setitem__(n, e.value),
            ).classes("w-72")
        else:
            ui.label(f"{label}: edit in the project YAML").classes("text-xs text-gray-500")


def list_form(
    items: list, model_class: type[BaseModel], component_names: list[str], title_key: str = "name"
) -> None:
    """Editable list of models (e.g. method steps): one expansion per item, add and remove."""

    @ui.refreshable
    def render():
        for index, item in enumerate(items):
            with ui.expansion(f"{index}: {item.get(title_key, '')}").classes("w-full"):
                model_form(item, model_class, component_names)
                ui.button(
                    "Remove",
                    icon="delete",
                    on_click=lambda i=index: (items.pop(i), render.refresh()),
                ).props("flat color=negative")
        ui.button(
            "Add (copy of last)",
            icon="add",
            on_click=lambda: (
                items.append(copy.deepcopy(items[-1]) if items else {}),
                render.refresh(),
            ),
        ).props("flat")

    render()
