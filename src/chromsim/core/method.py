"""Chromatography method: an ordered list of steps, and its conversion to an SI timeline.

A step has a duration (time, volume or column volumes), a flow (linear velocity or volumetric
flow) and inlet concentrations. A linear gradient is a step with ``inlet_end``: each component
listed there changes linearly from its ``inlet`` value to its ``inlet_end`` value over the step.
"""

from dataclasses import dataclass

from pydantic import BaseModel, model_validator

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.core.column import Column
from chromsim.core.components import Component
from chromsim.units import Concentration, Duration, FlowRate, Quantity


class Step(BaseModel):
    """One step of a chromatography method.

    Attributes
    ----------
    name : str
        Step name, e.g. ``"load"``.
    duration : Quantity
        Step length as time (``"10 min"``), volume (``"2 mL"``) or column volumes (``"5 CV"``).
    flow : Quantity
        Linear velocity (``"150 cm/h"``) or volumetric flow (``"0.5 mL/min"``).
    inlet : dict[str, Quantity]
        Component name → inlet concentration at the start of the step (mM or g/L).
        Components not listed are 0.
    inlet_end : dict[str, Quantity], optional
        For linear gradients: component name → concentration at the end of the step. Every
        component listed here must also be in ``inlet``.
    """

    model_config = MODEL_CONFIG

    name: str
    duration: Duration
    flow: FlowRate
    inlet: dict[str, Concentration]
    inlet_end: dict[str, Concentration] | None = None

    @model_validator(mode="after")
    def _gradient_has_start(self) -> "Step":
        for component_name in self.inlet_end or {}:
            if component_name not in self.inlet:
                raise ValueError(
                    f"Step '{self.name}': inlet_end gives '{component_name}' but inlet does not; "
                    f"add '{component_name}' to inlet with its start concentration."
                )
        return self

    @property
    def is_gradient(self) -> bool:
        """``True`` if any inlet concentration changes during the step."""
        return bool(self.inlet_end)


@dataclass(frozen=True)
class ResolvedStep:
    """A method step converted to SI units on an absolute time axis.

    Attributes
    ----------
    name : str
        Step name.
    start_s, end_s : float
        Start and end time from the beginning of the method, s.
    flow_rate_m3_per_s : float
        Volumetric flow rate, m³/s.
    inlet_start_mol_per_m3, inlet_end_mol_per_m3 : tuple of float
        Inlet concentration of every project component (in project order) at the start and end
        of the step, mol/m³. Equal for isocratic steps.
    """

    name: str
    start_s: float
    end_s: float
    flow_rate_m3_per_s: float
    inlet_start_mol_per_m3: tuple[float, ...]
    inlet_end_mol_per_m3: tuple[float, ...]

    @property
    def duration_s(self) -> float:
        """Step duration, s."""
        return self.end_s - self.start_s


def _inlet_vector(
    concentrations: dict[str, Quantity], components: list[Component]
) -> tuple[float, ...]:
    """Inlet concentrations of all components in project order, mol/m³ (unlisted → 0)."""
    vector = []
    for component in components:
        if component.name in concentrations:
            vector.append(
                units.concentration_to_si(
                    concentrations[component.name], component.molar_mass_kg_per_mol
                )
            )
        else:
            vector.append(0.0)
    return tuple(vector)


def resolve_method(
    steps: list[Step], components: list[Component], column: Column
) -> list[ResolvedStep]:
    """Convert method steps to SI start/end times, flow rates and inlet concentrations.

    Parameters
    ----------
    steps : list of Step
        The method, in order.
    components : list of Component
        Project components; defines the order of the concentration tuples.
    column : Column
        Column, for the cross-section (cm/h → m³/s) and column volume (CV → s).

    Returns
    -------
    list of ResolvedStep
        One entry per step; each step starts where the previous one ends, the first at 0 s.
    """
    resolved = []
    start_s = 0.0
    for step in steps:
        flow_rate = units.flow_rate_to_si(step.flow, column.cross_section_m2)
        duration_s = units.duration_to_si(step.duration, flow_rate, column.volume_m3)
        inlet_end = {**step.inlet, **(step.inlet_end or {})}
        resolved.append(
            ResolvedStep(
                name=step.name,
                start_s=start_s,
                end_s=start_s + duration_s,
                flow_rate_m3_per_s=flow_rate,
                inlet_start_mol_per_m3=_inlet_vector(step.inlet, components),
                inlet_end_mol_per_m3=_inlet_vector(inlet_end, components),
            )
        )
        start_s += duration_s
    return resolved
