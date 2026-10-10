"""Project: everything needed to describe one chromatography run (SPEC §4.2)."""

from pydantic import (
    BaseModel,
    Field,
    PositiveFloat,
    PositiveInt,
    field_validator,
    model_validator,
)

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.core.column import Column
from chromsim.core.components import Component
from chromsim.core.method import ResolvedStep, Step, resolve_method
from chromsim.core.resin import Resin
from chromsim.core.system import System

# Version of the project-file format; bump it and add a migration in chromsim/io/project_file.py
# whenever the structure of Project changes incompatibly.
CURRENT_SCHEMA_VERSION = 1


class IsothermConfig(BaseModel):
    """Which isotherm to use and its parameters.

    Attributes
    ----------
    model : str
        Isotherm name, e.g. ``"SMA"``. Checked against the isotherm registry by the backend.
    parameters : dict[str, dict[str, float]]
        Component name → {parameter name → value}. Every binding component has its own set.
        Values are in the units documented by the isotherm (Phase 2).
    """

    model_config = MODEL_CONFIG

    model: str
    parameters: dict[str, dict[str, float]] = Field(default_factory=dict)


class SolverSettings(BaseModel):
    """Numerical settings for the simulation.

    Attributes
    ----------
    backend : str
        Simulation backend: ``"cadet"`` (default) or ``"reference"``.
    n_col : int
        Number of axial cells.
    n_par : int
        Number of radial cells in the particle (GRM only).
    abstol : float
        Absolute solver tolerance, in the units of the state (mol/m³).
    """

    model_config = MODEL_CONFIG

    backend: str = "cadet"
    n_col: PositiveInt = 100
    n_par: PositiveInt = 25
    abstol: PositiveFloat = 1.0e-8


class Project(BaseModel):
    """A complete, self-contained description of a chromatography run.

    Attributes
    ----------
    schema_version : int
        Project-file format version.
    name : str
        Project name.
    components : list of Component
        Components in order. Salt is component 0 for SMA-type models.
    column : Column
        Column geometry, model and transport parameters.
    resin : Resin
        Resin properties.
    system : System
        Extra-column volumes. Default: none.
    isotherm : IsothermConfig
        Isotherm model and per-component parameters.
    method : list of Step
        Method steps in order.
    solver : SolverSettings
        Numerical settings. Default: CADET, 100 axial and 25 particle cells.

    Notes
    -----
    Constructing a Project raises ``pydantic.ValidationError`` listing every problem; use
    :func:`chromsim.io.project_file.load_project` to read files, which reports problems as
    :class:`chromsim.errors.ProjectFileError`. Cross-checks between sections (e.g. that inlet
    components exist) run when the Project is created, not when a single field is reassigned.
    """

    model_config = MODEL_CONFIG

    schema_version: int = CURRENT_SCHEMA_VERSION
    name: str
    components: list[Component] = Field(min_length=1)
    column: Column
    resin: Resin
    system: System = Field(default_factory=System)
    isotherm: IsothermConfig
    method: list[Step] = Field(min_length=1)
    solver: SolverSettings = Field(default_factory=SolverSettings)

    @property
    def component_names(self) -> list[str]:
        """Component names in order."""
        return [component.name for component in self.components]

    @field_validator("schema_version")
    @classmethod
    def _current_schema(cls, version: int) -> int:
        if version != CURRENT_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version {version} is not supported here (expected "
                f"{CURRENT_SCHEMA_VERSION}); read files with load_project() so older versions are "
                "migrated."
            )
        return version

    def resolved_method(self) -> list[ResolvedStep]:
        """The method as an SI timeline; see :func:`chromsim.core.method.resolve_method`."""
        return resolve_method(self.method, self.components, self.column)

    @model_validator(mode="after")
    def _check_cross_references(self) -> "Project":
        names = self.component_names
        known = ", ".join(names)

        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(
                f"Component names must be unique; {', '.join(repr(d) for d in duplicates)} "
                "appears more than once."
            )

        molar_masses = {c.name: c.molar_mass_kg_per_mol for c in self.components}
        for step in self.method:
            for concentrations in (step.inlet, step.inlet_end or {}):
                for component_name, concentration in concentrations.items():
                    if component_name not in molar_masses:
                        raise ValueError(
                            f"Step '{step.name}': unknown component '{component_name}' in inlet. "
                            f"Known components: {known}."
                        )
                    is_mass = (
                        units.check_dimension(
                            concentration, "molar_concentration", "mass_concentration"
                        )
                        == "mass_concentration"
                    )
                    if is_mass and molar_masses[component_name] is None:
                        raise ValueError(
                            f"Step '{step.name}': component '{component_name}' is given as "
                            f"'{concentration}' but has no molecular weight; set mw_kDa for "
                            f"'{component_name}' or give the concentration in mM."
                        )

        for component_name in self.isotherm.parameters:
            if component_name not in molar_masses:
                raise ValueError(
                    f"isotherm parameters: unknown component '{component_name}'. "
                    f"Known components: {known}."
                )

        for field_name in ("film_diffusion", "pore_diffusion"):
            for component_name in getattr(self.column, field_name) or {}:
                if component_name not in molar_masses:
                    raise ValueError(
                        f"column.{field_name}: unknown component '{component_name}'. "
                        f"Known components: {known}."
                    )
        return self
