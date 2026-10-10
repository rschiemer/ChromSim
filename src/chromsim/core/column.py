"""Column geometry, column model type and transport parameters."""

from enum import StrEnum

from pydantic import BaseModel

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.units import Diffusivity, FilmTransferCoefficient, Length


class ColumnModel(StrEnum):
    """Column model family (see SPEC §4; the reference backend supports LRM and LRMP only).

    - ``LRM``: lumped rate model without pores.
    - ``LRMP``: lumped rate model with pores (film transfer into lumped particles).
    - ``GRM``: general rate model (radial diffusion inside the particles). The default.
    """

    LRM = "LRM"
    LRMP = "LRMP"
    GRM = "GRM"


class Column(BaseModel):
    """A packed chromatography column.

    Attributes
    ----------
    model : ColumnModel
        Column model, default GRM.
    length : Quantity
        Packed-bed length L, e.g. ``"20 cm"``.
    diameter : Quantity
        Inner column diameter d, e.g. ``"0.5 cm"``.
    axial_dispersion : Quantity, optional
        Axial dispersion coefficient D_ax, e.g. ``"1e-7 m^2/s"``.
    film_diffusion : dict[str, Quantity], optional
        Film mass-transfer coefficient k_f per component, e.g. ``{"mAb": "1e-5 m/s"}``.
    pore_diffusion : dict[str, Quantity], optional
        Effective pore diffusion coefficient D_p per component, e.g. ``{"mAb": "1e-11 m^2/s"}``.

    Transport parameters are optional here so a project can be set up before calibration (U3);
    the simulation backend reports which ones are missing.
    """

    model_config = MODEL_CONFIG

    model: ColumnModel = ColumnModel.GRM
    length: Length
    diameter: Length
    axial_dispersion: Diffusivity | None = None
    film_diffusion: dict[str, FilmTransferCoefficient] | None = None
    pore_diffusion: dict[str, Diffusivity] | None = None

    @property
    def length_m(self) -> float:
        """Packed-bed length L, m."""
        return units.to_si(self.length, "length")

    @property
    def diameter_m(self) -> float:
        """Inner column diameter d, m."""
        return units.to_si(self.diameter, "length")

    @property
    def cross_section_m2(self) -> float:
        """Column cross-section A = π (d/2)², m²."""
        return units.column_cross_section(self.diameter_m)

    @property
    def volume_m3(self) -> float:
        """Packed-bed volume V_c = A L (one column volume, CV), m³."""
        return units.column_volume(self.diameter_m, self.length_m)
