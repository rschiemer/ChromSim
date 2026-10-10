"""Resin (stationary phase) properties."""

from pydantic import BaseModel, Field

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.units import Length, MolarConcentration


class Resin(BaseModel):
    """Chromatography resin packed in the column.

    Attributes
    ----------
    name : str
        Resin name, e.g. ``"Generic CEX"``.
    particle_radius : Quantity
        Mean particle radius r_p, e.g. ``"45 um"``.
    ionic_capacity : Quantity, optional
        Ionic capacity Λ for ion-exchange resins, molar concentration, e.g. ``"1200 mM"``.
    col_porosity : float
        Interstitial (column) porosity ε_c, dimensionless, 0 < ε_c < 1.
    par_porosity : float
        Particle porosity ε_p, dimensionless, 0 < ε_p < 1.
    """

    model_config = MODEL_CONFIG

    name: str
    particle_radius: Length
    ionic_capacity: MolarConcentration | None = None
    col_porosity: float = Field(gt=0.0, lt=1.0)
    par_porosity: float = Field(gt=0.0, lt=1.0)

    @property
    def particle_radius_m(self) -> float:
        """Particle radius r_p, m."""
        return units.to_si(self.particle_radius, "length")

    @property
    def ionic_capacity_mol_per_m3(self) -> float | None:
        """Ionic capacity Λ, mol/m³ (identical to mM), or ``None`` if not set."""
        if self.ionic_capacity is None:
            return None
        return units.to_si(self.ionic_capacity, "molar_concentration")

    @property
    def total_porosity(self) -> float:
        """Total porosity ε_t = ε_c + (1 − ε_c) ε_p, dimensionless.

        Reference: Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and
        Nonlinear Chromatography*, 2nd ed., Academic Press 2006.
        """
        eps_c = self.col_porosity
        eps_p = self.par_porosity
        return eps_c + (1.0 - eps_c) * eps_p
