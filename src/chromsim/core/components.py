"""Components of the feed: salt, pH pseudo-component, target protein and impurities.

Conventions (CLAUDE.md): salt is component 0 for SMA-type models; pH, when modelled, is a
non-binding pseudo-component placed directly after salt. Which components bind is decided by the
isotherm (Phase 2).
"""

from pydantic import BaseModel, Field, PositiveFloat

from chromsim import units
from chromsim.core import MODEL_CONFIG


class Component(BaseModel):
    """One component of the mobile phase.

    Attributes
    ----------
    name : str
        Unique name, used as key in method inlets and isotherm parameters.
    mw_kDa : float, optional
        Molecular weight, kDa. Needed to convert mass concentrations (g/L) to mol/m³.
    extinction_coeff_L_per_g_cm : float, optional
        Mass extinction coefficient at the UV wavelength, L/(g·cm). Used for UV prediction.
    charge_hint : float, optional
        Approximate net charge at the working pH, dimensionless. Informational only (e.g. a
        starting guess for the SMA characteristic charge).
    """

    model_config = MODEL_CONFIG

    name: str = Field(min_length=1)
    mw_kDa: PositiveFloat | None = None
    extinction_coeff_L_per_g_cm: PositiveFloat | None = None
    charge_hint: float | None = None

    @property
    def molar_mass_kg_per_mol(self) -> float | None:
        """Molar mass, kg/mol, or ``None`` if ``mw_kDa`` is not set."""
        if self.mw_kDa is None:
            return None
        return units.molar_mass_from_kda(self.mw_kDa)

    @property
    def extinction_coeff_m2_per_kg(self) -> float | None:
        """Mass extinction coefficient, m²/kg, or ``None`` if not set."""
        if self.extinction_coeff_L_per_g_cm is None:
            return None
        quantity = units.Quantity(magnitude=self.extinction_coeff_L_per_g_cm, unit="L/(g*cm)")
        return units.to_si(quantity, "specific_absorptivity")
