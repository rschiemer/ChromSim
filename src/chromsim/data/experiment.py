"""Experimental chromatograms: detector signals, fraction marks and run log (SPEC §7.2).

An :class:`Experiment` is what every chromatogram importer returns. Each signal keeps its own
axis, because ÄKTA systems record UV, conductivity and pH at different rates. All axes and values
are in ChromSim's internal units (volume in m³, time in s, absorbance in AU, conductivity in
S/m); the units of the export are kept for reference and display.
"""

from dataclasses import dataclass, field, replace

import numpy as np

from chromsim import units
from chromsim.errors import ChromSimError

SIGNAL_KINDS = (
    "uv",
    "conductivity",
    "ph",
    "concentration_b",
    "pressure",
    "concentration",
    "salt",
    "other",
)


@dataclass(frozen=True)
class Signal:
    """One recorded (or predicted) detector curve.

    Attributes
    ----------
    name : str
        Curve name as in the export, e.g. ``"UV1_280nm"``.
    kind : str
        What the curve measures: one of ``SIGNAL_KINDS``.
    x_axis : str
        ``"volume"`` (x in m³) or ``"time"`` (x in s).
    x : numpy.ndarray
        Axis values, m³ or s, increasing.
    values : numpy.ndarray
        Signal values in ``unit``.
    unit : str
        Internal unit of ``values``: ``"AU"``, ``"S/m"``, ``"pH"``, ``"1"`` (fraction),
        ``"mol/m**3"``, ``"kg/m**3"`` …
    original_unit : str
        Unit of the values in the export, e.g. ``"mAU"``; used as the default display unit.
    """

    name: str
    kind: str
    x_axis: str
    x: np.ndarray
    values: np.ndarray
    unit: str
    original_unit: str

    def values_in(self, unit: str) -> np.ndarray:
        """Values converted to ``unit`` for display, e.g. ``"mAU"`` or ``"mS/cm"``."""
        return units.signal_from_si(self.values, self.unit, unit)

    def shifted(self, delta_x: float) -> "Signal":
        """The same signal with its axis moved by ``delta_x`` (m³ or s)."""
        return replace(self, x=self.x + delta_x)


@dataclass(frozen=True)
class Fraction:
    """A collected fraction.

    Attributes
    ----------
    label : str
        Fraction id as written by the fraction collector, e.g. ``"A3"``.
    start_m3, end_m3 : float
        Start and end of the fraction on the run's volume axis, m³.
    """

    label: str
    start_m3: float
    end_m3: float

    @property
    def volume_m3(self) -> float:
        """Fraction volume, m³."""
        return self.end_m3 - self.start_m3


# Labels of fraction-collector marks that are not fractions.
NON_FRACTION_LABELS = frozenset({"waste", "frac"})


def fractions_from_marks(marks: list[tuple[float, str]], end_m3: float) -> list[Fraction]:
    """Turn fraction-collector marks into fractions.

    Each mark starts a fraction that ends at the next mark. Marks labelled ``Waste`` end the
    previous fraction without starting a new one. The last fraction ends at ``end_m3``.

    Parameters
    ----------
    marks : list of (float, str)
        (volume in m³, label), in any order.
    end_m3 : float
        Volume at which the last fraction ends (normally the end of the run), m³.

    Returns
    -------
    list of Fraction
        Fractions in collection order.
    """
    ordered = sorted(marks, key=lambda mark: mark[0])
    fractions = []
    for index, (start, label) in enumerate(ordered):
        if label.strip().strip('"').lower() in NON_FRACTION_LABELS:
            continue
        end = ordered[index + 1][0] if index + 1 < len(ordered) else end_m3
        fractions.append(Fraction(label=label.strip().strip('"'), start_m3=start, end_m3=end))
    return fractions


@dataclass
class Experiment:
    """An imported chromatography run.

    Attributes
    ----------
    name : str
        Run name (default: file name).
    signals : dict[str, Signal]
        Curve name → signal.
    fractions : list of Fraction
        Collected fractions, in order.
    injections_m3 : list of float
        Volumes of injection marks, m³.
    log : list of (float, str)
        Run-log entries: (volume in m³, text).
    metadata : dict
        Source file, importer, SHA-256 of the file, import time and user.
    """

    name: str
    signals: dict[str, Signal]
    fractions: list[Fraction] = field(default_factory=list)
    injections_m3: list[float] = field(default_factory=list)
    log: list[tuple[float, str]] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    def signal(self, name_or_kind: str) -> Signal:
        """Find a signal by its exact name, or else the first signal of that kind.

        Parameters
        ----------
        name_or_kind : str
            A curve name (``"UV1_280nm"``) or a kind (``"uv"``, ``"conductivity"`` …).

        Raises
        ------
        ChromSimError
            If nothing matches; the message lists the available signals.
        """
        if name_or_kind in self.signals:
            return self.signals[name_or_kind]
        for signal in self.signals.values():
            if signal.kind == name_or_kind:
                return signal
        available = ", ".join(f"{s.name} ({s.kind})" for s in self.signals.values())
        raise ChromSimError(
            f"Experiment '{self.name}' has no signal '{name_or_kind}'. Available: {available}."
        )

    def fraction(self, label: str) -> Fraction:
        """The fraction with this label."""
        for fraction in self.fractions:
            if fraction.label == label:
                return fraction
        raise ChromSimError(f"Experiment '{self.name}' has no fraction '{label}'.")

    @property
    def end_m3(self) -> float:
        """Largest volume recorded by any volume-based signal, m³."""
        volume_ends = [float(np.max(s.x)) for s in self.signals.values() if s.x_axis == "volume"]
        return max(volume_ends) if volume_ends else 0.0
