"""Auditable alkali isotope/line registry for solver scope and future extension.

Constants are from D. A. Steck's 85Rb and 87Rb D-line data, revisions 2.3.2
(2021). Registry presence does not imply that every high-fidelity solver supports
the entry; ``model_support`` states that boundary explicitly.
"""

from dataclasses import dataclass
from types import MappingProxyType
from scipy.constants import hbar, k as boltzmann
import numpy as np


@dataclass(frozen=True)
class AtomicLine:
    isotope: str
    line: str
    lower_fine_structure: str
    upper_fine_structure: str
    wavelength_m: float
    lifetime_s: float
    mass_kg: float
    nuclear_spin: float
    ground_hyperfine_hz: float
    natural_abundance: float
    model_support: str
    source: str

    @property
    def gamma_rad_s(self) -> float:
        return 1.0 / self.lifetime_s

    @property
    def wave_number_rad_m(self) -> float:
        return 2 * np.pi / self.wavelength_m

    @property
    def recoil_velocity_m_s(self) -> float:
        return hbar * self.wave_number_rad_m / self.mass_kg

    @property
    def recoil_temperature_k(self) -> float:
        momentum = hbar * self.wave_number_rad_m
        return momentum**2 / (2 * self.mass_kg * boltzmann)

    @property
    def doppler_temperature_k(self) -> float:
        return hbar * self.gamma_rad_s / (2 * boltzmann)


_STECK_87 = "Steck, Rubidium 87 D Line Data, rev. 2.3.2 (2021)"
_STECK_85 = "Steck, Rubidium 85 D Line Data, rev. 2.3.2 (2021)"

ATOMIC_LINES = MappingProxyType({
    ("87Rb", "D2"): AtomicLine(
        "87Rb", "D2", "5S1/2", "5P3/2", 780.241209686e-9,
        26.2348e-9, 1.44316060e-25, 1.5, 6.83468261090429e9, 0.2783,
        "Level A, Level B, and reduced two-level Level C",
        _STECK_87,
    ),
    ("87Rb", "D1"): AtomicLine(
        "87Rb", "D1", "5S1/2", "5P1/2", 794.978851156e-9,
        27.679e-9, 1.44316060e-25, 1.5, 6.83468261090429e9, 0.2783,
        "reduced generic two-level Level C only; hyperfine MOT not implemented",
        _STECK_87,
    ),
    ("85Rb", "D2"): AtomicLine(
        "85Rb", "D2", "5S1/2", "5P3/2", 780.241368271e-9,
        26.2348e-9, 1.409993199e-25, 2.5, 3.0357324390e9, 0.7217,
        "reduced generic two-level Level C only; hyperfine MOT not implemented",
        _STECK_85,
    ),
    ("85Rb", "D1"): AtomicLine(
        "85Rb", "D1", "5S1/2", "5P1/2", 794.9788509e-9,
        27.679e-9, 1.409993199e-25, 2.5, 3.0357324390e9, 0.7217,
        "reduced generic two-level Level C only; hyperfine MOT not implemented",
        _STECK_85,
    ),
})


def get_atomic_line(isotope: str, line: str) -> AtomicLine:
    """Return a registered isotope/line or reject an unsupported identifier."""
    try:
        return ATOMIC_LINES[(isotope, line.upper())]
    except KeyError as error:
        choices = ", ".join(f"{key[0]} {key[1]}" for key in ATOMIC_LINES)
        raise ValueError(f"unknown isotope/line; choose one of: {choices}") from error
