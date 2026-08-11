"""87Rb D2 constants with explicit provenance and SI units.

Values follow D. A. Steck, *Rubidium 87 D Line Data*, revision 2.3.2
(2021), unless noted.  CODATA exact constants are imported from SciPy.
The saturation intensity is the resonant circularly-polarized closed
|F=2,mF=2> -> |F'=3,mF'=3> convention; other conventions differ.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, k as boltzmann


@dataclass(frozen=True)
class Rb87D2:
    """Atomic data used by the Phase-1 effective cycling-transition model."""

    mass: float = 1.44316060e-25
    wavelength: float = 780.241209686e-9
    lifetime: float = 26.2348e-9
    saturation_intensity: float = 16.69
    ground_hyperfine_splitting_hz: float = 6.83468261090429e9
    # Frequencies relative to 5P3/2 F'=3; sufficient for locating manifolds.
    excited_hyperfine_offsets_hz: tuple[float, ...] = (
        -495.815e6,
        -423.600e6,
        -266.650e6,
        0.0,
    )
    ground_g_f: tuple[float, float] = (-0.501827, 0.499836)
    excited_g_f: tuple[float, float, float, float] = (0.0, 0.666, 0.666, 0.667)
    # Hyperfine line strengths, normalized within each ground F manifold.
    f2_strengths_to_fprime_1_2_3: tuple[float, ...] = (0.05, 0.25, 0.70)
    f1_strengths_to_fprime_0_1_2: tuple[float, ...] = (1 / 6, 5 / 12, 5 / 12)

    @property
    def gamma(self) -> float:
        """Natural linewidth Γ in rad/s (population decay rate)."""
        return 1.0 / self.lifetime

    @property
    def wave_number(self) -> float:
        return 2.0 * np.pi / self.wavelength

    @property
    def recoil_velocity(self) -> float:
        return hbar * self.wave_number / self.mass

    @property
    def recoil_temperature(self) -> float:
        return (hbar * self.wave_number) ** 2 / (2 * self.mass * boltzmann)

    @property
    def doppler_temperature(self) -> float:
        return hbar * self.gamma / (2 * boltzmann)

    def validate(self) -> None:
        if not np.isclose(sum(self.f2_strengths_to_fprime_1_2_3), 1.0):
            raise ValueError("F=2 transition strengths must normalize to one")
        if not np.isclose(sum(self.f1_strengths_to_fprime_0_1_2), 1.0):
            raise ValueError("F=1 transition strengths must normalize to one")
