"""Independent 3D Gaussian travelling-wave laser beams."""

from dataclasses import dataclass
import numpy as np
from .polarization import circular_polarization, unit


@dataclass(frozen=True)
class GaussianBeam:
    direction: np.ndarray
    origin: np.ndarray
    power: float
    waist: float
    detuning: float
    wavelength: float
    helicity: int
    frequency_offset: float = 0.0
    polarization_purity: float = 1.0
    label: str = "beam"

    def __post_init__(self) -> None:
        if self.power < 0 or self.waist <= 0 or self.wavelength <= 0:
            raise ValueError("power must be non-negative and waist/wavelength positive")
        if self.helicity not in (-1, 1):
            raise ValueError("helicity must be -1 or +1")
        if not 0 <= self.polarization_purity <= 1:
            raise ValueError("polarization_purity must be in [0,1]")
        object.__setattr__(self, "direction", unit(self.direction))
        object.__setattr__(self, "origin", np.asarray(self.origin, dtype=float))

    @property
    def k_vector(self) -> np.ndarray:
        return 2 * np.pi / self.wavelength * self.direction

    @property
    def polarization(self) -> np.ndarray:
        return circular_polarization(self.direction, self.helicity)

    @property
    def peak_intensity(self) -> float:
        return 2 * self.power / (np.pi * self.waist**2)

    def intensity(self, positions: np.ndarray) -> np.ndarray:
        """Transverse Gaussian intensity, neglecting diffraction over the MOT volume."""
        points = np.asarray(positions, dtype=float)
        displacement = points - self.origin
        axial = np.sum(displacement * self.direction, axis=-1, keepdims=True)
        transverse = displacement - axial * self.direction
        radius_squared = np.sum(transverse**2, axis=-1)
        return self.peak_intensity * np.exp(-2 * radius_squared / self.waist**2)


def six_beam_mot(power: float, waist: float, detuning: float, wavelength: float) -> list[GaussianBeam]:
    """Return six independent beams with helicity tied to each propagation vector.

    The signs correspond to a quadrupole gradient diag(+,+,-2) and the force
    module's signed Zeeman convention.  Each object remains independently
    replaceable for imbalance, displacement, or pointing studies.
    """
    beams = []
    gradient_signs = (1, 1, -1)
    axes = np.eye(3)
    for axis, gradient_sign in zip(axes, gradient_signs):
        for propagation_sign in (-1, 1):
            direction = propagation_sign * axis
            helicity = -propagation_sign * gradient_sign
            beams.append(GaussianBeam(direction, np.zeros(3), power, waist, detuning, wavelength, helicity, label=f"{propagation_sign:+d}{'xyz'[len(beams)//2]}"))
    return beams
