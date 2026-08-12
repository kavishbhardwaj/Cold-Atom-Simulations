"""Independent 3D Gaussian travelling-wave laser beams."""

from dataclasses import dataclass
import numpy as np
from .polarization import circular_polarization, unit, transverse_basis, propagate_jones


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
    phase: float = 0.0
    coherence_group: str | None = None
    linewidth: float = 0.0
    label: str = "beam"
    waist_y: float | None = None
    focus_offset: float = 0.0
    jones_vector: np.ndarray | None = None
    optical_elements: tuple = ()
    propagation_mode: str = "collimated"

    def __post_init__(self) -> None:
        if self.power < 0 or self.waist <= 0 or self.wavelength <= 0:
            raise ValueError("power must be non-negative and waist/wavelength positive")
        if self.helicity not in (-1, 1):
            raise ValueError("helicity must be -1 or +1")
        if not 0 <= self.polarization_purity <= 1 or self.linewidth < 0:
            raise ValueError("polarization_purity must be in [0,1]")
        object.__setattr__(self, "direction", unit(self.direction))
        object.__setattr__(self, "origin", np.asarray(self.origin, dtype=float))
        object.__setattr__(self, "waist_y", self.waist if self.waist_y is None else self.waist_y)
        if self.waist_y <= 0 or self.propagation_mode not in ("collimated", "gaussian"):
            raise ValueError("waist_y must be positive and propagation_mode collimated or gaussian")
        if self.jones_vector is not None:
            vector=np.asarray(self.jones_vector,complex)
            if vector.shape != (2,) or np.linalg.norm(vector)==0:
                raise ValueError("jones_vector must be a nonzero two-vector")
            object.__setattr__(self, "jones_vector", vector/np.linalg.norm(vector))

    @property
    def k_vector(self) -> np.ndarray:
        return 2 * np.pi / self.wavelength * self.direction

    @property
    def optical_frequency_hz(self) -> float:
        from scipy.constants import c
        return c/self.wavelength+self.frequency_offset/(2*np.pi)

    @property
    def polarization_ellipticity(self) -> float:
        """Normalized Stokes S3 in the propagation-relative Jones basis."""
        epsilon=self.polarization
        return float(np.dot(np.real(1j*np.cross(epsilon,epsilon.conjugate())),self.direction))

    @property
    def polarization(self) -> np.ndarray:
        if self.jones_vector is not None:
            e1, e2 = transverse_basis(self.direction)
            jones=propagate_jones(self.jones_vector,self.optical_elements)
            return jones[0]*e1+jones[1]*e2
        return circular_polarization(self.direction, self.helicity)

    @property
    def peak_intensity(self) -> float:
        return 2 * self.power / (np.pi * self.waist*self.waist_y)

    @property
    def rayleigh_range(self) -> float:
        return np.pi * self.waist**2 / self.wavelength

    @property
    def rayleigh_range_y(self) -> float:
        return np.pi*self.waist_y**2/self.wavelength

    def waist_at(self, axial_distance):
        return self.waist * np.sqrt(1 + (np.asarray(axial_distance) / self.rayleigh_range) ** 2)

    def intensity(self, positions: np.ndarray) -> np.ndarray:
        """Transverse Gaussian intensity, neglecting diffraction over the MOT volume."""
        points = np.asarray(positions, dtype=float)
        displacement = points - self.origin
        axial = np.sum(displacement * self.direction, axis=-1, keepdims=True)
        transverse = displacement - axial * self.direction
        e1, e2 = transverse_basis(self.direction)
        x, y = np.sum(transverse*e1, axis=-1), np.sum(transverse*e2, axis=-1)
        if self.propagation_mode == "collimated":
            wx, wy = self.waist, self.waist_y
        else:
            z = np.squeeze(axial, axis=-1)-self.focus_offset
            wx = self.waist*np.sqrt(1+(z/self.rayleigh_range)**2)
            wy = self.waist_y*np.sqrt(1+(z/self.rayleigh_range_y)**2)
        return 2*self.power/(np.pi*wx*wy)*np.exp(-2*(x*x/wx**2+y*y/wy**2))

    def complex_field(self, position: np.ndarray) -> np.ndarray:
        """Dimensionless Jones field; interference is handled by coherence group."""
        point = np.asarray(position, float)
        z=np.dot(self.direction,point-self.origin)-self.focus_offset
        extra=0.0
        if self.propagation_mode == "gaussian":
            e1,e2=transverse_basis(self.direction); d=point-self.origin-z*self.direction
            rx=np.inf if z==0 else z*(1+(self.rayleigh_range/z)**2)
            ry=np.inf if z==0 else z*(1+(self.rayleigh_range_y/z)**2)
            extra=self.k_vector.dot(self.direction)*(.5*((d@e1)**2/rx+(d@e2)**2/ry))-np.arctan(z/self.rayleigh_range)/2-np.arctan(z/self.rayleigh_range_y)/2
        return np.sqrt(self.intensity(point)) * self.polarization * np.exp(
            1j * (np.dot(self.k_vector, point - self.origin) + self.phase+extra))


def grouped_intensity(beams: list[GaussianBeam], position: np.ndarray) -> float:
    """Coherently sum within groups and incoherently sum between groups."""
    groups = {}
    total = 0.0
    for index, beam in enumerate(beams):
        if beam.coherence_group is None:
            total += float(beam.intensity(position))
        else:
            groups.setdefault(beam.coherence_group, np.zeros(3, complex))
            groups[beam.coherence_group] += beam.complex_field(position)
    return total + sum(float(np.vdot(field, field).real) for field in groups.values())


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
            # Counter-propagating beams use the same propagation-relative
            # helicity.  Reversing k then reverses their angular momentum in a
            # fixed laboratory/quantization basis, as required by a MOT pair.
            helicity = gradient_sign
            beams.append(GaussianBeam(direction, np.zeros(3), power, waist, detuning, wavelength, helicity, label=f"{propagation_sign:+d}{'xyz'[len(beams)//2]}"))
    return beams
