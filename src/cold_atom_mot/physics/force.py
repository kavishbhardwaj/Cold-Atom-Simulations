"""Effective two-level semiclassical MOT radiation pressure.

This effective model includes local Gaussian intensity, shared saturation, Doppler
shift, a signed scalar Zeeman shift, gravity, and per-beam momentum.  It does
not contain Zeeman populations, Clebsch-Gordan-resolved optical pumping,
coherences, stimulated-force interference, or sub-Doppler forces.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, physical_constants
from ..atomic.species import FineStructureLine
from ..laser.beam import GaussianBeam

BOHR_MAGNETON = physical_constants["Bohr magneton"][0]


@dataclass
class EffectiveMOTForce:
    atom: FineStructureLine
    beams: list[GaussianBeam]
    magnetic_field: object
    gravity: np.ndarray
    effective_magnetic_moment: float = BOHR_MAGNETON

    def scattering_rates(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        position = np.asarray(position, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        intensities = np.array([beam.intensity(position) for beam in self.beams])
        saturation = intensities / self.atom.saturation_intensity_w_m2
        shared_denominator = 1.0 + np.sum(saturation, axis=0)
        magnetic = np.asarray(self.magnetic_field.field(position, time))
        rates = []
        b_magnitude = np.linalg.norm(magnetic, axis=-1)
        for beam, s in zip(self.beams, saturation):
            b_hat = np.divide(magnetic, np.expand_dims(b_magnitude, -1),
                              out=np.zeros_like(magnetic), where=np.expand_dims(b_magnitude > 1e-15, -1))
            # Circular photon angular momentum projected on the local B axis.
            epsilon = beam.polarization
            spin = np.real(1j*np.cross(epsilon, epsilon.conjugate()))
            q_expectation = beam.polarization_purity * np.sum(spin * b_hat, axis=-1)
            zeeman = -self.effective_magnetic_moment * q_expectation * b_magnitude / hbar
            doppler = np.sum(beam.k_vector * velocity, axis=-1)
            delta = beam.detuning + beam.frequency_offset - doppler + zeeman
            # Lorentzian laser linewidth adds optical-coherence dephasing. The
            # configured linewidth is angular FWHM, consistently with detuning.
            width = self.atom.gamma_rad_s + beam.linewidth
            rates.append(0.5 * self.atom.gamma_rad_s * s * self.atom.gamma_rad_s/width /
                         (shared_denominator + (2 * delta / width) ** 2))
        return np.stack(rates, axis=-1)

    def per_beam_force(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        rates = self.scattering_rates(position, velocity, time)
        momenta = hbar * np.stack([beam.k_vector for beam in self.beams])
        return rates[..., :, None] * momenta

    def force(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        optical = np.sum(self.per_beam_force(position, velocity, time), axis=-2)
        return optical + self.atom.mass_kg * np.asarray(self.gravity)

    def linear_coefficients(self, position_step: float = 1e-6, velocity_step: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
        origin = np.zeros(3)
        damping = np.empty(3)
        restoring = np.empty(3)
        for axis in range(3):
            dx = np.zeros(3); dx[axis] = position_step
            dv = np.zeros(3); dv[axis] = velocity_step
            restoring[axis] = -(self.force(dx, origin)[axis] - self.force(-dx, origin)[axis]) / (2 * position_step)
            damping[axis] = -(self.force(origin, dv)[axis] - self.force(origin, -dv)[axis]) / (2 * velocity_step)
        return damping, restoring
