"""Incoherent hyperfine/Zeeman rate equations for supported D2 MOTs.

This model evolves populations only.  Off-diagonal density-matrix coherences,
coherent dark states, light shifts and polarization-gradient forces require the
future OBE/Level-C and phase-resolved/Level-D models.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, physical_constants
from ..atomic.species import AtomicBasis, FineStructureLine
from ..laser.beam import GaussianBeam
from ..laser.polarization import spherical_fractions

BOHR_MAGNETON = physical_constants["Bohr magneton"][0]


@dataclass(frozen=True)
class BeamFamily:
    beam: GaussianBeam
    ground_f: int
    target_excited_f: int
    name: str


@dataclass
class MultilevelRateEquationMOT:
    atom: FineStructureLine
    basis: AtomicBasis
    beam_families: list[BeamFamily]
    magnetic_field: object
    gravity: np.ndarray

    def _local_components(self, beam: GaussianBeam, magnetic: np.ndarray) -> dict[int, float]:
        magnitude = np.linalg.norm(magnetic)
        axis = magnetic / magnitude if magnitude > 1e-12 else np.array([0.0, 0.0, 1.0])
        ideal = spherical_fractions(beam.polarization, axis)
        return {q: beam.polarization_purity * ideal[q] + (1 - beam.polarization_purity) / 3 for q in (-1, 0, 1)}

    def stimulated_rates(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        """Return W[beam,transition] in s^-1 for one phase-space point."""
        position = np.asarray(position, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        if position.shape != (3,) or velocity.shape != (3,):
            raise ValueError("Level-B rates currently accept one 3D phase-space point")
        magnetic = np.asarray(self.magnetic_field.field(position, time), dtype=float)
        b_magnitude = np.linalg.norm(magnetic)
        saturation = np.array([family.beam.intensity(position) / self.atom.saturation_intensity_w_m2 for family in self.beam_families])
        output = np.zeros((len(self.beam_families), len(self.basis.transitions)))
        for bi, family in enumerate(self.beam_families):
            components = self._local_components(family.beam, magnetic)
            target_offset = (self.atom.hyperfine_energy_hz("excited", family.target_excited_f) -
                             self.atom.hyperfine_energy_hz("excited", max(self.atom.excited_f)))
            for ti, transition in enumerate(self.basis.transitions):
                ground = self.basis.ground[transition.ground_index]
                excited = self.basis.excited[transition.excited_index]
                if ground.F != family.ground_f:
                    continue
                hyperfine = -2 * np.pi * (excited.frequency_offset_hz - target_offset)
                zeeman = BOHR_MAGNETON / hbar * (excited.g_factor * excited.m - ground.g_factor * ground.m) * b_magnitude
                doppler = np.dot(family.beam.k_vector, velocity)
                delta = family.beam.detuning + family.beam.frequency_offset + hyperfine - zeeman - doppler
                effective_s = saturation[bi] * transition.strength * components[transition.q]
                # Saturation emerges from the bidirectional stimulated terms
                # and finite populations in the rate equations. Adding the
                # Level-A shared denominator here would count it twice.
                output[bi, ti] = 0.5 * self.atom.gamma_rad_s * effective_s / (1.0 + (2 * delta / self.atom.gamma_rad_s) ** 2)
        return output

    def generator(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        """Return population generator A with dp/dt=A p (column convention)."""
        ng = len(self.basis.ground)
        matrix = np.zeros((self.basis.state_count, self.basis.state_count))
        summed = self.stimulated_rates(position, velocity, time).sum(axis=0)
        for transition, rate in zip(self.basis.transitions, summed):
            gi = transition.ground_index
            ei = ng + transition.excited_index
            matrix[ei, gi] += rate
            matrix[gi, gi] -= rate
            matrix[gi, ei] += rate
            matrix[ei, ei] -= rate
        for ei, row in enumerate(self.basis.spontaneous_branching):
            source = ng + ei
            for gi, probability in enumerate(row):
                matrix[gi, source] += self.atom.gamma_rad_s * probability
            matrix[source, source] -= self.atom.gamma_rad_s
        return matrix

    def steady_state(self, position: np.ndarray, velocity: np.ndarray, time: float = 0.0) -> np.ndarray:
        matrix = self.generator(position, velocity, time).copy()
        rhs = np.zeros(self.basis.state_count)
        matrix[-1, :] = 1.0
        rhs[-1] = 1.0
        population = np.linalg.solve(matrix, rhs)
        population[np.abs(population) < 1e-14] = 0.0
        if np.min(population) < -1e-9:
            raise RuntimeError("rate-equation steady state has negative population")
        return np.maximum(population, 0) / np.maximum(population, 0).sum()

    def per_beam_force(self, position: np.ndarray, velocity: np.ndarray, population: np.ndarray | None = None, time: float = 0.0) -> np.ndarray:
        population = self.steady_state(position, velocity, time) if population is None else np.asarray(population)
        ng = len(self.basis.ground)
        rates = self.stimulated_rates(position, velocity, time)
        force = np.zeros((len(self.beam_families), 3))
        for bi, family in enumerate(self.beam_families):
            net_rate = 0.0
            for transition, rate in zip(self.basis.transitions, rates[bi]):
                net_rate += rate * (population[transition.ground_index] - population[ng + transition.excited_index])
            force[bi] = hbar * family.beam.k_vector * net_rate
        return force

    def force(self, position: np.ndarray, velocity: np.ndarray, population: np.ndarray | None = None, time: float = 0.0) -> np.ndarray:
        return self.per_beam_force(position, velocity, population, time).sum(axis=0) + self.atom.mass_kg * np.asarray(self.gravity)

    def manifold_populations(self, population: np.ndarray) -> dict[str, float]:
        ng = len(self.basis.ground)
        return {
            **{f"ground_F{f}": float(sum(population[i] for i, state in enumerate(self.basis.ground) if state.F == f)) for f in self.atom.ground_f},
            "excited": float(population[ng:].sum()),
        }
