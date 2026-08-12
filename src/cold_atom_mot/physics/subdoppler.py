"""Phase-resolved adiabatic polarization-gradient cooling.

This module implements a controlled, low-saturation adiabatic-elimination model
for the closed 87Rb D2 F=2 -> F'=3 manifold.  Five ground-state populations are
optically pumped through seven eliminated excited Zeeman states.  It includes
phase-resolved six-beam interference, scalar populations, spatial light shifts,
and their conservative dipole force.  Ground-state coherences are not retained.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, physical_constants
from scipy.integrate import solve_ivp

from ..atomic.species import AtomicBasis
from ..laser.polarization import circular_polarization, spherical_fractions, unit

MU_B = physical_constants["Bohr magneton"][0]


@dataclass(frozen=True)
class CoherentBeam:
    direction: np.ndarray
    polarization: np.ndarray
    saturation: float
    wave_number: float
    phase: float = 0.0
    coherence_group: str = "all"

    def __post_init__(self) -> None:
        if self.saturation < 0 or self.wave_number <= 0:
            raise ValueError("saturation must be non-negative and wave number positive")
        object.__setattr__(self, "direction", unit(self.direction))
        polarization = np.asarray(self.polarization, complex)
        if not np.isclose(np.vdot(polarization, self.direction), 0.0, atol=1e-12):
            raise ValueError("polarization must be transverse")
        object.__setattr__(self, "polarization", polarization / np.linalg.norm(polarization))

    def field(self, position: np.ndarray) -> np.ndarray:
        phase = self.wave_number * np.dot(self.direction, position) + self.phase
        return np.sqrt(self.saturation) * self.polarization * np.exp(1j * phase)


def coherent_six_beam_field(wave_number: float, saturation_per_beam: float, phases=None, coherence_groups=None):
    """Standard phase-resolved sigma+/sigma- MOT geometry.

    Helicity is defined relative to each propagation direction.  The z pair has
    the opposite helicity to x/y for a diag(+,+,-2) quadrupole convention.
    """
    phases = np.zeros(6) if phases is None else np.asarray(phases, float)
    if phases.shape != (6,):
        raise ValueError("six optical phases are required")
    coherence_groups = ["all"] * 6 if coherence_groups is None else list(coherence_groups)
    if len(coherence_groups) != 6:
        raise ValueError("six coherence groups are required")
    beams = []
    for axis_index, (axis, helicity) in enumerate(zip(np.eye(3), (1, 1, -1))):
        for direction_sign in (-1, 1):
            direction = direction_sign * axis
            index = 2 * axis_index + (direction_sign + 1) // 2
            beams.append(CoherentBeam(direction, circular_polarization(direction, helicity),
                                      saturation_per_beam, wave_number, phases[index], coherence_groups[index]))
    return tuple(beams)


class PolarizationGradientModel:
    """Adiabatic F=2 -> F'=3 light-shift and optical-pumping model."""

    ground_m = np.arange(-2, 3)
    excited_m = np.arange(-3, 4)

    def __init__(self, basis: AtomicBasis, ground_f, excited_f, detuning, beams,
                 magnetic_field_t=(0.0, 0.0, 0.0), quantization_axis=(0.0, 0.0, 1.0),
                 *, allow_projected_field=False):
        if detuning == 0:
            raise ValueError("detuning must be non-zero")
        self.basis, self.line = basis, basis.line
        self.gamma, self.detuning, self.wave_number = self.line.gamma_rad_s, detuning, self.line.wave_number_rad_m
        self.beams = tuple(beams)
        self.magnetic_field_t = np.asarray(magnetic_field_t, float)
        self.quantization_axis = unit(quantization_axis)
        transverse = self.magnetic_field_t - np.dot(self.magnetic_field_t, self.quantization_axis)*self.quantization_axis
        if np.linalg.norm(transverse) > 1e-15 and not allow_projected_field:
            raise ValueError("population PGC model supports only B parallel to its fixed quantization axis; use MultilevelOBE for vector B")
        self.allow_projected_field = allow_projected_field
        selected_ground = [(i,s) for i,s in enumerate(basis.ground) if s.F == ground_f]
        selected_excited = [(i,s) for i,s in enumerate(basis.excited) if s.F == excited_f]
        self.ground_m = np.array([s.m for _,s in selected_ground]); self.excited_m = np.array([s.m for _,s in selected_excited])
        self.ground_g_factor = selected_ground[0][1].g_factor; self.excited_g_factor = selected_excited[0][1].g_factor
        ground_map = {old:i for i,(old,_) in enumerate(selected_ground)}; excited_map = {old:i for i,(old,_) in enumerate(selected_excited)}
        self.strength = np.zeros((len(selected_ground), len(selected_excited), 3))
        for transition in basis.transitions:
            if transition.ground_index in ground_map and transition.excited_index in excited_map:
                self.strength[ground_map[transition.ground_index], excited_map[transition.excited_index], transition.q+1] = transition.strength
        self.branching = self.strength.sum(axis=2).T
        self.branching /= self.branching.sum(axis=1, keepdims=True)

    def transition_detunings(self):
        """Laser detuning including the projected linear Zeeman transition shift."""
        b_parallel = np.dot(self.magnetic_field_t, self.quantization_axis)
        magnetic_shift = MU_B * b_parallel / hbar * (
            self.excited_g_factor * self.excited_m[None, :] -
            self.ground_g_factor * self.ground_m[:, None]
        )
        return self.detuning - magnetic_shift

    def electric_field(self, position):
        return sum((beam.field(np.asarray(position, float)) for beam in self.beams), np.zeros(3, complex))

    def grouped_fields(self, position):
        fields = {}
        for beam in self.beams:
            fields.setdefault(beam.coherence_group, np.zeros(3, complex))
            fields[beam.coherence_group] += beam.field(np.asarray(position, float))
        return fields

    def polarization_components(self, position):
        fields = self.grouped_fields(position)
        intensity = sum(float(np.vdot(field, field).real) for field in fields.values())
        if intensity < 1e-30:
            return intensity, {-1: 1 / 3, 0: 1 / 3, 1: 1 / 3}
        weighted = {q: 0.0 for q in (-1,0,1)}
        for field in fields.values():
            group_intensity = float(np.vdot(field,field).real)
            if group_intensity > 1e-30:
                fractions = spherical_fractions(field,self.quantization_axis)
                for q in weighted: weighted[q] += group_intensity*fractions[q]
        return intensity, {q: weighted[q]/intensity for q in weighted}

    def _coupling_saturation(self, position):
        intensity, fractions = self.polarization_components(position)
        q_intensity = intensity * np.array([fractions[q] for q in (-1, 0, 1)])
        return np.sum(self.strength * q_intensity[None, None, :], axis=2)

    def light_shifts(self, position):
        """Ground shifts from adiabatic elimination, including axial Zeeman shifts."""
        detuning = self.transition_detunings()
        denominator = detuning**2 + (self.gamma / 2) ** 2
        ac_stark = hbar * self.gamma**2 * np.sum(detuning * self._coupling_saturation(position) /
                                                (8 * denominator), axis=1)
        b_parallel = np.dot(self.magnetic_field_t, self.quantization_axis)
        return ac_stark + MU_B * self.ground_g_factor * self.ground_m * b_parallel

    def pumping_generator(self, position):
        """Column-conservative ground-population generator in s^-1."""
        coupling = self._coupling_saturation(position)
        detuning = self.transition_detunings()
        denominator = detuning**2 + (self.gamma / 2) ** 2
        excitation = self.gamma**3 * coupling / (8 * denominator)
        generator = np.zeros((len(self.ground_m), len(self.ground_m)))
        for source in range(len(self.ground_m)):
            for excited in range(len(self.excited_m)):
                rate = excitation[source, excited]
                generator[:, source] += rate * self.branching[excited]
                generator[source, source] -= rate
        return generator

    def scattering_rates(self, position):
        """Total eliminated excited-state scattering rate from each ground m."""
        coupling = self._coupling_saturation(position)
        detuning = self.transition_detunings()
        return np.sum(self.gamma**3 * coupling /
                      (8 * (detuning**2 + (self.gamma / 2) ** 2)), axis=1)

    def stationary_populations(self, position):
        generator = self.pumping_generator(position).copy()
        generator[-1] = 1.0
        rhs = np.zeros(len(self.ground_m)); rhs[-1] = 1.0
        populations = np.linalg.solve(generator, rhs)
        return np.maximum(populations, 0) / np.maximum(populations, 0).sum()

    def state_forces(self, position, step=None):
        """Numerical -gradient(U_m), returned as a 5x3 array."""
        step = 2e-5 / self.wave_number if step is None else step
        force = np.empty((len(self.ground_m), 3))
        position = np.asarray(position, float)
        for axis in range(3):
            displacement = np.zeros(3); displacement[axis] = step
            force[:, axis] = -(self.light_shifts(position + displacement) -
                               self.light_shifts(position - displacement)) / (2 * step)
        return force

    def mean_force(self, position, populations):
        return np.asarray(populations) @ self.state_forces(position)

    def moving_average_force(self, velocity, *, periods=24, discard=12, steps_per_period=80,
                             initial_populations=None):
        """Cycle-averaged force for motion along x through the 3D light field."""
        if velocity == 0:
            samples = np.linspace(0, 2 * np.pi / self.wave_number, steps_per_period, endpoint=False)
            return np.mean([self.mean_force([x, 0, 0], self.stationary_populations([x, 0, 0]))[0]
                            for x in samples])
        # Cross-interference between an x beam and the four transverse beams
        # makes the general six-beam field periodic over lambda, not lambda/2.
        spatial_period = 2 * np.pi / self.wave_number
        duration = periods * spatial_period / abs(velocity)
        population0 = np.ones(5) / 5 if initial_populations is None else np.asarray(initial_populations, float)
        def rhs(time, population):
            return self.pumping_generator([velocity * time, 0, 0]) @ population

        samples = np.linspace(0, duration, periods * steps_per_period + 1)
        solution = solve_ivp(rhs, (0, duration), population0, t_eval=samples,
                             rtol=2e-7, atol=1e-10, max_step=duration / (periods * steps_per_period))
        first = discard * steps_per_period
        forces = [self.mean_force([velocity * t, 0, 0], p)[0]
                  for t, p in zip(solution.t[first:], solution.y.T[first:])]
        return float(np.trapezoid(forces, solution.t[first:]) / (solution.t[-1] - solution.t[first]))

    def friction_coefficient(self, velocity, **kwargs):
        force_plus = self.moving_average_force(abs(velocity), **kwargs)
        force_minus = self.moving_average_force(-abs(velocity), **kwargs)
        return -(force_plus - force_minus) / (2 * abs(velocity))

    def recoil_diffusion_tensor(self):
        """Matched recoil-event tensor for this population model.

        Includes absorption shot noise from the configured beam directions and
        isotropic spontaneous emission. It excludes internal-state and dipole-
        force fluctuations, so it must not be combined with coherent-OBE
        friction to claim a quantitative temperature.
        """
        positions = np.linspace(0, 2 * np.pi / self.wave_number, 80, endpoint=False)
        rates = []
        for x in positions:
            p = self.stationary_populations([x, 0, 0])
            rates.append(np.dot(self.scattering_rates([x, 0, 0]), p))
        scattering = np.mean(rates)
        weights = np.asarray([beam.saturation for beam in self.beams], float)
        weights /= weights.sum()
        absorption = sum(weight*np.outer(beam.direction, beam.direction)
                         for weight, beam in zip(weights, self.beams))
        return (hbar*self.wave_number)**2*scattering*(absorption+np.eye(3)/3)/2

    def diffusion_estimate(self):
        """Legacy scalar D_xx from :meth:`recoil_diffusion_tensor`."""
        return self.recoil_diffusion_tensor()[0, 0]
