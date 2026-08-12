"""Discrete photon-event Monte Carlo trajectories with isotropic recoil."""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar


@dataclass(frozen=True)
class EnsembleTrajectory:
    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    scattering_events: int
    seed: int


def isotropic_directions(rng: np.random.Generator, count: int) -> np.ndarray:
    """Sample directions uniformly on the unit sphere."""
    vectors = rng.normal(size=(count, 3))
    return vectors / np.linalg.norm(vectors, axis=1)[:, None]


def simulate_photon_events(force_model, position, velocity, duration, time_step, *, seed=0, store_every=1) -> EnsembleTrajectory:
    """Evolve atoms using Bernoulli absorption and isotropic emission events.

    At most one absorption is permitted per atom and step.  A run is rejected
    when the maximum total event probability exceeds 0.1, rather than silently
    entering the multiple-event regime.
    """
    positions = np.atleast_2d(np.asarray(position, dtype=float)).copy()
    velocities = np.atleast_2d(np.asarray(velocity, dtype=float)).copy()
    if positions.shape != velocities.shape or positions.shape[1] != 3:
        raise ValueError("position and velocity must have matching (N,3) shapes")
    steps = int(np.ceil(duration / time_step))
    if steps <= 0 or time_step <= 0:
        raise ValueError("duration and time_step must be positive")
    rng = np.random.default_rng(seed)
    times, saved_positions, saved_velocities = [], [], []
    events = 0
    recoil = hbar * force_model.atom.wave_number_rad_m / force_model.atom.mass_kg
    for step in range(steps + 1):
        time = min(step * time_step, duration)
        if step % store_every == 0 or step == steps:
            times.append(time); saved_positions.append(positions.copy()); saved_velocities.append(velocities.copy())
        if step == steps:
            break
        dt = min(time_step, duration - time)
        rates = force_model.scattering_rates(positions, velocities, time)
        total_rates = rates.sum(axis=1)
        event_probability = -np.expm1(-total_rates * dt)
        if np.max(event_probability) > 0.1:
            raise ValueError("time step too large: total scattering probability exceeds 0.1")
        scattered = rng.random(len(positions)) < event_probability
        indices = np.flatnonzero(scattered)
        for atom_index in indices:
            probabilities = rates[atom_index] / total_rates[atom_index]
            beam_index = rng.choice(len(force_model.beams), p=probabilities)
            velocities[atom_index] += recoil * force_model.beams[beam_index].direction
        if len(indices):
            # The emitted photon carries +hbar*k*n, so the atom receives -hbar*k*n.
            velocities[indices] -= recoil * isotropic_directions(rng, len(indices))
        events += len(indices)
        velocities += np.asarray(force_model.gravity) * dt
        positions += velocities * dt
    return EnsembleTrajectory(np.array(times), np.array(saved_positions), np.array(saved_velocities), events, seed)
