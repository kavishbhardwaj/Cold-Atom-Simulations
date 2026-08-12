"""Trajectory-derived capture metrics with an explicit numerical criterion."""
from dataclasses import dataclass
import numpy as np
from ..solvers.deterministic import integrate_trajectory
from ..vacuum import (
    one_sided_thermal_flux_m2_s,
    sample_flux_speeds_between,
    sample_spherical_inward_flux,
)

@dataclass(frozen=True)
class CaptureCriterion:
    radius_m: float
    speed_m_s: float
    duration_s: float
    dwell_s: float


@dataclass(frozen=True)
class CaptureEstimate:
    captured: np.ndarray
    capture_time_s: np.ndarray
    initial_speed_m_s: np.ndarray
    capture_probability: float
    capture_probability_standard_error: float
    incident_flux_s: float
    loading_rate_s: float
    seed: int
    sample_weights: np.ndarray
    omitted_high_speed_probability: float = 0.0

def evaluate_capture(force_model, positions, velocities, criterion, *, max_step):
    captured=[]; capture_time=[]
    for r,v in zip(np.asarray(positions),np.asarray(velocities)):
        tr=integrate_trajectory(force_model,r,v,criterion.duration_s,max_step=max_step)
        inside=(np.linalg.norm(tr.position,axis=1)<=criterion.radius_m)&(np.linalg.norm(tr.velocity,axis=1)<=criterion.speed_m_s)
        # Adaptive RK45 output is nonuniform: evaluate actual contiguous time,
        # not an assumed number of samples based on a median timestep.
        changes = np.diff(np.r_[False, inside, False].astype(int))
        starts = np.flatnonzero(changes == 1)
        stops = np.flatnonzero(changes == -1) - 1
        dwell = tr.time[stops] - tr.time[starts] if len(starts) else np.array([])
        qualifying = np.flatnonzero(dwell >= criterion.dwell_s)
        ok = bool(len(qualifying))
        captured.append(ok)
        capture_time.append(float(tr.time[starts[qualifying[0]]]) if ok else np.nan)
    return np.asarray(captured),np.asarray(capture_time)


def estimate_vapor_capture_rate(
    force_model,
    vapor_state,
    criterion,
    *,
    capture_surface_radius_m,
    atoms,
    max_step_s,
    seed,
):
    """Estimate isotope-specific MOT loading from thermal flux and trajectories.

    The spherical surface is an acceptance boundary, not a chamber-wall model.
    Every sampled trajectory has equal statistical weight because positions,
    incidence angles, and speeds are drawn from their flux distributions.
    """
    if atoms <= 0:
        raise ValueError("atoms must be positive")
    positions, velocities, speeds = sample_spherical_inward_flux(
        capture_surface_radius_m,
        vapor_state.temperature_k,
        force_model.atom.mass_kg,
        atoms,
        seed=seed,
    )
    captured, capture_time = evaluate_capture(
        force_model, positions, velocities, criterion, max_step=max_step_s
    )
    probability = float(captured.mean())
    standard_error = float(np.sqrt(probability * (1 - probability) / atoms))
    density = vapor_state.isotope_number_density_m3(force_model.atom.isotope)
    area = 4 * np.pi * capture_surface_radius_m**2
    incident_flux = area * one_sided_thermal_flux_m2_s(
        density, vapor_state.temperature_k, force_model.atom.mass_kg
    )
    return CaptureEstimate(
        captured,
        capture_time,
        speeds,
        probability,
        standard_error,
        float(incident_flux),
        float(incident_flux * probability),
        seed,
        np.full(atoms, 1 / atoms),
    )


def estimate_stratified_vapor_capture_rate(
    force_model,
    vapor_state,
    criterion,
    *,
    capture_surface_radius_m,
    speed_bin_edges_m_s,
    atoms_per_bin,
    max_step_s,
    seed,
):
    """Resolve the rare slow tail with flux-probability-weighted speed strata.

    Speeds above the last edge are conservatively assigned zero capture. The
    omitted probability is returned, so callers must demonstrate that their
    last edge exceeds the calculated capture range rather than hiding a cutoff.
    """
    edges = np.asarray(speed_bin_edges_m_s, float)
    if atoms_per_bin <= 0 or edges.ndim != 1 or len(edges) < 2:
        raise ValueError("positive atoms_per_bin and at least two speed edges required")
    if edges[0] != 0 or np.any(np.diff(edges) <= 0):
        raise ValueError("speed edges must start at zero and increase")
    rng = np.random.default_rng(seed)
    captured_parts, time_parts, speed_parts, weight_parts = [], [], [], []
    for low, high in zip(edges[:-1], edges[1:]):
        # Reuse the cosine-law geometry sampler, replacing only its speed draw.
        child_seed = int(rng.integers(0, 2**32 - 1))
        positions, directions, _ = sample_spherical_inward_flux(
            capture_surface_radius_m,
            vapor_state.temperature_k,
            force_model.atom.mass_kg,
            atoms_per_bin,
            seed=child_seed,
        )
        unit_directions = directions / np.linalg.norm(directions, axis=1)[:, None]
        speeds, bin_probability = sample_flux_speeds_between(
            low,
            high,
            vapor_state.temperature_k,
            force_model.atom.mass_kg,
            atoms_per_bin,
            rng,
        )
        caught, times = evaluate_capture(
            force_model,
            positions,
            speeds[:, None] * unit_directions,
            criterion,
            max_step=max_step_s,
        )
        captured_parts.append(caught)
        time_parts.append(times)
        speed_parts.append(speeds)
        weight_parts.append(np.full(atoms_per_bin, bin_probability / atoms_per_bin))
    captured = np.concatenate(captured_parts)
    times = np.concatenate(time_parts)
    speeds = np.concatenate(speed_parts)
    weights = np.concatenate(weight_parts)
    probability = float(np.dot(weights, captured))
    # Independent binomial strata: sum_j w_j² p_j(1-p_j)/n_j.
    variance = 0.0
    for caught, sample_weights in zip(captured_parts, weight_parts):
        bin_weight = sample_weights.sum()
        p_bin = caught.mean()
        variance += bin_weight**2 * p_bin * (1 - p_bin) / atoms_per_bin
    omitted = float(1 - weights.sum())
    density = vapor_state.isotope_number_density_m3(force_model.atom.isotope)
    area = 4 * np.pi * capture_surface_radius_m**2
    incident_flux = area * one_sided_thermal_flux_m2_s(
        density, vapor_state.temperature_k, force_model.atom.mass_kg
    )
    return CaptureEstimate(
        captured,
        times,
        speeds,
        probability,
        float(np.sqrt(variance)),
        float(incident_flux),
        float(incident_flux * probability),
        seed,
        weights,
        omitted,
    )
