"""Trajectory-derived capture metrics with an explicit numerical criterion."""
from dataclasses import dataclass
import numpy as np
from ..solvers.deterministic import integrate_trajectory
from ..vacuum import (
    one_sided_thermal_flux_m2_s,
    sample_flux_speeds_between,
    sample_spherical_inward_flux,
    wilson_interval,
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
    confidence_interval: tuple[float, float] = (0.0, 1.0)
    last_simulated_speed_m_s: float = 0.0
    omitted_capture_probability_upper: float = 1.0
    omitted_loading_rate_upper_s: float = np.inf
    tail_converged: bool = False

    @property
    def loading_rate_confidence_interval_s(self):
        """Flux-scaled finite-sample interval for the loading-rate estimate."""
        return tuple(self.incident_flux_s * value for value in self.confidence_interval)

def evaluate_capture(force_model, positions, velocities, criterion, *, max_step,
                     rtol=1e-7, atol=1e-10, sample_step=None):
    captured=[]; capture_time=[]
    for r,v in zip(np.asarray(positions),np.asarray(velocities)):
        tr=integrate_trajectory(force_model,r,v,criterion.duration_s,max_step=max_step,
                                rtol=rtol,atol=atol,
                                sample_step=sample_step or min(max_step/4, criterion.dwell_s/10))
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
        vapor_state.vapor_temperature_k,
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
        density, vapor_state.vapor_temperature_k, force_model.atom.mass_kg
    )
    confidence_interval = wilson_interval(int(captured.sum()), atoms)
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
        0.0,
        confidence_interval,
        float(np.max(speeds)),
        0.0,
        0.0,
        True,
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
    confidence_level=0.95,
    rtol=1e-7,
    atol=1e-10,
):
    """Resolve the rare slow tail with flux-probability-weighted speed strata.

    This fixed-edge primitive reports, but does not estimate, the omitted tail.
    Use :func:`estimate_adaptive_vapor_capture_rate` when a statistically bounded
    loading estimate is required; callers must never interpret omission as proof
    of exactly zero high-speed capture.
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
            vapor_state.vapor_temperature_k,
            force_model.atom.mass_kg,
            atoms_per_bin,
            seed=child_seed,
        )
        unit_directions = directions / np.linalg.norm(directions, axis=1)[:, None]
        speeds, bin_probability = sample_flux_speeds_between(
            low,
            high,
            vapor_state.vapor_temperature_k,
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
            rtol=rtol,
            atol=atol,
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
    # Independent stratified variance plus conservative weighted Wilson bounds.
    variance = 0.0
    confidence_low = 0.0
    confidence_high = 0.0
    for caught, sample_weights in zip(captured_parts, weight_parts):
        bin_weight = sample_weights.sum()
        p_bin = caught.mean()
        variance += bin_weight**2 * p_bin * (1 - p_bin) / atoms_per_bin
        low, high = wilson_interval(int(caught.sum()), atoms_per_bin, confidence_level)
        confidence_low += bin_weight * low
        confidence_high += bin_weight * high
    omitted = float(1 - weights.sum())
    density = vapor_state.isotope_number_density_m3(force_model.atom.isotope)
    area = 4 * np.pi * capture_surface_radius_m**2
    incident_flux = area * one_sided_thermal_flux_m2_s(
        density, vapor_state.vapor_temperature_k, force_model.atom.mass_kg
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
        (float(confidence_low), float(confidence_high)),
        float(edges[-1]),
        1.0,
        float(incident_flux * omitted),
        False,
    )


def estimate_adaptive_vapor_capture_rate(
    force_model, vapor_state, criterion, *, capture_surface_radius_m,
    initial_speed_edges_m_s, atoms_per_bin, maximum_speed_m_s,
    tail_relative_loading_tolerance, max_step_s, seed, confidence_level=0.95,
    rtol=1e-7, atol=1e-10,
):
    """Extend geometric speed strata until the unresolved loading tail is bounded.

    Stopping requires a zero-capture final stratum and a Wilson upper bound on
    its capture probability times the omitted flux to be below the requested
    fraction of the estimated captured probability. This is a statistical model
    bound, not a proof that faster atoms have zero capture.
    """
    edges=list(map(float,initial_speed_edges_m_s))
    # Preconstruct progressive geometric tail strata, then simulate each exactly
    # once. This is equivalent to sequential extension but avoids resampling all
    # lower-speed strata whenever a tail edge is added.
    while edges[-1] < maximum_speed_m_s:
        edges.append(min(maximum_speed_m_s,edges[-1]*2))
    estimate=estimate_stratified_vapor_capture_rate(
        force_model,vapor_state,criterion,
        capture_surface_radius_m=capture_surface_radius_m,
        speed_bin_edges_m_s=edges,atoms_per_bin=atoms_per_bin,
        max_step_s=max_step_s,seed=seed,confidence_level=confidence_level,
        rtol=rtol,atol=atol,
    )
    # Apply the sequential stopping rule to the independently generated strata.
    # Pre-generation is only a performance detail; samples above the selected
    # stopping edge are discarded from the reported estimator.
    stop = len(edges) - 1
    converged = False
    cumulative_probability = 0.0
    cumulative_weight = 0.0
    omitted_capture_upper = 1.0
    for index in range(len(edges) - 1):
        block = slice(index * atoms_per_bin, (index + 1) * atoms_per_bin)
        caught = estimate.captured[block]
        bin_weight = float(estimate.sample_weights[block].sum())
        cumulative_weight += bin_weight
        cumulative_probability += float(np.dot(estimate.sample_weights[block], caught))
        _, upper = wilson_interval(int(caught.sum()), atoms_per_bin, confidence_level)
        omitted_capture_upper = (1.0 - cumulative_weight) * upper
        relative = omitted_capture_upper / max(cumulative_probability, 1e-300)
        if caught.sum() == 0 and relative <= tail_relative_loading_tolerance:
            stop = index + 1
            converged = True
            break
    count = stop * atoms_per_bin
    selected_captured = estimate.captured[:count]
    selected_weights = estimate.sample_weights[:count]
    selected_times = estimate.capture_time_s[:count]
    selected_speeds = estimate.initial_speed_m_s[:count]
    probability = float(np.dot(selected_weights, selected_captured))
    variance = 0.0
    confidence_low = confidence_high = 0.0
    for index in range(stop):
        block = slice(index * atoms_per_bin, (index + 1) * atoms_per_bin)
        caught = selected_captured[block]
        weight = float(selected_weights[block].sum())
        p = caught.mean()
        variance += weight**2 * p * (1-p) / atoms_per_bin
        low, high = wilson_interval(int(caught.sum()), atoms_per_bin, confidence_level)
        confidence_low += weight * low
        confidence_high += weight * high
    omitted_probability = float(1.0 - selected_weights.sum())
    return CaptureEstimate(
        selected_captured,selected_times,selected_speeds,
        probability,float(np.sqrt(variance)),
        estimate.incident_flux_s,float(estimate.incident_flux_s*probability),estimate.seed,
        selected_weights,omitted_probability,
        (float(confidence_low),float(min(1.0, confidence_high + omitted_capture_upper))),
        float(edges[stop]),float(omitted_capture_upper),
        float(estimate.incident_flux_s*omitted_capture_upper),
        converged,
    )


def capture_response_map(force_model, criterion, speeds_m_s, impact_edges_m,
                         *, capture_surface_radius_m, samples_per_cell,
                         max_step_s, seed, rtol=1e-7, atol=1e-10):
    """Temperature-independent P_capture(speed, impact parameter) response."""
    speeds=np.asarray(speeds_m_s,float); impact_edges=np.asarray(impact_edges_m,float)
    rng=np.random.default_rng(seed)
    probability=np.empty((len(speeds),len(impact_edges)-1))
    low=np.empty_like(probability); high=np.empty_like(probability)
    for i,speed in enumerate(speeds):
        for j,(b0,b1) in enumerate(zip(impact_edges[:-1],impact_edges[1:])):
            # b=R*sin(theta); cosine-law flux makes b² uniform on [0,R²].
            b=np.sqrt(rng.uniform(b0*b0,b1*b1,samples_per_cell))
            mu=np.sqrt(1-(b/capture_surface_radius_m)**2)
            phi=2*np.pi*rng.random(samples_per_cell)
            positions=np.tile([-capture_surface_radius_m,0,0],(samples_per_cell,1))
            velocities=speed*np.column_stack([mu,np.sqrt(1-mu*mu)*np.cos(phi),
                                               np.sqrt(1-mu*mu)*np.sin(phi)])
            caught,_=evaluate_capture(force_model,positions,velocities,criterion,
                                      max_step=max_step_s,rtol=rtol,atol=atol)
            probability[i,j]=caught.mean()
            low[i,j],high[i,j]=wilson_interval(int(caught.sum()),samples_per_cell)
    # Impact area/flux weights are proportional to delta(b²).
    weights=np.diff(impact_edges**2)/capture_surface_radius_m**2
    return probability,low,high,weights
