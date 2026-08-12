"""Command-line entry points for reproducible fidelity-labelled runs."""

import argparse
import json
from pathlib import Path
import numpy as np
from . import __version__
from .io.config import build_effective_model, build_multilevel_model, build_subdoppler_model, build_vapor_state, load_config
from .solvers.deterministic import integrate_trajectory
from .solvers.monte_carlo import simulate_photon_events
from .atomic.species import get_atomic_line
from .physics.optical_bloch import TwoLevelOBE
from .simulation.capture import CaptureCriterion, estimate_adaptive_vapor_capture_rate
from .vacuum import (background_collision_loss_rate_s,
                     gaussian_two_body_effective_volume, loading_curve)


def simulate(config_path: str) -> None:
    config = load_config(config_path)
    model = build_effective_model(config)
    simulation = config["simulation"]
    output = Path(config["output"]["directory"])
    output.mkdir(parents=True, exist_ok=True)
    trajectory = integrate_trajectory(model, simulation["initial_position_m"], simulation["initial_velocity_m_per_s"], simulation["duration_s"], max_step=simulation["max_step_s"])
    mc = config["monte_carlo"]
    positions = np.tile(simulation["initial_position_m"], (mc["atoms"], 1))
    velocities = np.tile(simulation["initial_velocity_m_per_s"], (mc["atoms"], 1))
    stochastic = simulate_photon_events(model, positions, velocities, simulation["duration_s"], mc["time_step_s"], seed=mc["seed"], store_every=mc["store_every"])
    metadata = {"simulation_version": __version__, "config": config, "units": {"position": "m", "velocity": "m/s", "time": "s"}, "solver": "RK45 effective two-level + discrete photon Monte Carlo", "model_fidelity": "Effective model effective two-level semiclassical/stochastic", "random_seed": mc["seed"], "number_of_atoms": mc["atoms"], "time_step_s": mc["time_step_s"]}
    np.savez_compressed(output / "effective_mot_run.npz", deterministic_time=trajectory.time, deterministic_position=trajectory.position, deterministic_velocity=trajectory.velocity, monte_carlo_time=stochastic.time, monte_carlo_position=stochastic.position, monte_carlo_velocity=stochastic.velocity, metadata_json=json.dumps(metadata))
    print(output / "effective_mot_run.npz")


def rate_equation(config_path: str) -> None:
    """Evaluate reproducible rate-equation steady populations and a force profile."""
    config = load_config(config_path)
    model = build_multilevel_model(config)
    output = Path(config["output"]["directory"])
    output.mkdir(parents=True, exist_ok=True)
    positions = np.linspace(-config["rate_equation"]["position_extent_m"], config["rate_equation"]["position_extent_m"], config["rate_equation"]["points"])
    force = np.empty_like(positions)
    manifolds = np.empty((len(positions), 3))
    for index, x in enumerate(positions):
        population = model.steady_state(np.array([x, 0, 0]), np.zeros(3))
        force[index] = model.force(np.array([x, 0, 0]), np.zeros(3), population)[0]
        values = model.manifold_populations(population)
        manifolds[index] = [values["ground_F1"], values["ground_F2"], values["excited"]]
    metadata = {
        "simulation_version": __version__,
        "config": config,
        "units": {"position": "m", "force": "N", "population": "dimensionless"},
        "solver": "stationary linear hyperfine/Zeeman population rate equations",
        "model_fidelity": "Rate-equation model multilevel rate equation",
    }
    path = output / "multilevel_rate_equation.npz"
    np.savez_compressed(path, position_m=positions, force_x_n=force, manifold_population=manifolds, metadata_json=json.dumps(metadata))
    print(path)



def optical_bloch(config_path: str) -> None:
    """Run the explicitly reduced Two-level OBE two-level OBE configuration."""
    config = load_config(config_path)
    line = get_atomic_line(config["atom"]["isotope"], config["atom"]["line"])
    parameters = config["obe"]
    model = TwoLevelOBE.from_saturation(
        line.gamma_rad_s,
        parameters["detuning_gamma"] * line.gamma_rad_s,
        parameters["saturation"],
        dephasing_rate=parameters.get("pure_dephasing_gamma", 0.0) * line.gamma_rad_s,
    )
    duration = parameters["duration_lifetimes"] / line.gamma_rad_s
    time, density = model.evolve(
        np.array([[1, 0], [0, 0]], dtype=complex),
        duration,
        rtol=parameters["rtol"],
        atol=parameters["atol"],
        max_step=parameters["max_step_lifetimes"] / line.gamma_rad_s,
    )
    output = Path(config["output"]["directory"])
    output.mkdir(parents=True, exist_ok=True)
    metadata = {
        "simulation_version": __version__,
        "config": config,
        "units": {"time": "s", "density_matrix": "dimensionless"},
        "solver": "adaptive two-level Lindblad optical Bloch equation",
        "model_fidelity": "Coherent model reduced single-transition OBE",
    }
    path = output / "two_level_obe_run.npz"
    np.savez_compressed(path, time_s=time, density_matrix=density, steady_density_matrix=model.steady_state(), metadata_json=json.dumps(metadata))
    print(path)


def subdoppler(config_path: str) -> None:
    """Evaluate a reproducible polarization-gradient cycle-averaged force point."""
    config = load_config(config_path); model = build_subdoppler_model(config)
    parameters = config["simulation"]; velocity = parameters["velocity_m_per_s"]
    options = dict(periods=parameters["periods"], discard=parameters["discard_periods"],
                   steps_per_period=parameters["steps_per_period"])
    force = model.moving_average_force(velocity, **options)
    friction = model.friction_coefficient(velocity, **options)
    output = Path(config["output"]["directory"]); output.mkdir(parents=True, exist_ok=True)
    metadata = {"simulation_version": __version__, "config": config,
                "units": {"velocity": "m/s", "force": "N", "friction": "kg/s"},
                "solver": "adiabatically eliminated F=2 to Fprime=3 optical-pumping trajectory",
                "model_fidelity": "phase-resolved adiabatic population Sisyphus model"}
    path = output / "polarization_gradient_run.npz"
    np.savez_compressed(path, velocity_m_per_s=velocity, force_x_n=force,
                        friction_kg_per_s=friction, metadata_json=json.dumps(metadata))
    print(path)


def vapor_loading(config_path: str) -> None:
    """Estimate trajectory-linked vapour loading and apply configured losses."""
    config = load_config(config_path)
    mot_path = Path(config_path).resolve().parents[1] / config["mot_config"]
    force_model = build_effective_model(load_config(mot_path))
    vapor = build_vapor_state(config)
    capture = config["capture"]
    criterion = CaptureCriterion(**capture["criterion"])
    estimate = estimate_adaptive_vapor_capture_rate(
        force_model,
        vapor,
        criterion,
        capture_surface_radius_m=capture["surface_radius_m"],
        initial_speed_edges_m_s=capture["speed_bin_edges_m_s"],
        atoms_per_bin=capture["atoms_per_bin"],
        maximum_speed_m_s=capture["maximum_speed_m_s"],
        tail_relative_loading_tolerance=capture["tail_relative_loading_tolerance"],
        max_step_s=capture["max_step_s"],
        seed=capture["seed"],
        confidence_level=capture["confidence_level"],
        rtol=capture["rtol"],atol=capture["atol"],
    )
    loss = config["loading"]
    background_loss = loss["background_one_body_loss_s"]
    collision = loss.get("background_collision_model")
    components = loss.get("background_gas_components")
    if collision is not None and components:
        raise ValueError("choose aggregate background collision model or components, not both")
    if collision is not None:
        background_loss = background_collision_loss_rate_s(
            vapor.background_gas_pressure_pa,
            vapor.background_temperature_k,
            force_model.atom.mass_kg,
            collision["particle_mass_kg"],
            collision["effective_loss_cross_section_m2"],
        )
    elif components:
        background_loss=0.0
        for component in components:
            background_loss += background_collision_loss_rate_s(
                component["partial_pressure_pa"],vapor.background_temperature_k,
                force_model.atom.mass_kg,component["particle_mass_kg"],
                component["effective_loss_cross_section_m2"],
            )
    one_body_loss = background_loss + loss["hot_rb_one_body_loss_s"]
    time = np.linspace(0, loss["curve_duration_s"], loss["curve_points"])
    effective_volume=loss["effective_volume_m3"]
    if loss.get("gaussian_cloud_sigma_m") is not None:
        effective_volume=gaussian_two_body_effective_volume(**loss["gaussian_cloud_sigma_m"])
    atom_number = loading_curve(
        time,
        estimate.loading_rate_s,
        one_body_loss,
        two_body_coefficient=loss["two_body_loss_m3_s"],
        effective_volume_m3=effective_volume,
    )
    metadata = {
        "simulation_version": __version__,
        "config": config,
        "units": {"speed": "m/s", "time": "s", "loading_rate": "atoms/s"},
        "solver": "stratified thermal-flux deterministic capture plus loading ODE",
        "limitations": "capture surface is an acceptance boundary; loss inputs require calibration",
    }
    output = Path(config["output"]["directory"])
    output.mkdir(parents=True, exist_ok=True)
    path = output / "vapor_loading_run.npz"
    np.savez_compressed(
        path,
        initial_speed_m_s=estimate.initial_speed_m_s,
        captured=estimate.captured,
        sample_weight=estimate.sample_weights,
        capture_time_s=estimate.capture_time_s,
        capture_probability=estimate.capture_probability,
        capture_probability_standard_error=estimate.capture_probability_standard_error,
        incident_flux_s=estimate.incident_flux_s,
        loading_rate_s=estimate.loading_rate_s,
        omitted_high_speed_probability=estimate.omitted_high_speed_probability,
        capture_probability_confidence_interval=estimate.confidence_interval,
        loading_rate_confidence_interval_s=estimate.loading_rate_confidence_interval_s,
        last_simulated_speed_m_s=estimate.last_simulated_speed_m_s,
        omitted_capture_probability_upper=estimate.omitted_capture_probability_upper,
        omitted_loading_rate_upper_s=estimate.omitted_loading_rate_upper_s,
        tail_converged=estimate.tail_converged,
        time_s=time,
        atom_number=atom_number,
        rb_partial_pressure_pa=vapor.rb_partial_pressure_pa,
        background_gas_pressure_pa=vapor.background_gas_pressure_pa,
        background_one_body_loss_s=background_loss,
        total_one_body_loss_s=one_body_loss,
        metadata_json=json.dumps(metadata, sort_keys=True),
    )
    print(path)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="cold-atom-mot")
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("simulate")
    command.add_argument("config")
    rate_parser = subparsers.add_parser("rate-equation")
    rate_parser.add_argument("config")
    obe_parser = subparsers.add_parser("obe")
    obe_parser.add_argument("config")
    pgc_parser = subparsers.add_parser("subdoppler")
    pgc_parser.add_argument("config")
    loading_parser = subparsers.add_parser("loading")
    loading_parser.add_argument("config")
    args = parser.parse_args(argv)
    if args.command == "simulate":
        simulate(args.config)
    elif args.command == "rate-equation":
        rate_equation(args.config)
    elif args.command == "obe":
        optical_bloch(args.config)
    elif args.command == "subdoppler":
        subdoppler(args.config)
    elif args.command == "loading":
        vapor_loading(args.config)
