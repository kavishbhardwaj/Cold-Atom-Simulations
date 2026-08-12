"""Command-line entry points for reproducible fidelity-labelled runs."""

import argparse
import json
from pathlib import Path
import numpy as np
from . import __version__
from .io.config import build_effective_model, build_multilevel_model, build_subdoppler_model, load_config
from .solvers.deterministic import integrate_trajectory
from .solvers.monte_carlo import simulate_photon_events
from .atomic.species import get_atomic_line
from .physics.optical_bloch import TwoLevelOBE


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
                "model_fidelity": "Polarization-gradient model phase-resolved population Sisyphus model"}
    path = output / "polarization_gradient_run.npz"
    np.savez_compressed(path, velocity_m_per_s=velocity, force_x_n=force,
                        friction_kg_per_s=friction, metadata_json=json.dumps(metadata))
    print(path)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="cold-atom-mot")
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("simulate")
    command.add_argument("config")
    level_b = subparsers.add_parser("rate-equation")
    level_b.add_argument("config")
    level_c = subparsers.add_parser("obe")
    level_c.add_argument("config")
    level_d = subparsers.add_parser("subdoppler")
    level_d.add_argument("config")
    args = parser.parse_args(argv)
    if args.command == "simulate":
        simulate(args.config)
    elif args.command == "rate-equation":
        rate_equation(args.config)
    elif args.command == "obe":
        optical_bloch(args.config)
    elif args.command == "subdoppler":
        subdoppler(args.config)
