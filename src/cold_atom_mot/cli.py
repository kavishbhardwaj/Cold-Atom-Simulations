"""Command-line entry points for reproducible Phase-1 runs."""

import argparse
import json
from pathlib import Path
import numpy as np
from . import __version__
from .io.config import build_effective_model, build_multilevel_model, load_config
from .solvers.deterministic import integrate_trajectory
from .solvers.monte_carlo import simulate_photon_events


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
    metadata = {"simulation_version": __version__, "config": config, "units": {"position": "m", "velocity": "m/s", "time": "s"}, "solver": "RK45 effective two-level + discrete photon Monte Carlo", "model_fidelity": "Level A effective two-level semiclassical/stochastic", "random_seed": mc["seed"], "number_of_atoms": mc["atoms"], "time_step_s": mc["time_step_s"]}
    np.savez_compressed(output / "phase1_run.npz", deterministic_time=trajectory.time, deterministic_position=trajectory.position, deterministic_velocity=trajectory.velocity, monte_carlo_time=stochastic.time, monte_carlo_position=stochastic.position, monte_carlo_velocity=stochastic.velocity, metadata_json=json.dumps(metadata))
    print(output / "phase1_run.npz")


def rate_equation(config_path: str) -> None:
    """Evaluate reproducible Level-B steady populations and a force profile."""
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
        "model_fidelity": "Level B multilevel rate equation",
    }
    path = output / "phase2_rate_equation.npz"
    np.savez_compressed(path, position_m=positions, force_x_n=force, manifold_population=manifolds, metadata_json=json.dumps(metadata))
    print(path)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="cold-atom-mot")
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("simulate")
    command.add_argument("config")
    level_b = subparsers.add_parser("rate-equation")
    level_b.add_argument("config")
    args = parser.parse_args(argv)
    if args.command == "simulate":
        simulate(args.config)
    elif args.command == "rate-equation":
        rate_equation(args.config)
