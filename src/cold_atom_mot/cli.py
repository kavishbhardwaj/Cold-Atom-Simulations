"""Command-line entry points for reproducible Phase-1 runs."""

import argparse
import json
from pathlib import Path
import numpy as np
from . import __version__
from .io.config import build_effective_model, load_config
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


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="cold-atom-mot")
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("simulate")
    command.add_argument("config")
    args = parser.parse_args(argv)
    if args.command == "simulate":
        simulate(args.config)
