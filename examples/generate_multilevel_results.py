"""Generate reproducible rate-equation multilevel rate-equation diagnostics."""

import copy
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model, build_multilevel_model, load_config

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "rb87_d2_multilevel.yaml"
OUTPUT = ROOT / "results" / "multilevel"


def save(name: str) -> None:
    plt.tight_layout()
    plt.savefig(OUTPUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    config = load_config(CONFIG_PATH)
    model = build_multilevel_model(config)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-multilevel"

    position = np.linspace(-3e-3, 3e-3, 81)
    force_multilevel = np.empty_like(position)
    manifold = np.empty((len(position), 3))
    for index, x in enumerate(position):
        point = np.array([x, 0.0, 0.0])
        population = model.steady_state(point, np.zeros(3))
        force_multilevel[index] = model.force(point, np.zeros(3), population)[0]
        values = model.manifold_populations(population)
        manifold[index] = [values["ground_F1"], values["ground_F2"], values["excited"]]

    level_a_config = copy.deepcopy(config)
    level_a_config["model"] = "effective_mot"
    level_a = build_effective_model(level_a_config)
    force_effective = np.array([level_a.force([x, 0, 0], [0, 0, 0])[0] for x in position])
    plt.figure(figsize=(7.2, 4.8))
    plt.plot(position * 1e3, force_effective * 1e21, label="Effective model: effective two-level", linewidth=2)
    plt.plot(position * 1e3, force_multilevel * 1e21, label="Rate-equation model: 24-state rate equation", linewidth=2)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xlabel("x (mm)"); plt.ylabel("$F_x$ ($10^{-21}$ N)")
    plt.title("Restoring-force fidelity comparison")
    plt.legend()
    save("effective_vs_multilevel_force")

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(position * 1e3, manifold[:, 0], label="5S$_{1/2}$ F=1")
    plt.plot(position * 1e3, manifold[:, 1], label="5S$_{1/2}$ F=2")
    plt.plot(position * 1e3, manifold[:, 2], label="5P$_{3/2}$ total")
    plt.xlabel("x (mm)"); plt.ylabel("Steady-state population")
    plt.title("Hyperfine-manifold populations across the MOT")
    plt.legend()
    save("manifold_populations")

    velocity = np.linspace(-1.0, 1.0, 81)
    force_velocity = np.array([model.force(np.zeros(3), [v, 0, 0])[0] for v in velocity])
    plt.figure(figsize=(7.2, 4.8))
    plt.plot(velocity, force_velocity * 1e21, linewidth=2)
    plt.axhline(0, color="black", linewidth=0.8)
    plt.xlabel("$v_x$ (m/s)"); plt.ylabel("$F_x$ ($10^{-21}$ N)")
    plt.title("rate-equation damping force near zero velocity")
    save("multilevel_force_velocity")

    repump_power = np.geomspace(1e-7, 2e-3, 40)
    f2_population = []
    excited_population = []
    for power in repump_power:
        varied = copy.deepcopy(config)
        varied["repump"]["power_per_beam_w"] = float(power)
        varied_model = build_multilevel_model(varied)
        population = varied_model.steady_state(np.zeros(3), np.zeros(3))
        values = varied_model.manifold_populations(population)
        f2_population.append(values["ground_F2"])
        excited_population.append(values["excited"])
    plt.figure(figsize=(7.2, 4.8))
    plt.semilogx(repump_power * 1e3, f2_population, label="ground F=2")
    plt.semilogx(repump_power * 1e3, excited_population, label="excited total")
    plt.axvline(config["repump"]["power_per_beam_w"] * 1e3, color="black", linestyle="--", label="reference")
    plt.xlabel("Repump power per beam (mW)"); plt.ylabel("Steady-state population")
    plt.title("Optical-pumping response to repump power")
    plt.legend()
    save("repump_power_scan")

    centre = model.steady_state(np.zeros(3), np.zeros(3))
    ground_labels = [f"F={state.F}, m={state.m:+d}" for state in model.basis.ground]
    plt.figure(figsize=(9, 4.8))
    plt.bar(np.arange(len(ground_labels)), centre[:len(ground_labels)])
    plt.xticks(np.arange(len(ground_labels)), ground_labels, rotation=45, ha="right")
    plt.ylabel("Steady-state population")
    plt.title("Ground-state Zeeman populations at the field zero")
    save("zeeman_populations")

    metadata = {
        "configuration": config,
        "simulation_version": __version__,
        "units": {"position": "m", "velocity": "m/s", "force": "N", "population": "dimensionless", "power": "W"},
        "solver": "stationary 24-state hyperfine/Zeeman population rate equations",
        "model_fidelity": "Rate-equation model multilevel rate equation",
        "basis": "87Rb 5S1/2 F=1,2 and 5P3/2 F'=0,1,2,3 Zeeman states",
        "coherences_included": False,
    }
    np.savez_compressed(
        OUTPUT / "multilevel_reference.npz",
        position_m=position,
        force_effective_n=force_effective,
        force_multilevel_n=force_multilevel,
        manifold_population=manifold,
        velocity_m_per_s=velocity,
        force_velocity_n=force_velocity,
        repump_power_per_beam_w=repump_power,
        f2_population=np.array(f2_population),
        excited_population=np.array(excited_population),
        centre_ground_population=centre[:len(model.basis.ground)],
        ground_labels=np.array(ground_labels),
        metadata_json=json.dumps(metadata),
    )


if __name__ == "__main__":
    main()
