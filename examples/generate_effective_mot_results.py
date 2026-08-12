"""Generate reproducible Effective-model MOT data and figures from documented models."""

import json
import copy
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model, load_config
from cold_atom_mot.magnetic.coils import AntiHelmholtzPair
from cold_atom_mot.solvers.deterministic import integrate_trajectory
from cold_atom_mot.solvers.monte_carlo import simulate_photon_events

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "effective_mot"
CONFIG_PATH = ROOT / "configs" / "rb87_d2_mot.yaml"
COLORS = {"field": "#315a8c", "force": "#a33b20", "deterministic": "#315a8c", "monte_carlo": "#3a8f65"}


def save(name: str) -> None:
    """Save PNG and SVG from the same current figure."""
    plt.tight_layout()
    plt.savefig(OUTPUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    config = load_config(CONFIG_PATH)
    model = build_effective_model(config)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-effective"

    # Apparatus geometry: independent beam axes and the configured physical pair.
    pair = AntiHelmholtzPair.symmetric(0.04, 0.04, 2.0, 50, segments=256)
    theta = np.linspace(0, 2 * np.pi, 200)
    fig = plt.figure(figsize=(7.5, 6))
    axis = fig.add_subplot(projection="3d")
    for beam in model.beams:
        line = np.vstack([-0.045 * beam.direction, 0.045 * beam.direction])
        axis.plot(*line.T, linewidth=3, alpha=0.75, label=beam.label)
        axis.quiver(*(-0.035 * beam.direction), *(0.012 * beam.direction), arrow_length_ratio=0.3)
    for z in (-0.02, 0.02):
        axis.plot(0.04 * np.cos(theta), 0.04 * np.sin(theta), np.full_like(theta, z), color="black", linewidth=2)
    axis.set(xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", title="Six independent MOT beams and anti-Helmholtz coils")
    axis.set_box_aspect((1, 1, 1))
    save("apparatus_geometry")

    # Physical field map in the x-z plane.
    coordinate = np.linspace(-0.018, 0.018, 45)
    x_grid, z_grid = np.meshgrid(coordinate, coordinate)
    field_points = np.stack([x_grid, np.zeros_like(x_grid), z_grid], axis=-1)
    field = pair.field(field_points)
    field_magnitude = np.linalg.norm(field, axis=-1)
    plt.figure(figsize=(7, 5.5))
    image = plt.pcolormesh(x_grid * 1e3, z_grid * 1e3, field_magnitude * 1e3, shading="auto", cmap="viridis")
    stride = 4
    plt.quiver(x_grid[::stride, ::stride] * 1e3, z_grid[::stride, ::stride] * 1e3, field[::stride, ::stride, 0], field[::stride, ::stride, 2], color="white", alpha=0.8)
    plt.colorbar(image, label="|B| (mT)")
    plt.xlabel("x (mm)"); plt.ylabel("z (mm)"); plt.title("Segmented Biot–Savart anti-Helmholtz field (y=0)")
    save("antihelmholtz_field")

    # Force-position/velocity slices and the genuine 2D F(x,vx) scan.
    positions = np.linspace(-0.004, 0.004, 121)
    velocities = np.linspace(-8, 8, 121)
    force_map = np.empty((len(velocities), len(positions)))
    for i, velocity in enumerate(velocities):
        r = np.column_stack([positions, np.zeros_like(positions), np.zeros_like(positions)])
        v = np.tile([velocity, 0, 0], (len(positions), 1))
        force_map[i] = model.force(r, v)[:, 0]
    plt.figure(figsize=(7, 5.5))
    limit = np.max(abs(force_map))
    image = plt.pcolormesh(positions * 1e3, velocities, force_map / 1e-20, shading="auto", cmap="coolwarm", vmin=-limit / 1e-20, vmax=limit / 1e-20)
    plt.colorbar(image, label="$F_x$ ($10^{-20}$ N)")
    plt.contour(positions * 1e3, velocities, force_map, levels=[0], colors="black", linewidths=1)
    plt.xlabel("x (mm)"); plt.ylabel("$v_x$ (m/s)"); plt.title("effective-model deterministic force map")
    save("force_map_x_vx")

    # Three-dimensional deterministic trajectories, without claiming capture.
    initial_conditions = [([0.002, 0.001, 0.001], [0, 0, 0]), ([-0.002, 0.001, -0.001], [0.3, -0.1, 0]), ([0.001, -0.002, 0.0015], [-0.2, 0.2, -0.1])]
    trajectories = [integrate_trajectory(model, position, velocity, 0.004, max_step=1e-5) for position, velocity in initial_conditions]
    fig = plt.figure(figsize=(7.5, 6))
    axis = fig.add_subplot(projection="3d")
    for index, trajectory in enumerate(trajectories, 1):
        axis.plot(*(trajectory.position * 1e3).T, label=f"trajectory {index}")
        axis.scatter(*(trajectory.position[0] * 1e3), marker="o")
    axis.scatter(0, 0, 0, marker="+", color="black", s=100, label="field zero")
    axis.set(xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)", title="Adaptive mean-force 3D trajectories")
    axis.legend()
    save("deterministic_trajectories")

    # Photon-event ensemble and atom-number convergence of the sample mean.
    sizes = np.array([32, 64, 128, 256, 512])
    mean_velocity = []
    standard_error = []
    final_cloud = None
    for count in sizes:
        initial_position = np.zeros((count, 3))
        initial_velocity = np.tile([0.05, 0, 0], (count, 1))
        run = simulate_photon_events(model, initial_position, initial_velocity, 2e-6, 5e-9, seed=20260811, store_every=400)
        final = run.velocity[-1, :, 0]
        mean_velocity.append(final.mean())
        standard_error.append(final.std(ddof=1) / np.sqrt(count))
        if count == sizes[-1]:
            final_cloud = run.velocity[-1]
    plt.figure(figsize=(7, 4.8))
    plt.errorbar(sizes, mean_velocity, yerr=standard_error, marker="o", capsize=4, color=COLORS["monte_carlo"])
    plt.xscale("log", base=2)
    plt.xlabel("Monte Carlo atoms"); plt.ylabel("Final mean $v_x$ (m/s)")
    plt.title("Photon-event ensemble convergence (2 µs, fixed seed)")
    save("monte_carlo_convergence")

    plt.figure(figsize=(7, 4.8))
    plt.hist(final_cloud[:, 0], bins=25, alpha=0.75, label="$v_x$")
    plt.hist(final_cloud[:, 1], bins=25, alpha=0.55, label="$v_y$")
    plt.xlabel("Velocity (m/s)"); plt.ylabel("Atoms per bin")
    plt.title("Photon-recoil velocity distribution (512 atoms, 2 µs)")
    plt.legend()
    save("recoil_velocity_distribution")

    # Sensitivities are recomputed from the force/coil models, not scaling laws.
    detunings = np.linspace(-4.0, -0.6, 24)
    powers = np.linspace(0.002, 0.020, 24)
    waists = np.linspace(0.004, 0.014, 24)
    gradients = np.linspace(0.03, 0.20, 24)
    damping_detuning = []
    damping_power = []
    off_axis_force = []
    restoring_gradient = []
    for value in detunings:
        varied = copy.deepcopy(config)
        varied["laser"]["detuning_gamma"] = float(value)
        damping_detuning.append(build_effective_model(varied).linear_coefficients()[0][0])
    for value in powers:
        varied = copy.deepcopy(config)
        varied["laser"]["power_per_beam_w"] = float(value)
        damping_power.append(build_effective_model(varied).linear_coefficients()[0][0])
    for value in waists:
        varied = copy.deepcopy(config)
        varied["laser"]["waist_m"] = float(value)
        varied_model = build_effective_model(varied)
        off_axis_force.append(varied_model.force([0.002, 0, 0], [0, 0, 0])[0])
    for value in gradients:
        varied = copy.deepcopy(config)
        varied["magnetic_field"]["radial_gradient_t_per_m"] = float(value)
        restoring_gradient.append(build_effective_model(varied).linear_coefficients()[1][0])

    bias = np.linspace(-50e-6, 50e-6, 25)
    ideal_zero_x = -bias / config["magnetic_field"]["radial_gradient_t_per_m"]
    tilt_degrees = np.array([0.0, 0.25, 0.5, 1.0, 2.0])
    tilted_zero = np.array([
        AntiHelmholtzPair.symmetric(
            0.04,
            0.04,
            2.0,
            50,
            segments=128,
            tilt_y=np.deg2rad(angle),
            lateral_offset=5e-4,
        ).field_zero()
        for angle in tilt_degrees
    ])

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5))
    axes[0, 0].plot(detunings, np.array(damping_detuning) * 1e22, color=COLORS["force"])
    axes[0, 0].set(xlabel="Detuning (Γ)", ylabel="βx ($10^{-22}$ kg/s)", title="Damping vs detuning")
    axes[0, 1].plot(powers * 1e3, np.array(damping_power) * 1e22, color=COLORS["force"])
    axes[0, 1].set(xlabel="Power per beam (mW)", ylabel="βx ($10^{-22}$ kg/s)", title="Damping vs power")
    axes[0, 2].plot(waists * 1e3, np.array(off_axis_force) * 1e20, color=COLORS["force"])
    axes[0, 2].set(xlabel="Beam waist (mm)", ylabel="$F_x$(x=2 mm) ($10^{-20}$ N)", title="Waist dependence")
    axes[1, 0].plot(gradients, np.array(restoring_gradient) * 1e18, color=COLORS["field"])
    axes[1, 0].set(xlabel="Radial gradient (T/m)", ylabel="κx ($10^{-18}$ N/m)", title="Restoring coefficient")
    axes[1, 1].plot(bias * 1e6, ideal_zero_x * 1e3, color=COLORS["field"])
    axes[1, 1].set(xlabel="Uniform $B_x$ (µT)", ylabel="Field-zero x (mm)", title="Stray-field displacement")
    axes[1, 2].plot(tilt_degrees, np.linalg.norm(tilted_zero, axis=1) * 1e3, marker="o", color=COLORS["field"])
    axes[1, 2].set(xlabel="Second-coil tilt (degree)", ylabel="|field-zero shift| (mm)", title="Tilt with 0.5 mm offset")
    fig.suptitle("Effective-model parameter sensitivities (model calculations, not fitted scalings)")
    save("parameter_sensitivities")

    metadata = {
        "configuration": config,
        "random_seed": 20260811,
        "simulation_version": __version__,
        "units": {"position": "m", "velocity": "m/s", "force": "N", "magnetic_field": "T", "time": "s"},
        "solver": "RK45 mean force; Bernoulli photon events with isotropic spontaneous recoil",
        "model_fidelity": "Effective model effective two-level semiclassical/stochastic",
        "monte_carlo_time_step_s": 5e-9,
        "maximum_monte_carlo_atoms": 512,
        "coil_segments_per_loop": 256,
    }
    np.savez_compressed(
        OUTPUT / "effective_mot_reference.npz",
        x_grid_m=x_grid,
        z_grid_m=z_grid,
        magnetic_field_t=field,
        force_positions_m=positions,
        force_velocities_m_per_s=velocities,
        force_x_n=force_map,
        convergence_atoms=sizes,
        convergence_mean_velocity_m_per_s=np.array(mean_velocity),
        convergence_standard_error_m_per_s=np.array(standard_error),
        final_velocity_m_per_s=final_cloud,
        detuning_gamma=detunings,
        damping_vs_detuning_kg_per_s=np.array(damping_detuning),
        power_per_beam_w=powers,
        damping_vs_power_kg_per_s=np.array(damping_power),
        beam_waist_m=waists,
        off_axis_force_n=np.array(off_axis_force),
        radial_gradient_t_per_m=gradients,
        restoring_vs_gradient_n_per_m=np.array(restoring_gradient),
        uniform_bias_x_t=bias,
        ideal_field_zero_x_m=ideal_zero_x,
        second_coil_tilt_degree=tilt_degrees,
        tilted_field_zero_m=tilted_zero,
        metadata_json=json.dumps(metadata),
    )


if __name__ == "__main__":
    main()
