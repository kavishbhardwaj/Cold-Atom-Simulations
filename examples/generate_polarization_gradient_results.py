"""Generate reproducible polarization-gradient polarization-gradient diagnostics."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar, k as k_B

from cold_atom_mot import __version__
from cold_atom_mot.atomic.species import build_atomic_basis, get_atomic_line
from cold_atom_mot.physics.subdoppler import coherent_six_beam_field, PolarizationGradientModel

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "polarization_gradient"
ATOM = get_atomic_line("87Rb", "D2")
PHASES = np.array([0, 0, 0, np.pi / 2, 0, np.pi / 4])


def build(saturation=0.08, detuning_gamma=-3, bias_t=0):
    beams = coherent_six_beam_field(ATOM.wave_number, saturation, PHASES)
    return PolarizationGradientModel(build_atomic_basis("87Rb", "D2"), 2, 3,
                                     detuning_gamma * ATOM.gamma, beams, magnetic_field_t=[0, 0, bias_t])


def save(name):
    plt.tight_layout(); plt.savefig(OUTPUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None}); plt.close()


def doppler_force(velocity, saturation=0.08, detuning_gamma=-3):
    k, gamma = ATOM.wave_number, ATOM.gamma
    plus = detuning_gamma * gamma - k * velocity
    minus = detuning_gamma * gamma + k * velocity
    return hbar * k * gamma / 2 * saturation * (
        1 / (1 + 2 * saturation + (2 * plus / gamma) ** 2) -
        1 / (1 + 2 * saturation + (2 * minus / gamma) ** 2))


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True); plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-pgc"
    model = build(); period = 2 * np.pi / ATOM.wave_number
    x = np.linspace(0, period, 161)
    fractions = {q: [] for q in (-1, 0, 1)}; shifts = []; populations = []
    for value in x:
        _, local = model.polarization_components([value, 0, 0])
        for q in fractions: fractions[q].append(local[q])
        shifts.append(model.light_shifts([value, 0, 0])); populations.append(model.stationary_populations([value, 0, 0]))
    shifts, populations = np.asarray(shifts), np.asarray(populations)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for q, label in ((-1, "$\\sigma^-$"), (0, "$\\pi$"), (1, "$\\sigma^+$")):
        axes[0].plot(x / period, fractions[q], label=label)
        axes[0].set(xlabel="Position x / λ", ylabel="Local polarization fraction", title="Phase-resolved spherical components"); axes[0].legend()
    intensity = [np.vdot(model.electric_field([value, 0, 0]), model.electric_field([value, 0, 0])).real for value in x]
    axes[1].plot(x / period, intensity)
    axes[1].set(xlabel="Position x / λ", ylabel="Coherent saturation sum", title="Six-beam interference intensity")
    save("polarization_lattice")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for index, m in enumerate(model.ground_m):
        axes[0].plot(x / period, shifts[:, index] / k_B * 1e6, label=f"m={m:+d}")
        axes[1].plot(x / period, populations[:, index], label=f"m={m:+d}")
    axes[0].set(xlabel="Position x / λ", ylabel="Light shift U/kB (µK)", title="Ground-state light-shift potentials")
    axes[1].set(xlabel="Position x / λ", ylabel="Stationary population", title="Position-dependent optical pumping")
    axes[0].legend(ncol=2); axes[1].legend(ncol=2); save("light_shifts_pumping")

    velocities = np.r_[-np.geomspace(0.004, 0.12, 10)[::-1], np.geomspace(0.004, 0.12, 10)]
    pg_force = np.array([model.moving_average_force(v, periods=12, discard=6, steps_per_period=40) for v in velocities])
    d_force = np.array([doppler_force(v) for v in velocities])
    plt.figure(figsize=(7.5, 4.8)); plt.plot(velocities, pg_force * 1e22, "o-", label="Polarization-gradient model dipole/pumping")
    plt.plot(velocities, d_force * 1e22, label="two-level Doppler pair")
    plt.axhline(0, color="black", linewidth=.7); plt.xlabel("Velocity vx (m/s)"); plt.ylabel("Cycle-averaged force Fx ($10^{-22}$ N)")
    plt.title("Sub-Doppler structure versus Doppler radiation pressure"); plt.legend(); save("subdoppler_force_velocity")

    probe_velocity = 0.02
    biases = np.array([0, 2, 5, 10, 20, 50]) * 1e-6
    bias_beta = np.array([build(bias_t=b).friction_coefficient(probe_velocity, periods=10, discard=5, steps_per_period=32) for b in biases])
    detunings = np.array([-1.5, -2, -3, -4, -5])
    saturations = np.array([0.02, 0.04, 0.08, 0.16, 0.30])
    beta_grid = np.empty((len(saturations), len(detunings)))
    for i, saturation in enumerate(saturations):
        for j, detuning in enumerate(detunings):
            beta_grid[i, j] = build(saturation, detuning).friction_coefficient(probe_velocity, periods=8, discard=4, steps_per_period=28)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    axes[0].plot(biases * 1e6, bias_beta * 1e23, "o-"); axes[0].set(xlabel="Axial residual field (µT)", ylabel="Friction β ($10^{-23}$ kg/s)", title="Residual-field sensitivity")
    image = axes[1].pcolormesh(detunings, saturations, beta_grid * 1e23, shading="nearest", cmap="coolwarm")
    fig.colorbar(image, ax=axes[1], label="β ($10^{-23}$ kg/s)"); axes[1].set(xlabel="Detuning (Γ)", ylabel="Saturation per beam", title="Intensity–detuning dependence")
    save("subdoppler_sensitivities")

    steps = np.array([20, 28, 40, 56, 80]); convergence = np.array([
        model.moving_average_force(probe_velocity, periods=12, discard=6, steps_per_period=int(n)) for n in steps])
    plt.figure(figsize=(7.3, 4.6)); plt.plot(steps, convergence * 1e22, "o-")
    plt.xlabel("Integration samples per λ period"); plt.ylabel("Mean Fx ($10^{-22}$ N)"); plt.title("Polarization-gradient force-grid convergence")
    save("subdoppler_convergence")

    beta = model.friction_coefficient(probe_velocity, periods=14, discard=7, steps_per_period=48)
    diffusion = model.diffusion_estimate()
    temperature = diffusion / (k_B * beta) if beta > 0 else np.nan
    metadata = {"simulation_version": __version__, "model_fidelity": "Polarization-gradient model adiabatic F=2 to F'=3 population Sisyphus model",
                "isotope": "87Rb", "line": "D2", "detuning_gamma": -3, "saturation_per_beam": .08,
                "phases_rad": PHASES.tolist(), "periods": 12, "discard_periods": 6,
                "temperature_warning": "Einstein estimate only; diffusion is isotropic recoil approximation and not an experimental prediction."}
    np.savez_compressed(OUTPUT / "polarization_gradient_reference.npz", x_m=x, polarization_minus=fractions[-1], polarization_pi=fractions[0],
                        polarization_plus=fractions[1], light_shifts_j=shifts, populations=populations,
                        velocity_m_per_s=velocities, pg_force_n=pg_force, doppler_force_n=d_force,
                        bias_t=biases, bias_friction_kg_per_s=bias_beta, detuning_gamma=detunings,
                        saturation_per_beam=saturations, friction_grid_kg_per_s=beta_grid,
                        convergence_steps=steps, convergence_force_n=convergence, diffusion_kg2_m2_per_s3=diffusion,
                        einstein_temperature_k=temperature, metadata_json=json.dumps(metadata, sort_keys=True))


if __name__ == "__main__": main()
