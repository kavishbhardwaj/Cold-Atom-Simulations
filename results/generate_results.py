"""Generate reproducible figures and a numerical summary for the example trap."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cold_atom import (  # noqa: E402
    BOLTZMANN,
    axial_trap_frequency,
    gaussian_beam_waist,
    gravitational_sag,
    optical_dipole_potential,
    radial_trap_frequency,
    thermal_velocity_sigma,
)

MASS_RB87 = 1.44316060e-25
TEMPERATURE = 20e-6
TRAP_DEPTH = BOLTZMANN * 1e-3
WAIST = 50e-6
WAVELENGTH = 1064e-9
GRAVITY = 9.80665
COLORS = ("#315a8c", "#d45b36", "#3a8f65")


def save_figure(output_dir: Path, name: str) -> None:
    """Save the current figure as display-ready PNG and vector SVG files."""

    plt.tight_layout()
    plt.savefig(
        output_dir / f"{name}.png",
        format="png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.savefig(
        output_dir / f"{name}.svg",
        format="svg",
        bbox_inches="tight",
        metadata={"Date": None},
    )
    plt.close()


def generate_figures(output_dir: Path, sigma_v: float) -> None:
    """Generate six complementary views of the representative simulation."""

    radius = np.linspace(-100e-6, 100e-6, 500)
    potential = np.array(
        [
            optical_dipole_potential(
                abs(r), 0.0, trap_depth=TRAP_DEPTH, waist=WAIST, wavelength=WAVELENGTH
            )
            for r in radius
        ]
    )
    plt.figure(figsize=(7, 4.5))
    plt.plot(radius * 1e6, potential / BOLTZMANN * 1e3, color=COLORS[0], linewidth=2)
    plt.title("Gaussian optical-dipole potential")
    plt.xlabel("Radial position (µm)")
    plt.ylabel("Potential (mK × $k_B$)")
    save_figure(output_dir, "trap_potential")

    z = np.linspace(-20e-3, 20e-3, 500)
    beam_radius = np.array(
        [gaussian_beam_waist(value, WAIST, WAVELENGTH) for value in z]
    )
    plt.figure(figsize=(7, 4.5))
    plt.plot(z * 1e3, beam_radius * 1e6, color=COLORS[1], linewidth=2)
    plt.title("Gaussian beam propagation")
    plt.xlabel("Axial position (mm)")
    plt.ylabel("Beam radius (µm)")
    save_figure(output_dir, "beam_waist")

    time = np.linspace(0, 30e-3, 300)
    cloud_sigma = np.sqrt((8e-6) ** 2 + (sigma_v * time) ** 2)
    plt.figure(figsize=(7, 4.5))
    plt.plot(time * 1e3, cloud_sigma * 1e6, color=COLORS[1], linewidth=2)
    plt.title("Ballistic cloud expansion")
    plt.xlabel("Time of flight (ms)")
    plt.ylabel("Cloud width σ (µm)")
    save_figure(output_dir, "time_of_flight")

    velocity = np.linspace(-4 * sigma_v, 4 * sigma_v, 500)
    density = np.exp(-0.5 * (velocity / sigma_v) ** 2)
    density /= sigma_v * np.sqrt(2 * np.pi)
    plt.figure(figsize=(7, 4.5))
    plt.plot(velocity * 1e3, density / 1e3, color=COLORS[2], linewidth=2)
    plt.fill_between(velocity * 1e3, density / 1e3, color=COLORS[2], alpha=0.18)
    plt.title("Thermal velocity distribution")
    plt.xlabel("Velocity (mm/s)")
    plt.ylabel("Probability density (s/mm)")
    save_figure(output_dir, "thermal_velocity")

    depths = BOLTZMANN * np.linspace(0.05e-3, 2e-3, 300)
    radial = np.array([radial_trap_frequency(d, MASS_RB87, WAIST) for d in depths])
    axial = np.array([axial_trap_frequency(d, MASS_RB87, WAIST, WAVELENGTH) for d in depths])
    plt.figure(figsize=(7, 4.5))
    plt.plot(
        depths / BOLTZMANN * 1e3,
        radial / (2 * np.pi),
        label="Radial",
        color=COLORS[0],
    )
    plt.plot(
        depths / BOLTZMANN * 1e3,
        axial / (2 * np.pi),
        label="Axial",
        color=COLORS[1],
    )
    plt.yscale("log")
    plt.title("Trap frequencies versus depth")
    plt.xlabel("Trap depth (mK × $k_B$)")
    plt.ylabel("Frequency (Hz)")
    plt.legend()
    save_figure(output_dir, "trap_frequencies")

    sag = np.array([gravitational_sag(GRAVITY, omega) for omega in radial])
    plt.figure(figsize=(7, 4.5))
    plt.plot(depths / BOLTZMANN * 1e3, sag * 1e6, color=COLORS[2], linewidth=2)
    plt.title("Gravitational sag versus trap depth")
    plt.xlabel("Trap depth (mK × $k_B$)")
    plt.ylabel("Radial sag (µm)")
    save_figure(output_dir, "gravitational_sag")


def main() -> None:
    """Calculate observables and write all result artifacts."""

    output_dir = Path(__file__).resolve().parent
    sigma_v = thermal_velocity_sigma(TEMPERATURE, MASS_RB87)
    omega_r = radial_trap_frequency(TRAP_DEPTH, MASS_RB87, WAIST)
    omega_z = axial_trap_frequency(TRAP_DEPTH, MASS_RB87, WAIST, WAVELENGTH)
    sag = gravitational_sag(GRAVITY, omega_r)
    cloud_sigma_30_ms = np.sqrt((8e-6) ** 2 + (sigma_v * 30e-3) ** 2)

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-simulations"
    generate_figures(output_dir, sigma_v)

    summary = f"""# Representative simulation results

These deterministic results use rubidium-87 atoms in an ideal 1064 nm Gaussian
optical dipole trap. They can be regenerated with `python results/generate_results.py`.

| Quantity | Value |
| --- | ---: |
| Atom temperature | {TEMPERATURE * 1e6:.1f} µK |
| Trap depth | {TRAP_DEPTH / BOLTZMANN * 1e3:.1f} mK × k_B |
| Beam waist | {WAIST * 1e6:.1f} µm |
| 1D thermal velocity σ | {sigma_v * 1e3:.2f} mm/s |
| Radial trap frequency ω_r / 2π | {omega_r / (2 * np.pi):.1f} Hz |
| Axial trap frequency ω_z / 2π | {omega_z / (2 * np.pi):.1f} Hz |
| Radial gravitational sag | {sag * 1e6:.3f} µm |
| Cloud σ after 30 ms time of flight | {cloud_sigma_30_ms * 1e3:.3f} mm |

The calculations use the idealized models in `cold_atom.py`; they omit atom-atom
interactions, photon scattering, trap anharmonicity, and technical noise.
"""
    (output_dir / "summary.md").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
