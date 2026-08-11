"""Generate reproducible plots and a numerical summary for the example trap."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cold_atom import (  # noqa: E402
    BOLTZMANN,
    axial_trap_frequency,
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


def main() -> None:
    """Calculate the example observables and write the result artifacts."""

    output_dir = Path(__file__).resolve().parent
    sigma_v = thermal_velocity_sigma(TEMPERATURE, MASS_RB87)
    omega_r = radial_trap_frequency(TRAP_DEPTH, MASS_RB87, WAIST)
    omega_z = axial_trap_frequency(TRAP_DEPTH, MASS_RB87, WAIST, WAVELENGTH)
    sag = gravitational_sag(GRAVITY, omega_r)

    radius = np.linspace(-100e-6, 100e-6, 500)
    potential = np.array(
        [
            optical_dipole_potential(
                abs(r), 0.0, trap_depth=TRAP_DEPTH, waist=WAIST, wavelength=WAVELENGTH
            )
            for r in radius
        ]
    )

    time = np.linspace(0, 30e-3, 300)
    initial_cloud_sigma = 8e-6
    cloud_sigma = np.sqrt(initial_cloud_sigma**2 + (sigma_v * time) ** 2)

    velocity = np.linspace(-4 * sigma_v, 4 * sigma_v, 500)
    probability_density = np.exp(-0.5 * (velocity / sigma_v) ** 2)
    probability_density /= sigma_v * np.sqrt(2 * np.pi)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    axes[0].plot(
        radius * 1e6,
        potential / BOLTZMANN * 1e3,
        color="#315a8c",
        linewidth=2,
    )
    axes[0].set(
        title="Gaussian dipole trap",
        xlabel="Radial position (µm)",
        ylabel="Potential (mK × $k_B$)",
    )

    axes[1].plot(time * 1e3, cloud_sigma * 1e6, color="#d45b36", linewidth=2)
    axes[1].set(
        title="Ballistic cloud expansion",
        xlabel="Time of flight (ms)",
        ylabel="Cloud width σ (µm)",
    )

    axes[2].plot(
        velocity * 1e3,
        probability_density / 1e3,
        color="#3a8f65",
        linewidth=2,
    )
    axes[2].fill_between(
        velocity * 1e3,
        probability_density / 1e3,
        color="#3a8f65",
        alpha=0.18,
    )
    axes[2].set(
        title="Thermal velocity distribution",
        xlabel="Velocity (mm/s)",
        ylabel="Probability density (s/mm)",
    )

    fig.suptitle("Representative $^{87}$Rb cold-atom simulation (T = 20 µK)", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_dir / "cold_atom_results.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

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
| Cloud σ after 30 ms time of flight | {cloud_sigma[-1] * 1e3:.3f} mm |

The calculations use the idealized models in `cold_atom.py`; they omit atom-atom
interactions, photon scattering, trap anharmonicity, and technical noise.
"""
    (output_dir / "summary.md").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
