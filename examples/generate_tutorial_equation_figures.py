"""Generate the pedagogical equation figures used by docs/tutorial.

These figures are not experimental data. They are direct visualizations of the
analytical equations or the repository's validated low-level solvers, intended
to help a student connect a formula to its physical shape.

Run from the repository root:

    python examples/generate_tutorial_equation_figures.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import k as k_B

from cold_atom_mot.atomic.species import get_atomic_line
from cold_atom_mot.physics.optical_bloch import TwoLevelOBE
from cold_atom_mot.vacuum import loading_curve


OUT = Path("docs/tutorial/figures")
OUT.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.svg")
    plt.close(fig)


def hyperfine_levels():
    line = get_atomic_line("87Rb", "D2")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    x_ground, x_excited = 0.25, 0.72

    g_ref = line.hyperfine_energy_hz("ground", max(line.ground_f))
    e_ref = line.hyperfine_energy_hz("excited", max(line.excited_f))
    ground = {f: (line.hyperfine_energy_hz("ground", f) - g_ref) / 1e9
              for f in line.ground_f}
    excited = {f: (line.hyperfine_energy_hz("excited", f) - e_ref) / 1e6
               for f in line.excited_f}

    # Use separate normalized vertical coordinates so GHz ground and MHz excited
    # splittings are both readable. Labels retain the physical offsets.
    gy = np.linspace(0.3, 0.7, len(ground))
    ey = np.linspace(0.18, 0.82, len(excited))
    for (f, value), y in zip(ground.items(), gy):
        ax.hlines(y, x_ground - 0.12, x_ground + 0.12)
        ax.text(x_ground + 0.15, y, f"F={f}: {value:.6f} GHz", va="center")
    for (f, value), y in zip(excited.items(), ey):
        ax.hlines(y, x_excited - 0.12, x_excited + 0.12)
        ax.text(x_excited + 0.15, y, f"F'={f}: {value:.3f} MHz", va="center")
    ax.text(x_ground, 0.08, "5S1/2", ha="center")
    ax.text(x_excited, 0.08, "5P3/2", ha="center")
    ax.set_title("87Rb D2 hyperfine offsets generated from $A_{hfs}$ and $B_{hfs}$")
    ax.set_xlim(0, 1.15)
    ax.set_ylim(0, 1)
    ax.axis("off")
    save(fig, "hyperfine_energy_levels")


def gaussian_profile():
    u = np.linspace(-2.2, 2.2, 500)
    intensity = np.exp(-2 * u**2)
    fig, ax = plt.subplots()
    ax.plot(u, intensity)
    ax.axvline(1, linestyle="--")
    ax.axhline(np.exp(-2), linestyle="--")
    ax.set_xlabel(r"$r/w$")
    ax.set_ylabel(r"$I/I_0$")
    ax.set_title(r"Gaussian beam: $I/I_0=e^{-2r^2/w^2}$")
    save(fig, "gaussian_beam_profile")


def obe_steady_state():
    detuning = np.linspace(-5, 5, 1001)
    fig, ax = plt.subplots()
    for saturation in (0.1, 1.0, 5.0):
        rho_ee = (saturation / 2) / (
            1 + saturation + (2 * detuning) ** 2
        )
        ax.plot(detuning, rho_ee, label=f"s={saturation:g}")
    ax.set_xlabel(r"detuning $\delta/\Gamma$")
    ax.set_ylabel(r"excited population $\rho_{ee}$")
    ax.set_title("Two-level OBE steady state")
    ax.legend()
    save(fig, "obe_steady_state_lorentzian")


def lindblad_decay():
    gt = np.linspace(0, 6, 500)
    excited = np.exp(-gt)
    fig, ax = plt.subplots()
    ax.plot(gt, excited, label=r"$\rho_{ee}$")
    ax.plot(gt, 1 - excited, label=r"$\rho_{gg}$")
    ax.set_xlabel(r"normalized time $\Gamma t$")
    ax.set_ylabel("population")
    ax.set_title("Lindblad spontaneous decay with the laser off")
    ax.legend()
    save(fig, "lindblad_spontaneous_decay")


def rabi_oscillations():
    # Work in normalized units Gamma=1. The plotted x-axis is therefore Gamma*t.
    model = TwoLevelOBE(gamma=1.0, detuning=0.0, rabi_frequency=3.0)
    initial = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
    time, rho = model.evolve(initial, 12.0, rtol=1e-10, atol=1e-12,
                             max_step=0.01)
    fig, ax = plt.subplots()
    ax.plot(time, rho[:, 1, 1].real)
    ax.set_xlabel(r"normalized time $\Gamma t$")
    ax.set_ylabel(r"excited population $\rho_{ee}$")
    ax.set_title(r"Damped Rabi oscillations: $\Omega=3\Gamma$, $\delta=0$")
    save(fig, "obe_rabi_oscillations")


def normalized_two_beam_force(variable, *, detuning=-2.0, saturation=0.05):
    """Two-beam force/(hbar*k*Gamma), Gamma=k=1, shared denominator."""
    shared = 1 + 2 * saturation
    delta_plus = detuning - variable
    delta_minus = detuning + variable
    r_plus = 0.5 * saturation / (shared + (2 * delta_plus) ** 2)
    r_minus = 0.5 * saturation / (shared + (2 * delta_minus) ** 2)
    return r_plus, r_minus, r_plus - r_minus


def doppler_force():
    velocity = np.linspace(-3, 3, 1001)
    plus, minus, net = normalized_two_beam_force(velocity)
    fig, ax = plt.subplots()
    ax.plot(velocity, plus, label="+k beam")
    ax.plot(velocity, -minus, label="-k beam")
    ax.plot(velocity, net, label="net")
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel(r"velocity $kv/\Gamma$")
    ax.set_ylabel(r"force $/(\hbar k\Gamma)$")
    ax.set_title(r"Doppler damping for two red-detuned beams ($\delta=-2\Gamma$)")
    ax.legend()
    save(fig, "doppler_force_vs_velocity")


def mot_restoring_force():
    # The same Lorentzian imbalance occurs if the dimensionless variable is a
    # position-dependent Zeeman shift rather than kv/Gamma.
    zeeman = np.linspace(-3, 3, 1001)
    _, _, net = normalized_two_beam_force(zeeman)
    fig, ax = plt.subplots()
    ax.plot(zeeman, net)
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel(r"position-induced Zeeman shift $\delta_Z/\Gamma$ (proportional to $x$)")
    ax.set_ylabel(r"force $/(\hbar k\Gamma)$")
    ax.set_title("MOT restoring force from a Zeeman-dependent resonance")
    save(fig, "mot_restoring_force")


def thermal_flux_distribution():
    # u = v/sqrt(2 kBT/m). Bulk Maxwell p(u)=4/sqrt(pi) u^2 exp(-u^2).
    # The normalized surface-flux law is p(u)=2 u^3 exp(-u^2).
    u = np.linspace(0, 4, 1000)
    bulk = 4 / np.sqrt(np.pi) * u**2 * np.exp(-u**2)
    flux = 2 * u**3 * np.exp(-u**2)
    fig, ax = plt.subplots()
    ax.plot(u, bulk, label="bulk Maxwell")
    ax.plot(u, flux, label="surface-crossing flux")
    ax.set_xlabel(r"$u=v/\sqrt{2k_BT/m}$")
    ax.set_ylabel("normalized probability density")
    ax.set_title("Thermal speed distribution seen by a capture surface")
    ax.legend()
    save(fig, "thermal_flux_distribution")


def loading_dynamics():
    time = np.linspace(0, 20, 500)
    one_body = loading_curve(time, loading_rate_s=1.0, one_body_loss_s=0.15)
    two_body = loading_curve(
        time, loading_rate_s=1.0, one_body_loss_s=0.15,
        two_body_coefficient=0.06, effective_volume_m3=1.0,
    )
    fig, ax = plt.subplots()
    ax.plot(time, one_body, label="one-body loss only")
    ax.plot(time, two_body, label="one + two-body loss")
    ax.set_xlabel("time (illustrative units)")
    ax.set_ylabel("atom number (illustrative units)")
    ax.set_title(r"Loading: $\dot N=R-\gamma N-(\beta_2/V)N^2$")
    ax.legend()
    save(fig, "loading_loss_dynamics")


def main():
    hyperfine_levels()
    gaussian_profile()
    obe_steady_state()
    lindblad_decay()
    rabi_oscillations()
    doppler_force()
    mot_restoring_force()
    thermal_flux_distribution()
    loading_dynamics()
    print(f"Wrote tutorial equation figures to {OUT}")


if __name__ == "__main__":
    main()
