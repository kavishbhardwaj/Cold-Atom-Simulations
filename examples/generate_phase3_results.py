"""Generate Phase-3 reduced-OBE validation and conditioned trend figures."""

import copy
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.atomic.rb87 import Rb87D2
from cold_atom_mot.atomic.species import get_atomic_line
from cold_atom_mot.io.config import build_effective_model, load_config
from cold_atom_mot.physics.optical_bloch import TwoLevelOBE

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "phase3"
PHASE1 = ROOT / "configs" / "rb87_standard_mot.yaml"
PHASE3 = ROOT / "configs" / "rb87_phase3_obe.yaml"


def save(name: str) -> None:
    plt.tight_layout()
    plt.savefig(OUTPUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close()


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-phase3"
    line = get_atomic_line("87Rb", "D2")
    gamma = line.gamma_rad_s

    detuning = np.linspace(-4, 4, 101)
    saturation = np.geomspace(0.01, 30, 70)
    excited = np.empty((len(saturation), len(detuning)))
    analytic = np.empty_like(excited)
    for i, s in enumerate(saturation):
        for j, delta in enumerate(detuning):
            model = TwoLevelOBE.from_saturation(gamma, delta * gamma, s)
            excited[i, j] = model.steady_state()[1, 1].real
            analytic[i, j] = model.analytic_excited_population()
    plt.figure(figsize=(7.4, 5.2))
    image = plt.pcolormesh(detuning, saturation, excited, shading="auto", cmap="magma")
    plt.yscale("log"); plt.colorbar(image, label="Steady excited population")
    plt.xlabel("Detuning (Γ)"); plt.ylabel("Saturation parameter s")
    plt.title("Reduced OBE steady state and power broadening")
    save("obe_steady_state")

    initial = np.array([[1, 0], [0, 0]], complex)
    plt.figure(figsize=(7.4, 4.8))
    transient_data = {}
    for s in (0.2, 2.0, 10.0):
        model = TwoLevelOBE.from_saturation(gamma, -gamma, s)
        time, density = model.evolve(initial, 12 / gamma, max_step=0.03 / gamma)
        plt.plot(time * gamma, density[:, 1, 1].real, label=f"s={s:g}")
        transient_data[str(s)] = (time, density[:, 1, 1].real)
    plt.xlabel("Time (Γt)"); plt.ylabel("Excited population")
    plt.title("Coherent transients with spontaneous damping (δ=−Γ)")
    plt.legend()
    save("obe_transients")

    phase1 = load_config(PHASE1)
    powers = np.geomspace(1e-4, 0.05, 55)
    detuning_scan = np.linspace(-5.0, -0.25, 60)
    damping = np.empty((len(powers), len(detuning_scan)))
    for i, power in enumerate(powers):
        for j, delta in enumerate(detuning_scan):
            varied = copy.deepcopy(phase1)
            varied["laser"]["power_per_beam_w"] = float(power)
            varied["laser"]["detuning_gamma"] = float(delta)
            damping[i, j] = build_effective_model(varied).linear_coefficients()[0][0]
    plt.figure(figsize=(7.5, 5.2))
    image = plt.pcolormesh(detuning_scan, powers * 1e3, damping * 1e22, shading="auto", cmap="viridis")
    plt.yscale("log"); plt.colorbar(image, label="βx ($10^{-22}$ kg/s)")
    plt.xlabel("Detuning (Γ)"); plt.ylabel("Power per beam (mW)")
    plt.title("Damping is non-monotonic under saturation and power broadening")
    save("damping_power_detuning")

    waists = np.linspace(2e-3, 20e-3, 80)
    x_probe = 2e-3
    fixed_power_force = []
    fixed_intensity_force = []
    reference_waist = phase1["laser"]["waist_m"]
    reference_power = phase1["laser"]["power_per_beam_w"]
    for waist in waists:
        fixed_power = copy.deepcopy(phase1)
        fixed_power["laser"]["waist_m"] = float(waist)
        fixed_power_force.append(build_effective_model(fixed_power).force([x_probe, 0, 0], [0, 0, 0])[0])
        fixed_intensity = copy.deepcopy(fixed_power)
        fixed_intensity["laser"]["power_per_beam_w"] = reference_power * (waist / reference_waist) ** 2
        fixed_intensity_force.append(build_effective_model(fixed_intensity).force([x_probe, 0, 0], [0, 0, 0])[0])
    plt.figure(figsize=(7.4, 4.8))
    plt.plot(waists * 1e3, -np.array(fixed_power_force) * 1e20, label="fixed 10 mW per beam")
    plt.plot(waists * 1e3, -np.array(fixed_intensity_force) * 1e20, label="fixed peak intensity")
    plt.axvline(reference_waist * 1e3, color="black", linestyle="--", label="reference waist")
    plt.xlabel("Gaussian waist (mm)"); plt.ylabel("Restoring magnitude −$F_x$(x=2 mm) ($10^{-20}$ N)")
    plt.title("Beam-waist trends depend on what is held fixed")
    plt.legend()
    save("waist_conditioned_force")

    metadata = {
        "simulation_version": __version__,
        "model_fidelity": "Level C reduced single-transition OBE plus conditioned Level-A scans",
        "atomic_system": "87Rb D2 effective stretched transition",
        "obe_basis": ["ground", "excited"],
        "collapse_operator": "sqrt(Gamma)|g><e|",
        "units": {"power": "W", "waist": "m", "force": "N", "damping": "kg/s"},
        "maximum_obe_analytic_error": float(np.max(abs(excited - analytic))),
    }
    np.savez_compressed(
        OUTPUT / "phase3_reference.npz",
        obe_detuning_gamma=detuning,
        obe_saturation=saturation,
        obe_excited_population=excited,
        obe_analytic_population=analytic,
        damping_power_w=powers,
        damping_detuning_gamma=detuning_scan,
        damping_kg_per_s=damping,
        waist_m=waists,
        fixed_power_force_n=np.array(fixed_power_force),
        fixed_peak_intensity_force_n=np.array(fixed_intensity_force),
        metadata_json=json.dumps(metadata),
    )


if __name__ == "__main__":
    main()
