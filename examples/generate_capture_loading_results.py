"""Generate trajectory-linked vapour capture and loading diagnostics."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model, build_vapor_state, load_config
from cold_atom_mot.simulation.capture import (
    CaptureCriterion,
    estimate_stratified_vapor_capture_rate,
)
from cold_atom_mot.vacuum import (
    VaporState,
    flux_speed_cdf,
    loading_curve,
    one_sided_thermal_flux_m2_s,
    rubidium_vapor_pressure_pa,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "rb_vapor_loading.yaml"
OUTPUT = ROOT / "results" / "capture_loading"


def save(name):
    plt.tight_layout()
    plt.savefig(OUTPUT / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close()


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-capture-loading"
    config = load_config(CONFIG_PATH)
    force_model = build_effective_model(load_config(ROOT / config["mot_config"]))
    vapor = build_vapor_state(config)
    capture = config["capture"]
    criterion = CaptureCriterion(**capture["criterion"])
    estimate = estimate_stratified_vapor_capture_rate(
        force_model,
        vapor,
        criterion,
        capture_surface_radius_m=capture["surface_radius_m"],
        speed_bin_edges_m_s=capture["speed_bin_edges_m_s"],
        atoms_per_bin=capture["atoms_per_bin"],
        max_step_s=capture["max_step_s"],
        seed=capture["seed"],
    )
    edges = np.asarray(capture["speed_bin_edges_m_s"])
    bin_centres = 0.5 * (edges[:-1] + edges[1:])
    bin_capture = np.array([
        estimate.captured[(estimate.initial_speed_m_s >= low) &
                          (estimate.initial_speed_m_s < high)].mean()
        for low, high in zip(edges[:-1], edges[1:])
    ])

    temperatures = np.linspace(285, 330, 91)
    pressure = np.array([rubidium_vapor_pressure_pa(t) for t in temperatures])
    density = pressure / (1.380649e-23 * temperatures)
    weighted_capture = []
    loading_rate = []
    area = 4 * np.pi * capture["surface_radius_m"]**2
    isotope_fraction = vapor.isotope_fractions[force_model.atom.isotope]
    for temperature, number, fractions in zip(
        temperatures,
        density,
        [np.diff(flux_speed_cdf(edges, t, force_model.atom.mass_kg)) for t in temperatures],
    ):
        probability = float(np.dot(fractions, bin_capture))
        weighted_capture.append(probability)
        flux = area * one_sided_thermal_flux_m2_s(
            isotope_fraction * number, temperature, force_model.atom.mass_kg
        )
        loading_rate.append(flux * probability)
    weighted_capture = np.asarray(weighted_capture)
    loading_rate = np.asarray(loading_rate)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    axes[0].semilogy(temperatures, pressure)
    axes[0].set(xlabel="Cell temperature (K)", ylabel="Rb partial pressure (Pa)",
                title="Sourced natural-Rb vapour pressure")
    axes[1].step(bin_centres, bin_capture, where="mid")
    axes[1].set(xlabel="Incident speed bin centre (m/s)", ylabel="Captured fraction",
                ylim=(-0.03, 1.03), title="Trajectory capture by speed")
    axes[2].semilogy(temperatures, loading_rate)
    axes[2].set(xlabel="Cell temperature (K)", ylabel="Estimated 87Rb loading (atoms/s)",
                title="Flux × simulated capture probability")
    fig.suptitle("Rubidium vapour, rare slow-tail capture, and loading")
    save("vapor_capture_loading")

    loss = config["loading"]
    time = np.linspace(0, loss["curve_duration_s"], loss["curve_points"])
    one_body_rates = np.array([0.05, 0.1, 0.25])
    plt.figure(figsize=(7.4, 4.8))
    curves = []
    for gamma in one_body_rates:
        curve = loading_curve(time, estimate.loading_rate_s, gamma)
        curves.append(curve)
        plt.plot(time, curve / 1e6, label=f"calibrated γ={gamma:g} s⁻¹")
    plt.xlabel("Loading time (s)")
    plt.ylabel("MOT population (millions)")
    plt.title("Loading curves require independently calibrated loss")
    plt.legend()
    save("loading_loss_sensitivity")

    samples_per_bin = np.array([4, 8, 16, capture["atoms_per_bin"]])
    convergence_probability = []
    convergence_error = []
    bin_probabilities = np.diff(
        flux_speed_cdf(edges, vapor.temperature_k, force_model.atom.mass_kg)
    )
    for count in samples_per_bin:
        probabilities = []
        variance = 0.0
        for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
            mask = ((estimate.initial_speed_m_s >= low) &
                    (estimate.initial_speed_m_s < high))
            outcomes = estimate.captured[mask][:count]
            p_bin = outcomes.mean()
            probabilities.append(p_bin)
            variance += bin_probabilities[index]**2 * p_bin * (1 - p_bin) / count
        convergence_probability.append(np.dot(bin_probabilities, probabilities))
        convergence_error.append(np.sqrt(variance))
    convergence_probability = np.asarray(convergence_probability)
    convergence_error = np.asarray(convergence_error)
    plt.figure(figsize=(7.2, 4.6))
    plt.errorbar(
        samples_per_bin,
        convergence_probability,
        yerr=convergence_error,
        marker="o",
        capsize=4,
    )
    plt.xlabel("Trajectories per speed stratum")
    plt.ylabel("Weighted capture probability")
    plt.title("Rare-tail capture convergence and binomial uncertainty")
    save("capture_sampling_convergence")

    metadata = {
        "simulation_version": __version__,
        "model": "effective MOT trajectories plus equilibrium surface flux",
        "configuration": config,
        "criterion": criterion.__dict__,
        "isotope": force_model.atom.isotope,
        "line": force_model.atom.line,
        "limitations": (
            "capture sphere is an acceptance boundary; speeds above the final edge are assigned zero capture; "
            "loss curves use explicitly assumed calibration scenarios"
        ),
    }
    np.savez_compressed(
        OUTPUT / "capture_loading_reference.npz",
        initial_speed_m_s=estimate.initial_speed_m_s,
        captured=estimate.captured,
        sample_weight=estimate.sample_weights,
        capture_time_s=estimate.capture_time_s,
        speed_bin_edges_m_s=edges,
        speed_bin_capture_probability=bin_capture,
        temperatures_k=temperatures,
        rb_partial_pressure_pa=pressure,
        rb_number_density_m3=density,
        weighted_capture_probability=weighted_capture,
        loading_rate_s=loading_rate,
        reference_capture_probability=estimate.capture_probability,
        reference_capture_standard_error=estimate.capture_probability_standard_error,
        reference_incident_flux_s=estimate.incident_flux_s,
        omitted_high_speed_probability=estimate.omitted_high_speed_probability,
        loading_time_s=time,
        assumed_one_body_loss_s=one_body_rates,
        loading_curves=np.asarray(curves),
        convergence_atoms_per_bin=samples_per_bin,
        convergence_capture_probability=convergence_probability,
        convergence_standard_error=convergence_error,
        metadata_json=json.dumps(metadata, sort_keys=True),
    )


if __name__ == "__main__":
    main()
