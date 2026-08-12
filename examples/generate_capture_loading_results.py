"""Generate trajectory-linked vapour capture and loading diagnostics."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model, build_multilevel_model, build_vapor_state, load_config
from cold_atom_mot.solvers.deterministic import integrate_trajectory
from cold_atom_mot.simulation.capture import (
    CaptureCriterion,
    capture_response_map,
    estimate_adaptive_vapor_capture_rate,
)
from cold_atom_mot.vacuum import (
    flux_speed_pdf,
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
    estimate = estimate_adaptive_vapor_capture_rate(
        force_model,
        vapor,
        criterion,
        capture_surface_radius_m=capture["surface_radius_m"],
        initial_speed_edges_m_s=capture["speed_bin_edges_m_s"],
        maximum_speed_m_s=capture["maximum_speed_m_s"],
        tail_relative_loading_tolerance=capture["tail_relative_loading_tolerance"],
        atoms_per_bin=capture["atoms_per_bin"],
        max_step_s=capture["max_step_s"],
        seed=capture["seed"], confidence_level=capture["confidence_level"],
        rtol=capture["rtol"], atol=capture["atol"],
    )
    edges = np.asarray(capture["speed_bin_edges_m_s"])
    bin_centres = 0.5 * (edges[:-1] + edges[1:])
    bin_capture = np.array([
        estimate.captured[(estimate.initial_speed_m_s >= low) &
                          (estimate.initial_speed_m_s < high)].mean()
        for low, high in zip(edges[:-1], edges[1:])
    ])

    reservoir_temperatures = np.linspace(298.15, 330, 65)
    vapor_temperatures = np.array([285.0, 300.0, 330.0])
    pressure = np.array([rubidium_vapor_pressure_pa(t) for t in reservoir_temperatures])
    density = pressure / (1.380649e-23 * vapor.vapor_temperature_k)
    # Temperature-independent response map; thermal integration changes with T.
    response_speeds = np.linspace(0.5, 50.0, 12)
    impact_edges = capture["surface_radius_m"] * np.sqrt(np.linspace(0, 1, 4))
    response, response_low, response_high, impact_weights = capture_response_map(
        force_model, criterion, response_speeds, impact_edges,
        capture_surface_radius_m=capture["surface_radius_m"], samples_per_cell=4,
        max_step_s=capture["max_step_s"], seed=capture["seed"]+101,
        rtol=capture["rtol"], atol=capture["atol"],
    )
    speed_response = response @ impact_weights
    speed_grid = np.linspace(0, response_speeds[-1], 800)
    response_grid = np.interp(speed_grid, response_speeds, speed_response, left=speed_response[0], right=0)
    thermal_capture = np.array([np.trapezoid(response_grid*flux_speed_pdf(speed_grid,t,force_model.atom.mass_kg),speed_grid) for t in vapor_temperatures])
    coarse_grid = speed_grid[::4]
    # Quadrature-grid refinement holds the independently simulated response map
    # fixed; response-node convergence is separately visible in the Wilson bands.
    coarse_response = np.interp(coarse_grid, response_speeds, speed_response, left=speed_response[0], right=0)
    coarse_capture = np.trapezoid(coarse_response * flux_speed_pdf(coarse_grid, vapor.vapor_temperature_k, force_model.atom.mass_kg), coarse_grid)
    velocity_grid_relative_change = abs(coarse_capture - thermal_capture[1]) / max(thermal_capture[1], 1e-30)
    reference_response_probability = float(np.interp(vapor.vapor_temperature_k,vapor_temperatures,thermal_capture))
    area = 4*np.pi*capture["surface_radius_m"]**2
    isotope_fraction=vapor.isotope_fractions[force_model.atom.isotope]
    reference_flux_per_density = area*one_sided_thermal_flux_m2_s(1.0,vapor.vapor_temperature_k,force_model.atom.mass_kg)
    loading_rate = isotope_fraction*density*reference_flux_per_density*reference_response_probability

    fig,axes=plt.subplots(2,2,figsize=(12,9))
    pressure_line = axes[0,0].semilogy(reservoir_temperatures,pressure,label="Rb partial pressure",color="tab:blue")
    density_axis = axes[0,0].twinx()
    density_line = density_axis.semilogy(reservoir_temperatures,density,label="Rb number density",color="tab:orange")
    axes[0,0].set(xlabel="Rb reservoir/cold-spot temperature (K)",ylabel="Pressure (Pa)",title="Validated Alcock pressure and ideal-gas density")
    density_axis.set_ylabel("Number density (m⁻³)")
    axes[0,0].legend(pressure_line+density_line,[line.get_label() for line in pressure_line+density_line],loc="upper left")
    bulk_speed=np.linspace(0,800,500); a=force_model.atom.mass_kg/(2*1.380649e-23*vapor.vapor_temperature_k)
    bulk=4/np.sqrt(np.pi)*a**1.5*bulk_speed**2*np.exp(-a*bulk_speed**2)
    incident=flux_speed_pdf(bulk_speed,vapor.vapor_temperature_k,force_model.atom.mass_kg)
    axes[0,1].plot(bulk_speed,bulk,label="bulk Maxwell"); axes[0,1].plot(bulk_speed,incident,label="surface flux")
    axes[0,1].set(xlabel="Speed (m/s)",ylabel="Probability density (s/m)",title="Bulk versus incident atoms"); axes[0,1].legend()
    axes[1,0].plot(response_speeds,speed_response,"o-")
    axes[1,0].set(xlabel="Incident speed (m/s)",ylabel="Impact-averaged capture",title="Trajectory response Pcapture(v)")
    axes[1,1].semilogy(reservoir_temperatures,loading_rate)
    axes[1,1].set(xlabel="Reservoir temperature (K)",ylabel="Estimated loading (atoms/s)",title="Density × flux × thermal response")
    fig.suptitle("Reservoir pressure, vapour kinetics, capture response, and loading")
    save("vapor_capture_loading")

    plt.figure(figsize=(7.4, 4.8))
    impact_centres = 0.5 * (impact_edges[:-1] + impact_edges[1:])
    mesh = plt.pcolormesh(impact_centres * 1e3, response_speeds, response, shading="nearest", vmin=0, vmax=1, cmap="viridis")
    plt.colorbar(mesh, label="Capture probability")
    plt.xlabel("Impact-parameter bin edge (mm)")
    plt.ylabel("Incident speed (m/s)")
    plt.title("Trajectory response versus speed and impact parameter")
    save("capture_response_map")

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

    # Independent seeds: no prefixes from a single realization.
    samples_per_bin=np.array([4,8,16,32,64]); repetitions=2
    convergence_probability=np.empty((len(samples_per_bin),repetitions))
    convergence_loading=np.empty_like(convergence_probability)
    convergence_interval=np.empty((len(samples_per_bin),repetitions,2))
    for i,count in enumerate(samples_per_bin):
        for j in range(repetitions):
            trial=estimate_adaptive_vapor_capture_rate(
                force_model,vapor,criterion,capture_surface_radius_m=capture["surface_radius_m"],
                initial_speed_edges_m_s=capture["speed_bin_edges_m_s"],atoms_per_bin=int(count),
                maximum_speed_m_s=capture["maximum_speed_m_s"],
                tail_relative_loading_tolerance=capture["tail_relative_loading_tolerance"],
                max_step_s=capture["max_step_s"],seed=capture["seed"]+1000*i+j,
                confidence_level=capture["confidence_level"],rtol=capture["rtol"],atol=capture["atol"])
            convergence_probability[i,j]=trial.capture_probability
            convergence_loading[i,j]=trial.loading_rate_s
            convergence_interval[i,j]=trial.confidence_interval
    mean_probability=convergence_probability.mean(axis=1); spread_probability=convergence_probability.std(axis=1,ddof=1)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    axes[0].errorbar(samples_per_bin,mean_probability,yerr=spread_probability,marker="o",capsize=4,label="between-seed spread")
    axes[0].fill_between(samples_per_bin,convergence_interval[:,:,0].mean(axis=1),
                         convergence_interval[:,:,1].mean(axis=1),alpha=.18,label="mean 95% Wilson envelope")
    axes[0].set(xlabel="Trajectories per speed stratum",ylabel="Weighted capture probability",title="Independent-seed convergence")
    axes[0].legend()
    axes[1].errorbar(samples_per_bin,convergence_loading.mean(axis=1),yerr=convergence_loading.std(axis=1,ddof=1),marker="o",capsize=4)
    axes[1].set(xlabel="Trajectories per speed stratum",ylabel="Loading rate (atoms/s)",title="Loading-rate convergence")
    save("capture_sampling_convergence")

    # Matched selected trajectories: effective versus multilevel force.
    multilevel=build_multilevel_model(load_config(ROOT/"configs"/"rb87_d2_multilevel.yaml"))
    comparison_speeds=np.array([5.0,10.0,20.0])
    fig,axes=plt.subplots(1,2,figsize=(11,4.5))
    effective_final=[]; multilevel_final=[]
    for speed in comparison_speeds:
        effective_tr=integrate_trajectory(force_model,[-capture["surface_radius_m"],0,0],[speed,0,0],.002,max_step=1e-4,rtol=1e-6,atol=1e-9,sample_step=2e-5)
        multilevel_tr=integrate_trajectory(multilevel,[-capture["surface_radius_m"],0,0],[speed,0,0],.002,max_step=1e-4,rtol=1e-5,atol=1e-8,sample_step=2e-5)
        axes[0].plot(effective_tr.time*1e3,effective_tr.position[:,0]*1e3,label=f"effective {speed:g} m/s")
        axes[0].plot(multilevel_tr.time*1e3,multilevel_tr.position[:,0]*1e3,"--",label=f"multilevel {speed:g} m/s")
        effective_final.append(np.linalg.norm(effective_tr.velocity[-1])); multilevel_final.append(np.linalg.norm(multilevel_tr.velocity[-1]))
    axes[0].set(xlabel="Time (ms)",ylabel="x (mm)",title="Representative matched trajectories"); axes[0].legend(ncol=2,fontsize=8)
    axes[1].plot(comparison_speeds,effective_final,"o-",label="effective")
    axes[1].plot(comparison_speeds,multilevel_final,"s--",label="multilevel rate equation")
    axes[1].set(xlabel="Initial speed (m/s)",ylabel="Final speed after 2 ms (m/s)",title="Model-fidelity sensitivity"); axes[1].legend()
    save("effective_multilevel_capture_comparison")

    # Beam-waist dependence with held-fixed 10 mW/beam and all other MOT inputs.
    waist_values=np.array([4,6,8,12,16])*1e-3; waist_loading=[]
    mot_base=load_config(ROOT/config["mot_config"])
    for waist in waist_values:
        varied=dict(mot_base); varied["laser"]=dict(mot_base["laser"]); varied["laser"]["waist_m"]=float(waist)
        varied_model=build_effective_model(varied)
        local_response,_,_,local_weights=capture_response_map(varied_model,criterion,response_speeds,impact_edges,capture_surface_radius_m=capture["surface_radius_m"],samples_per_cell=4,max_step_s=capture["max_step_s"],seed=capture["seed"]+int(waist*1e6),rtol=capture["rtol"],atol=capture["atol"])
        local_speed=local_response@local_weights
        local_grid=np.interp(speed_grid,response_speeds,local_speed,left=local_speed[0],right=0)
        local_probability=np.trapezoid(local_grid*flux_speed_pdf(speed_grid,vapor.vapor_temperature_k,force_model.atom.mass_kg),speed_grid)
        waist_loading.append(estimate.incident_flux_s*local_probability)
    plt.figure(figsize=(7.2,4.5)); plt.plot(waist_values*1e3,waist_loading,"o-")
    plt.xlabel("Beam waist (mm)"); plt.ylabel("Estimated loading (atoms/s)")
    plt.title("Fixed 10 mW/beam: intensity–capture-volume tradeoff")
    save("loading_vs_beam_waist")

    metadata = {
        "simulation_version": __version__,
        "model": "effective MOT trajectories plus equilibrium surface flux",
        "configuration": config,
        "criterion": criterion.__dict__,
        "isotope": force_model.atom.isotope,
        "line": force_model.atom.line,
        "limitations": (
            "the capture sphere is an acceptance boundary; the adaptive tail is bounded statistically, not "
            "declared identically zero; loss curves use explicitly assumed calibration scenarios"
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
        reservoir_temperatures_k=reservoir_temperatures,
        vapor_temperatures_k=vapor_temperatures,
        rb_partial_pressure_pa=pressure,
        rb_number_density_m3=density,
        thermal_capture_probability=thermal_capture,
        velocity_grid_relative_change=velocity_grid_relative_change,
        capture_response_speeds_m_s=response_speeds,
        capture_impact_edges_m=impact_edges,
        capture_response_map=response,
        capture_response_low=response_low,
        capture_response_high=response_high,
        loading_rate_vs_reservoir_temperature_s=loading_rate,
        reference_capture_probability=estimate.capture_probability,
        reference_capture_standard_error=estimate.capture_probability_standard_error,
        reference_incident_flux_s=estimate.incident_flux_s,
        omitted_high_speed_probability=estimate.omitted_high_speed_probability,
        last_simulated_speed_m_s=estimate.last_simulated_speed_m_s,
        omitted_capture_probability_upper=estimate.omitted_capture_probability_upper,
        omitted_loading_rate_upper_s=estimate.omitted_loading_rate_upper_s,
        tail_converged=estimate.tail_converged,
        reference_confidence_interval=estimate.confidence_interval,
        reference_loading_rate_confidence_interval_s=estimate.loading_rate_confidence_interval_s,
        loading_time_s=time,
        assumed_one_body_loss_s=one_body_rates,
        loading_curves=np.asarray(curves),
        convergence_atoms_per_bin=samples_per_bin,
        convergence_capture_probability=convergence_probability,
        convergence_loading_rate_s=convergence_loading,
        convergence_confidence_interval=convergence_interval,
        comparison_initial_speed_m_s=comparison_speeds,
        effective_final_speed_m_s=effective_final,
        multilevel_final_speed_m_s=multilevel_final,
        waist_m=waist_values,
        loading_vs_waist_s=waist_loading,
        metadata_json=json.dumps(metadata, sort_keys=True),
    )


if __name__ == "__main__":
    main()
