"""Generate illustrative 87Rb D2 experimental-sequence timing diagnostics."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar, k as k_B

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model, build_multilevel_model, load_config
from cold_atom_mot.simulation.sequence import EddyCurrentResponse, ExperimentalSequence, Ramp, Stage

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT/"results"/"sequence"
BASE = load_config(ROOT/"configs/rb87_d2_mot.yaml")
ATOM = build_effective_model(BASE).atom


def make_sequence(delay=.002, tau=.002, residual=1e-6):
    gamma = ATOM.gamma_rad_s
    def s(name, duration, detuning, power, repump, gradient, kind="step"):
        return Stage(name, duration, Ramp(detuning[0]*gamma, detuning[1]*gamma, kind),
                     Ramp(*power, kind), Ramp(*repump, kind), Ramp(*gradient, kind))
    stages = (s("vapor MOT load", .020, (-2,-2), (.010,.010),(.0005,.0005),(.10,.10)),
              s("compressed MOT", .005, (-2,-3),(.010,.006),(.0005,.0003),(.10,.16),"smooth"),
              s("gradient switch-off", .002,(-3,-4),(.006,.003),(.0003,.0002),(.16,0),"smooth"),
              s("field settling",delay,(-4,-4),(0,0),(0,0),(0,0)),
              s("PGC ramp",.002,(-4,-6),(.003,.0015),(.0002,.0001),(0,0),"smooth"),
              s("molasses hold",.006,(-6,-6),(.0015,.0015),(.0001,.0001),(0,0)),
              s("release / TOF",.010,(-6,-6),(0,0),(0,0),(0,0)))
    response = EddyCurrentResponse(.025,.16,.0005,[0,0,1e-5],tau,[0,0,residual])
    return ExperimentalSequence(stages,response)


def cooling_proxy(delay, tau, residual, hold=.006):
    """Effective Doppler/recoil proxy, explicitly not coherent PGC temperature."""
    sequence=make_sequence(delay,tau,residual); model=build_effective_model(BASE)
    probe_time=.027+delay+.002+min(hold,.006)
    snapshot=sequence.snapshot(min(probe_time,sequence.duration))
    model.beams=sequence.apply_beams(model.beams,min(probe_time,sequence.duration))
    model.magnetic_field=sequence.magnetic_response
    dv=1e-3; beta=-(model.force(np.zeros(3),[dv,0,0],probe_time)[0]-model.force(np.zeros(3),[-dv,0,0],probe_time)[0])/(2*dv)
    rates=model.scattering_rates(np.zeros(3),np.zeros(3),probe_time).sum()
    diffusion=(hbar*ATOM.wave_number_rad_m)**2*rates*(2/3)/2
    temperature=diffusion/(k_B*beta) if beta>0 else np.nan
    return beta,temperature,np.linalg.norm(sequence.magnetic_response.field(np.zeros(3),probe_time)),snapshot


def main():
    OUT.mkdir(parents=True,exist_ok=True); sequence=make_sequence()
    time=np.linspace(0,sequence.duration,800); controls=[sequence.snapshot(t) for t in time]
    power=np.array([x.cooling_power for x in controls]); detuning=np.array([x.cooling_detuning for x in controls])/ATOM.gamma_rad_s
    gradient=np.array([x.quadrupole_gradient for x in controls]); field=np.array([np.linalg.norm(sequence.magnetic_response.field(np.zeros(3),t)) for t in time])
    model=build_effective_model(BASE); sample_index=np.arange(0,len(time),8)
    beta_sample=[]; scatter_sample=[]
    for index in sample_index:
        t=time[index]; model.beams=sequence.apply_beams(build_effective_model(BASE).beams,t); model.magnetic_field=sequence.magnetic_response
        dv=1e-3; beta_sample.append(-(model.force(np.zeros(3),[dv,0,0],t)[0]-model.force(np.zeros(3),[-dv,0,0],t)[0])/(2*dv))
        scatter_sample.append(model.scattering_rates(np.zeros(3),np.zeros(3),t).sum())
    beta=np.interp(time,time[sample_index],beta_sample); scattering=np.interp(time,time[sample_index],scatter_sample)
    variance=np.empty(len(time)); variance[0]=k_B*300e-6/ATOM.mass_kg
    for i,dt in enumerate(np.diff(time),1):
        diffusion=(hbar*ATOM.wave_number_rad_m/ATOM.mass_kg)**2*scattering[i-1]*(2/3)
        variance[i]=max(0,variance[i-1]+(-2*beta[i-1]/ATOM.mass_kg*variance[i-1]+diffusion)*dt)
    temperature=ATOM.mass_kg*variance/k_B; rms_velocity=np.sqrt(variance)
    rms_position=np.zeros(len(time))
    rms_position[1:]=np.cumsum((rms_velocity[:-1]+rms_velocity[1:])*np.diff(time)/2)
    rate=build_multilevel_model(load_config(ROOT/"configs/rb87_d2_multilevel.yaml")); population_sample=[]
    population_index=np.arange(0,len(time),16)
    for index in population_index:
        t=time[index]; rate.beam_families=sequence.apply_beam_families(rate.beam_families,t); rate.magnetic_field=sequence.magnetic_response
        if sum(f.beam.power for f in rate.beam_families) > 0:
            population_sample.append(rate.steady_state(np.zeros(3),np.zeros(3),t))
        else:
            population_sample.append(population_sample[-1].copy())
    population_sample=np.asarray(population_sample)
    fig,axes=plt.subplots(6,1,figsize=(10,12),sharex=True)
    axes[0].plot(time*1e3,power*1e3);axes[0].set_ylabel("Cooling power/beam (mW)")
    axes[1].plot(time*1e3,detuning);axes[1].set_ylabel("Detuning (Gamma)")
    axes[2].plot(time*1e3,gradient);axes[2].set_ylabel("Gradient (T/m)")
    axes[3].semilogy(time*1e3,np.maximum(field*1e6,1e-4));axes[3].set(ylabel="|B(0)| (uT)",xlabel="Sequence time (ms)")
    axes[4].plot(time*1e3,rms_velocity);axes[4].set_ylabel("RMS vx proxy (m/s)")
    axes[5].semilogy(time*1e3,np.maximum(temperature*1e6,1e-3));axes[5].set(ylabel="Kinetic proxy (uK)",xlabel="Sequence time (ms)")
    fig.suptitle("Illustrative 87Rb D2 sequence controls (not an experimental optimum)");fig.tight_layout();fig.savefig(OUT/"sequence_timeline.png",dpi=220);plt.close(fig)
    delays=np.linspace(0,.012,25); taus=np.array([.0005,.002,.005]); residuals=np.array([0,1e-6,5e-6])
    proxy=np.empty((len(taus),len(residuals),len(delays),3))
    for i,tau in enumerate(taus):
      for j,residual in enumerate(residuals):
       for k,delay in enumerate(delays): proxy[i,j,k,:3]=cooling_proxy(delay,tau,residual)[:3]
    durations=np.linspace(.001,.012,20); duration_proxy=np.array([cooling_proxy(.002,.002,1e-6,d)[:3] for d in durations])
    fig,axes=plt.subplots(1,2,figsize=(12,4.6))
    for i,tau in enumerate(taus): axes[0].plot(delays*1e3,proxy[i,1,:,1]*1e6,label=f"tau_eddy={tau*1e3:g} ms")
    axes[0].set(xlabel="Wait after switch-off (ms)",ylabel="Effective recoil/friction proxy (uK)",title="Timing sensitivity; residual DC=10 mG");axes[0].legend()
    axes[1].plot(durations*1e3,duration_proxy[:,1]*1e6);axes[1].set(xlabel="Molasses duration (ms)",ylabel="Effective proxy (uK)",title="Duration diagnostic")
    fig.tight_layout();fig.savefig(OUT/"sequence_timing_sensitivity.png",dpi=220);plt.close(fig)
    metadata={"simulation_version":__version__,"fidelity":"effective two-level local damping plus recoil-event diffusion proxy",
              "warning":"Not a coherent PGC temperature or experimental optimum; no internal-state/dipole-force diffusion.","mains_amplitude_t":0.0,"random_seed":None}
    np.savez_compressed(OUT/"sequence_reference.npz",time_s=time,cooling_power_w=power,detuning_gamma=detuning,gradient_t_m=gradient,field_t=field,
      rms_velocity_m_s=rms_velocity,rms_position_upper_proxy_m=rms_position,
      kinetic_temperature_proxy_k=temperature,scattering_rate_s=scattering,
      population_time_s=time[population_index],rate_equation_population=population_sample,
      delay_s=delays,eddy_tau_s=taus,residual_t=residuals,delay_proxy=proxy,molasses_duration_s=durations,duration_proxy=duration_proxy,metadata_json=json.dumps(metadata,sort_keys=True))


if __name__=="__main__": main()
