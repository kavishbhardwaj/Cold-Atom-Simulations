"""Generate constrained damping/power and beam-waist/capture studies."""
import copy,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar
from cold_atom_mot import __version__
from cold_atom_mot.io.config import load_config,build_effective_model
from cold_atom_mot.simulation.capture import CaptureCriterion,evaluate_capture

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results"/"parameter_studies"; CONFIG=ROOT/"configs"/"rb87_d2_mot.yaml"
def save(name):
 plt.tight_layout(); plt.savefig(OUT/f"{name}.png",dpi=220,bbox_inches="tight"); plt.savefig(OUT/f"{name}.svg",bbox_inches="tight",metadata={"Date":None}); plt.close()

def main():
 OUT.mkdir(parents=True,exist_ok=True); plt.style.use("seaborn-v0_8-whitegrid"); plt.rcParams["svg.hashsalt"]="mot-parameter-studies"
 base=load_config(CONFIG); model=build_effective_model(base); atom=model.atom
 s=np.geomspace(.002,100,90); detunings=np.array([-1.,-2.,-4.]); v=np.linspace(-2,2,201)
 beta=np.empty((len(detunings),len(s))); beta_analytic=np.empty_like(beta); scatter=np.empty_like(beta); excited=np.empty_like(beta)
 for j,delta in enumerate(detunings):
  for i,si in enumerate(s):
   c=copy.deepcopy(base); w=c["laser"]["waist_m"]; c["laser"]["power_per_beam_w"]=si*atom.saturation_intensity_w_m2*np.pi*w*w/2; c["laser"]["detuning_gamma"]=float(delta)
   m=build_effective_model(c); beta[j,i]=m.linear_coefficients(velocity_step=2e-4)[0][0]
   denominator=1+6*si+4*delta*delta; beta_analytic[j,i]=-8*hbar*atom.wave_number_rad_m**2*si*delta/denominator**2
   rates=m.scattering_rates(np.zeros(3),np.zeros(3)); scatter[j,i]=rates.sum(); excited[j,i]=scatter[j,i]/atom.gamma_rad_s
 fig,axes=plt.subplots(1,3,figsize=(14,4.5))
 for j,d in enumerate(detunings): axes[0].loglog(s,beta[j]*1e22,label=f"δ={d:g}Γ"); axes[1].semilogx(s,scatter[j]/1e6,label=f"δ={d:g}Γ")
 axes[0].loglog(s,beta_analytic[1]*1e22,"k--",label="analytic δ=−2Γ"); axes[0].set(xlabel="On-axis saturation per beam",ylabel="β ($10^{-22}$ kg/s)",title="Velocity-slope damping"); axes[0].legend()
 axes[1].set(xlabel="Saturation per beam",ylabel="Total scattering (MHz)",title="Scattering keeps saturating"); axes[1].legend()
 powers=(.02,.2,2,20)
 for si in powers:
  c=copy.deepcopy(base); w=c["laser"]["waist_m"]; c["laser"]["power_per_beam_w"]=si*atom.saturation_intensity_w_m2*np.pi*w*w/2
  m=build_effective_model(c); axes[2].plot(v,[m.force([0,0,0],[vv,0,0])[0]*1e21 for vv in v],label=f"s={si:g}")
 axes[2].set(xlabel="vx (m/s)",ylabel="Fx ($10^{-21}$ N)",title="Power broadens F(v)"); axes[2].legend(); fig.suptitle("More scattering is not necessarily more damping"); save("damping_power_physics")

 waists=np.array([3,5,8,12,16])*1e-3; speeds=np.linspace(5,100,12); fractions=[]; capture_velocity=[]; peak=[]
 criterion=CaptureCriterion(.002,.5,.006,.0005)
 for waist in waists:
  c=copy.deepcopy(base); c["laser"]["waist_m"]=float(waist); m=build_effective_model(c); peak.append(m.beams[0].peak_intensity)
  positions=np.tile([-min(waist,.01)*.7,0,0],(len(speeds),1)); velocities=np.column_stack([speeds,np.zeros((len(speeds),2))])
  caught,_=evaluate_capture(m,positions,velocities,criterion,max_step=2e-5); fractions.append(caught.mean()); capture_velocity.append(speeds[np.flatnonzero(caught)[-1]] if caught.any() else 0)
 fig,axes=plt.subplots(1,3,figsize=(13,4.3)); axes[0].plot(waists*1e3,np.array(peak),"o-"); axes[0].set(xlabel="Waist (mm)",ylabel="Peak intensity (W/m²)",title="Fixed 10 mW per beam")
 axes[1].plot(waists*1e3,capture_velocity,"o-"); axes[1].set(xlabel="Waist (mm)",ylabel="Captured-speed threshold (m/s)",title="Trajectory criterion")
 axes[2].plot(waists*1e3,fractions,"o-"); axes[2].set(xlabel="Waist (mm)",ylabel="Captured fraction of speed grid",title="Local acceptance proxy")
 fig.suptitle("Beam-waist tradeoff: local intensity versus interaction distance"); save("beam_waist_capture")
 metadata={"version":__version__,"model":"effective semiclassical MOT","held_fixed":"10 mW per beam; detuning and gradient from rb87_d2_mot.yaml","capture_criterion":criterion.__dict__,"limitation":"speed-grid acceptance proxy, not vapour loading or atom number"}
 np.savez_compressed(OUT/"parameter_studies.npz",saturation=s,detuning_gamma=detunings,beta_kg_s=beta,beta_analytic_kg_s=beta_analytic,scattering_s=scatter,excited_proxy=excited,waist_m=waists,peak_intensity_w_m2=peak,capture_velocity_m_s=capture_velocity,capture_fraction=fractions,metadata_json=json.dumps(metadata,sort_keys=True))
if __name__=="__main__": main()
