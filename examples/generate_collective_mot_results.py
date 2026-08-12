"""Generate optional Gaussian mean-field collective-MOT diagnostics."""
import csv,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar,pi

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_effective_model,load_config
from cold_atom_mot.physics.collective import GaussianCloud,MultipleScatteringModel,CollectiveLoading

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results"/"collective_mot"


def main():
    OUT.mkdir(parents=True,exist_ok=True); force=build_effective_model(load_config(ROOT/"configs/rb87_d2_mot.yaml")); atom=force.atom
    saturation=sum(b.peak_intensity for b in force.beams)/atom.saturation_intensity_w_m2; detuning=-2
    sigma0=3*atom.wavelength_m**2/(2*pi); sigma_laser=sigma0/(1+saturation+4*detuning**2)
    # Effective reabsorption cross section is explicitly an input scenario. The
    # resonant sigma0 choice is an upper, not fitted, reference case.
    scattering=MultipleScatteringModel(sigma_laser,sigma0,sum(b.peak_intensity for b in force.beams),atom.wave_number_rad_m,1e6)
    populations=np.geomspace(1e3,1e9,100); kappa=2e-19; temperature=150e-6
    sigma=np.array([scattering.equilibrium_sigma(n,kappa,temperature,atom.mass_kg) for n in populations])
    clouds=[GaussianCloud(n,[s]*3,[temperature]*3,atom.mass_kg) for n,s in zip(populations,sigma)]
    density=np.array([c.peak_density_m3 for c in clouds]);od=np.array([c.optical_depth(sigma_laser)[0] for c in clouds]);reabs=np.array([scattering.reabsorption_probability(c) for c in clouds])
    radius=np.geomspace(1e-6,5e-3,150); reference=clouds[70];repulsion=scattering.repulsive_force(radius,reference);shadow=np.array([scattering.attenuated_pair_force(1e-20,reference) for _ in radius])
    time=np.linspace(0,20,201);common=dict(loading_rate_s=2e7,background_loss_s=.08,hot_rb_loss_s=.04,temperature_k=temperature,restoring_coefficient_n_m=kappa,atom_mass_kg=atom.mass_kg,scattering=scattering)
    independent=CollectiveLoading(**common,two_body_coefficient_m3_s=0)
    # This is an explicit user-specified sensitivity scenario, not a package
    # default or claimed 87Rb literature value.
    lossy=CollectiveLoading(**common,two_body_coefficient_m3_s=1e-16,beta_source="user-specified sensitivity scenario: 1e-16 m3/s")
    a,b=independent.evolve(time),lossy.evolve(time)
    vapor_scale=np.geomspace(.1,10,40);steady=[]
    for scale in vapor_scale:
        result=CollectiveLoading(**{**common,"loading_rate_s":2e7*scale},two_body_coefficient_m3_s=1e-16,beta_source="user-specified sensitivity scenario").evolve(np.linspace(0,100,301))
        steady.append((result["population"][-1],result["peak_density_m3"][-1],result["sigma_m"][-1,0]))
    steady=np.asarray(steady)
    fig,axes=plt.subplots(2,2,figsize=(12,9));axes[0,0].plot(time,a["population"],label="beta=0");axes[0,0].plot(time,b["population"],label="beta=1e-16 m3/s (user scenario)");axes[0,0].set(xlabel="time (s)",ylabel="N",title="Loading with separated one-/two-body losses");axes[0,0].legend()
    axes[0,1].loglog(populations,sigma*1e3);axes[0,1].set(xlabel="N",ylabel="RMS radius (mm)",title="Thermal to density-limited expansion")
    axes[1,0].loglog(populations,density,label="peak density");ax=axes[1,0].twinx();ax.loglog(populations,od,"C1",label="central OD");axes[1,0].set(xlabel="N",ylabel="peak n (m-3)");ax.set_ylabel("OD")
    axes[1,1].loglog(radius*1e3,repulsion,label="multiple scattering");axes[1,1].axhline(shadow[0],color="C1",label="Beer-Lambert shadow scale");axes[1,1].set(xlabel="radius (mm)",ylabel="force (N)",title="Mean-field force scales");axes[1,1].legend();fig.tight_layout();fig.savefig(OUT/"collective_mot_diagnostics.png",dpi=220);plt.close(fig)
    metadata={"simulation_version":__version__,"calculated":{"sigma0_m2":sigma0,"laser_cross_section_m2":sigma_laser,"saturation_total":saturation},"user_specified":{"temperature_k":temperature,"restoring_n_m":kappa,"loading_s":2e7,"background_loss_s":.08,"hot_rb_loss_s":.04,"beta_scenario_m3_s":1e-16,"reabsorption_cross_section_m2":"sigma0 upper reference"},"theory":"Walker, Sesko, and Wieman, Phys. Rev. Lett. 64, 408 (1990), DOI 10.1103/PhysRevLett.64.408; Coulomb-like multiple scattering and constant-density trend","limitations":"Gaussian/quasi-static size, effective scalar cross sections, single reabsorption, no exact radiative transport or coherent temperature"}
    np.savez_compressed(OUT/"collective_mot.npz",population_grid=populations,sigma_m=sigma,peak_density_m3=density,optical_depth=od,reabsorption_probability=reabs,radius_m=radius,multiple_scattering_force_n=repulsion,shadow_force_n=shadow,time_s=time,without_two_body_population=a["population"],with_two_body_population=b["population"],with_two_body_sigma_m=b["sigma_m"],vapor_density_scale=vapor_scale,steady_population_density_sigma=steady,metadata_json=json.dumps(metadata,sort_keys=True))
    with open(OUT/"parameter_provenance.csv","w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(("parameter","value","provenance"));writer.writerows((("sigma0",sigma0,"calculated atomic wavelength"),("sigma_laser",sigma_laser,"calculated two-level detuning/saturation"),("sigma_reabsorption",sigma0,"user scenario: resonant upper reference"),("beta",1e-16,"user sensitivity scenario; not default"),("temperature",temperature,"user supplied"),("kappa",kappa,"user supplied")))


if __name__=="__main__":main()
