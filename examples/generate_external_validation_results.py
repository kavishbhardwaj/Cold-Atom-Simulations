"""Reproducible QuTiP, PyLCP, analytical, and literature audit matrix."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pylcp,qutip

from cold_atom_mot.physics.optical_bloch import TwoLevelOBE
from cold_atom_mot.atomic.species import get_atomic_line,MU_B
from cold_atom_mot.atomic.zeeman import hyperfine_zeeman_hamiltonian
from scipy.constants import hbar

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results"/"validation"


def qutip_benchmark():
    gamma=1.7;saturations=np.array([.01,.05,.2,1,5]);detunings=np.array([-3,-1,0,1])
    internal=[];external=[];analytic=[]
    for s in saturations:
      for d in detunings:
        model=TwoLevelOBE.from_saturation(gamma,d*gamma,s);omega=model.rabi_frequency
        h=qutip.Qobj([[0,np.conjugate(omega)/2],[omega/2,-model.detuning]])
        c=np.sqrt(gamma)*qutip.basis(2,0)*qutip.basis(2,1).dag()
        internal.append(model.steady_state()[1,1].real);external.append(qutip.steadystate(h,[c]).full()[1,1].real);analytic.append(model.analytic_excited_population())
    return np.asarray(internal),np.asarray(external),np.asarray(analytic),saturations,detunings


def pylcp_benchmark():
    saturation=.05;detuning=-2.;velocities=np.linspace(-.5,.5,101)
    params=[dict(kvec=np.array([1.,0,0]),pol=1,s=saturation,delta=detuning),dict(kvec=np.array([-1.,0,0]),pol=1,s=saturation,delta=detuning)]
    equation=pylcp.heuristiceq(pylcp.laserBeams(params,beam_type=pylcp.infinitePlaneWaveBeam),pylcp.constantMagneticField([0,0,0]),gamma=1,k=1,mass=1)
    external=np.array([equation.force(np.zeros(3),[v,0,0],0)[0][0] for v in velocities])
    internal=np.array([.5*saturation/(1+2*saturation+(2*(detuning-v))**2)-.5*saturation/(1+2*saturation+(2*(detuning+v))**2) for v in velocities])
    return velocities,internal,external


def main():
    OUT.mkdir(parents=True,exist_ok=True);internal,external,analytic,saturations,detunings=qutip_benchmark();v,ours,pylcp_force=pylcp_benchmark()
    qutip_abs=np.max(np.abs(internal-external));qutip_rel=np.max(np.abs(internal-external)/np.maximum(np.abs(external),1e-30));analytic_abs=np.max(np.abs(internal-analytic))
    pylcp_rms=np.sqrt(np.mean((ours-pylcp_force)**2));pylcp_rel=np.max(np.abs(ours-pylcp_force)/np.maximum(np.abs(pylcp_force),1e-15))
    line=get_atomic_line("87Rb","D2");h0,mu=pylcp.hamiltonians.hyperfine_coupled(.5,1.5,line.ground_g_j,-line.species.nuclear_g_factor,line.species.ground_hyperfine_a_hz,muB=MU_B/hbar/(2*np.pi));fields=np.array([0,1e-6,1e-4]);zeeman_internal=[];zeeman_pylcp=[]
    for field in fields:
        zeeman_internal.append(np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line,"ground",[0,0,field]))/(2*np.pi));zeeman_pylcp.append(np.linalg.eigvalsh(h0+mu[1]*field))
    zeeman_internal=np.asarray(zeeman_internal);zeeman_pylcp=np.asarray(zeeman_pylcp);zeeman_error=np.max(np.abs(zeeman_internal-zeeman_pylcp))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));axes[0].plot(external,internal,"o");limit=max(external);axes[0].plot([0,limit],[0,limit],"k--");axes[0].set(xlabel="QuTiP rho_ee",ylabel="Internal rho_ee",title="Matched two-level steady states")
    axes[1].plot(v,ours,label="internal formula");axes[1].plot(v,pylcp_force,"--",label="PyLCP public API");axes[1].set(xlabel="v (Gamma/k)",ylabel="F (hbar k Gamma)",title="Matched 1D two-beam molasses");axes[1].legend();fig.tight_layout();fig.savefig(OUT/"independent_software_comparison.png",dpi=220);plt.close(fig)
    matrix=[
      {"claim":"Two-level OBE steady state, Liouvillian, Rabi/decay","A":"ANALYTICALLY VERIFIED","B":"INTERNAL TESTED","C":"INDEPENDENT-SOFTWARE VERIFIED","D":"not needed","E":"standard OBE","F":"NOT YET VALIDATED"},
      {"claim":"1D low-saturation two-beam force","A":"ANALYTICALLY VERIFIED","B":"INTERNAL TESTED","C":"not needed","D":"INDEPENDENT-SOFTWARE VERIFIED","E":"standard Doppler theory","F":"NOT YET VALIDATED"},
      {"claim":"24-state moving multilevel OBE","A":"trace/Hermiticity limits","B":"INTERNAL TESTED","C":"small-system only; full case NOT YET VALIDATED","D":"NOT YET VALIDATED","E":"NOT YET VALIDATED","F":"NOT YET VALIDATED"},
      {"claim":"Population PGC model","A":"controlled adiabatic model","B":"INTERNAL TESTED","C":"not compared","D":"not compared","E":"LITERATURE-TREND VERIFIED","F":"NOT YET VALIDATED"},
      {"claim":"MOT/loading predictions","A":"component formulas","B":"INTERNAL TESTED","C":"not compared","D":"simple 1D force only","E":"qualitative trends","F":"EXPERIMENTALLY COMPARED: parameters incomplete"}]
    metadata={"qutip_version":qutip.__version__,"pylcp_version":"1.0.2","conventions":{"gamma":1.7,"hamiltonian":"[[0,Omega*/2],[Omega/2,-Delta]]","collapse":"sqrt(Gamma)|g><e|","saturation":"2|Omega|^2/Gamma^2"},
      "errors":{"qutip_max_absolute_population":qutip_abs,"qutip_max_relative_population":qutip_rel,"analytic_max_absolute_population":analytic_abs,"pylcp_force_rms":pylcp_rms,"pylcp_force_max_relative":pylcp_rel,"pylcp_rb87_zeeman_max_absolute_hz":zeeman_error},
      "literature":[{"citation":"J. Dalibard and C. Cohen-Tannoudji, JOSA B 6, 2023 (1989)","doi":"10.1364/JOSAB.6.002023","status":"LITERATURE-TREND VERIFIED only","matched":"red detuning, polarization-gradient optical pumping/light shifts","unmatched":"paper geometries, reduced manifolds, diffusion and dimensionality"},{"citation":"P. D. Lett et al., JOSA B 6, 2084 (1989)","doi":"10.1364/JOSAB.6.002084","status":"EXPERIMENTALLY COMPARED, not quantitatively reproduced","matched":"sub-Doppler trend","unmatched":"87Rb requirement (paper is sodium), apparatus phases, intensity calibration"}],
      "warning":"No parameter was fitted. PyLCP comparison is its heuristic two-level public API, not a full 87Rb multilevel validation."}
    np.savez_compressed(OUT/"external_validation.npz",qutip_internal_population=internal,qutip_population=external,analytic_population=analytic,saturations=saturations,detunings_gamma=detunings,velocity_gamma_over_k=v,internal_force=ours,pylcp_force=pylcp_force,zeeman_field_t=fields,internal_zeeman_hz=zeeman_internal,pylcp_zeeman_hz=zeeman_pylcp,matrix_json=json.dumps(matrix),metadata_json=json.dumps(metadata,sort_keys=True))


if __name__=="__main__":main()
