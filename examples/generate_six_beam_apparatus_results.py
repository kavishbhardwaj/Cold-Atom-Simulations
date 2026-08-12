"""Physical six-beam apparatus and controlled imperfection diagnostics."""
from dataclasses import replace
import csv,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root

from cold_atom_mot.io.config import build_effective_model,load_config
from cold_atom_mot.laser.apparatus import Retroreflection,SixBeamApparatus
from cold_atom_mot.laser.polarization import JonesElement

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results"/"laser_apparatus"


def center_and_coefficients(model):
    center=root(lambda x:model.force(x,np.zeros(3)),np.zeros(3)).x
    step=1e-6; restoring=np.empty(3); damping=np.empty(3)
    for axis in range(3):
        dx=np.zeros(3);dx[axis]=step;dv=np.zeros(3);dv[axis]=1e-3
        restoring[axis]=-(model.force(center+dx,np.zeros(3))[axis]-model.force(center-dx,np.zeros(3))[axis])/(2*step)
        damping[axis]=-(model.force(center,dv)[axis]-model.force(center,-dv)[axis])/.002
    return center,restoring,damping


def main():
    OUT.mkdir(parents=True,exist_ok=True);model=build_effective_model(load_config(ROOT/"configs/rb87_d2_mot.yaml"))
    # Attach explicit Jones trains equivalent to the ideal circular constructor.
    beams=[]
    for beam in model.beams:
        elements=(JonesElement("quarter_wave",np.pi/4*beam.helicity),)
        beams.append(replace(beam,jones_vector=np.array([1,0]),optical_elements=elements))
    apparatus=SixBeamApparatus(tuple(beams));model.beams=list(apparatus.beams)
    fractions=apparatus.local_polarizations([0,0,1]);rows=[]
    for beam,local in zip(apparatus.beams,fractions):
        rows.append((beam.label,beam.power,beam.waist,beam.waist_y,*beam.direction,
                     local[-1],local[0],local[1],beam.coherence_group))
    with open(OUT/"beam_diagnostics.csv","w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(("label","power_W","wx_m","wy_m","kx_hat","ky_hat","kz_hat","sigma_minus","pi","sigma_plus","coherence_group"));writer.writerows(rows)
    fig=plt.figure(figsize=(7,7));ax=fig.add_subplot(111,projection="3d");colors=plt.cm.tab10(np.arange(6))
    for beam,color in zip(apparatus.beams,colors):
        start=-.012*beam.direction;ax.quiver(*start,*beam.direction,length=.009,color=color,label=f"{beam.label}: {beam.power*1e3:.1f} mW")
    ax.set(xlim=(-.015,.015),ylim=(-.015,.015),zlim=(-.015,.015),xlabel="x (m)",ylabel="y (m)",zlabel="z (m)",title="Three counterpropagating physical beam pairs");ax.legend(fontsize=7);fig.tight_layout();fig.savefig(OUT/"six_beam_apparatus.png",dpi=220);plt.close(fig)
    imbalances=np.array([-10,-5,-1,0,1,5,10]);centers=[];restoring=[];damping=[]
    for percent in imbalances:
        trial=replace(apparatus.beams[0],power=apparatus.beams[0].power*(1+percent/100));model.beams=[trial,*apparatus.beams[1:]]
        c,k,b=center_and_coefficients(model);centers.append(c);restoring.append(k);damping.append(b)
    errors=np.linspace(-8,8,33);qwp_force=[]
    for degrees in errors:
        qwp=JonesElement("quarter_wave",np.pi/4+np.deg2rad(degrees));model.beams=[replace(apparatus.beams[0],jones_vector=[1,0],optical_elements=(qwp,)),*apparatus.beams[1:]]
        qwp_force.append(model.force([1e-4,0,0],np.zeros(3))[0])
    pointing=np.linspace(0,5,21);point_center=[]
    for mrad in pointing:
        direction=apparatus.beams[0].direction+np.array([0,mrad*1e-3,0]);model.beams=[replace(apparatus.beams[0],direction=direction),*apparatus.beams[1:]]
        point_center.append(center_and_coefficients(model)[0])
    fig,axes=plt.subplots(1,3,figsize=(15,4.5));axes[0].plot(imbalances,np.array(centers)[:,0]*1e3,"o-");axes[0].set(xlabel="x-beam power error (%)",ylabel="MOT centre x (mm)")
    axes[1].plot(errors,np.array(qwp_force)*1e21);axes[1].set(xlabel="QWP angle error (deg)",ylabel="Fx at x=0.1 mm ($10^{-21}$ N)")
    axes[2].plot(pointing,np.linalg.norm(point_center,axis=1)*1e3);axes[2].set(xlabel="Pointing error (mrad)",ylabel="MOT centre displacement (mm)")
    fig.tight_layout();fig.savefig(OUT/"apparatus_imperfections.png",dpi=220);plt.close(fig)
    retro=SixBeamApparatus.three_retroreflected(apparatus.beams[::2],tuple(Retroreflection(.9,coherence_group=f"pair{i}") for i in range(3)))
    metadata={"fidelity":"effective MOT force apparatus sensitivity; not coherent PGC/capture","quantization_axis":[0,0,1],"retroreflection_power_transmission":.9,
              "warning":"Individual sweeps hold all other parameters fixed; material effects are recipe-specific."}
    np.savez_compressed(OUT/"six_beam_apparatus.npz",imbalance_percent=imbalances,center_m=centers,restoring_n_m=restoring,damping_kg_s=damping,qwp_error_deg=errors,qwp_force_n=qwp_force,pointing_mrad=pointing,pointing_center_m=point_center,
                        retro_directions=np.array([b.direction for b in retro.beams]),metadata_json=json.dumps(metadata,sort_keys=True))


if __name__=="__main__":main()
