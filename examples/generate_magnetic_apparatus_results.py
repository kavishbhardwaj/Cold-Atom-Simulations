"""Reference three-axis compensation-coil and imperfection diagnostics."""
from dataclasses import replace
import csv,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot.io.config import build_effective_model,build_multilevel_model,load_config
from cold_atom_mot.magnetic.coils import HelmholtzPair,ThreeAxisBiasCoils,AntiHelmholtzPair
from cold_atom_mot.magnetic.fields import CompositeField,ResidualField
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results"/"magnetic_apparatus"


def pair(axis,**kwargs): return HelmholtzPair.imperfect(axis,.04,.04,1,20,segments=128,**kwargs)


def centre_shift(field):
    gradient=np.diag([.1,.1,-.2]); return -np.linalg.solve(gradient,field)


def coherence(field):
    rate=build_multilevel_model(load_config(ROOT/"configs/rb87_d2_multilevel.yaml"))
    families=[replace(f,beam=replace(f.beam,coherence_group="fixed")) for f in rate.beam_families]
    obe=MultilevelOBE(rate.basis,families,ResidualField(uniform=field))
    rho=obe.steady_state_realization(); off=rho.copy();np.fill_diagonal(off,0)
    return np.linalg.norm(off)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    pairs=(pair([1,0,0]),pair([0,1,0]),pair([0,0,1]))
    # Representative measured-like nonorthogonal calibration, explicitly synthetic.
    matrix=np.array([[3.55,.08,-.04],[.05,3.48,.06],[-.03,.09,3.62]])*1e-5
    offset=np.array([2,-1,5])*1e-7;coils=ThreeAxisBiasCoils(pairs,matrix,offset)
    background=np.array([18,-7,46])*1e-6
    currents=coils.compensation_currents(background);residual=coils.calibrated_field(currents,background)
    grid=np.linspace(-.01,.01,31);x,y,z=np.meshgrid(grid,grid,[0.0],indexing="ij");points=np.stack([x[...,0],y[...,0],z[...,0]],axis=-1)
    physical=coils.physical_field(points,currents)+background
    fig,axes=plt.subplots(2,2,figsize=(10,8));values=(np.linalg.norm(physical,axis=-1),physical[...,0],physical[...,1],physical[...,2]);labels=("|B|","Bx","By","Bz")
    for ax,value,label in zip(axes.flat,values,labels):
        image=ax.pcolormesh(grid*1e3,grid*1e3,value*1e6,shading="auto");fig.colorbar(image,ax=ax,label="uT");ax.set(xlabel="x (mm)",ylabel="y (mm)",title=label)
    fig.tight_layout();fig.savefig(OUT/"compensated_field_maps.png",dpi=220);plt.close(fig)
    # Imperfections are controlled examples, not tolerance claims.
    cases=[("ideal",0,np.zeros(3)),("1% current imbalance",.01,pair([0,0,1],current_imbalance=.01).field([0,0,0])-pairs[2].field([0,0,0])),
           ("1 mm lateral displacement",.001,pair([0,0,1],lateral_displacement=[.001,0,0]).field([0,0,0])-pairs[2].field([0,0,0])),
           ("1 degree tilt",1,pair([0,0,1],tilt=[np.deg2rad(1),0,0]).field([0,0,0])-pairs[2].field([0,0,0])),
           ("10 mG background",10,np.array([0,0,1e-6])),
           ("50 Hz snapshot amplitude",10,np.array([1e-6,0,0]))]
    zero_coherence=coherence(np.zeros(3));rows=[]
    for name,magnitude,field in cases:
        shift=centre_shift(field); rows.append((name,magnitude,*field,*shift,
            coherence(field)/zero_coherence-1,"not calculated: no matched OBE diffusion"))
    with open(OUT/"imperfection_table.csv","w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(("imperfection","configured_magnitude","Bx_T","By_T","Bz_T","mot_dx_m","mot_dy_m","mot_dz_m","OBE_coherence_relative_change","PGC_temperature"));writer.writerows(rows)
    anti=AntiHelmholtzPair.symmetric(.04,.04,2,20,segments=128)
    imperfections=np.linspace(-.02,.02,21);zero=[]
    for imbalance in imperfections:
        candidate=AntiHelmholtzPair.symmetric(.04,.04,2,20,segments=128,current_imbalance=imbalance,lateral_offset=5e-4*abs(imbalance)/.02)
        zero.append(candidate.field_zero())
    zero=np.asarray(zero);plt.figure(figsize=(7,4.5));plt.plot(imbalances:=imperfections*100,zero*1e3);plt.xlabel("Anti-Helmholtz current imbalance (%)");plt.ylabel("Field-zero coordinate (mm)");plt.legend(("x","y","z"));plt.tight_layout();plt.savefig(OUT/"field_zero_imperfections.png",dpi=220);plt.close()
    jac=anti.jacobian(np.zeros(3)); divergence=np.trace(jac)
    metadata={"calibration_matrix_t_per_a":matrix.tolist(),"offset_t":offset.tolist(),"background_t":background.tolist(),"compensation_currents_a":currents.tolist(),"residual_t":residual.tolist(),
              "divergence_t_per_m":divergence,"fidelity":"Biot-Savart coil maps; effective linear MOT-centre shift; static 24-state OBE coherence diagnostic",
              "warning":"No PGC friction/diffusion or 10% tolerance claimed; requires moving OBE plus matched force-noise calculation."}
    np.savez_compressed(OUT/"magnetic_apparatus.npz",grid_m=grid,field_t=physical,calibration_matrix_t_a=matrix,compensation_currents_a=currents,residual_t=residual,antihelmholtz_jacobian_t_m=jac,imbalance=imperfections,field_zero_m=zero,metadata_json=json.dumps(metadata,sort_keys=True))


if __name__=="__main__":main()
