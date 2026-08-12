"""Vector-field OBE diagnostics for the configured 87Rb D2 light field.

This deliberately does not manufacture a sub-Doppler temperature: the present
24-state OBE point solver cannot yet average force and force-noise over the many
spatial lattice periods needed for a converged PGC friction/diffusion pair.
"""
from dataclasses import replace
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import physical_constants, hbar

from cold_atom_mot import __version__
from cold_atom_mot.io.config import build_multilevel_model, load_config
from cold_atom_mot.magnetic.fields import ResidualField
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results"/"polarization_gradient"
MU_B = physical_constants["Bohr magneton"][0]
PHASES = [0, 0, 0, np.pi/2, 0, np.pi/4]


def solver(field):
    rate = build_multilevel_model(load_config(ROOT/"configs/rb87_d2_multilevel.yaml"))
    families = []
    for index, family in enumerate(rate.beam_families):
        phase = PHASES[index] if index < 6 else 0.0
        # One explicitly phase-resolved realization. Cooling and repump address
        # separate rotating F blocks; this scan does not claim phase averaging.
        beam = replace(family.beam, phase=phase, coherence_group="fixed-realization")
        families.append(replace(family, beam=beam))
    return MultilevelOBE(rate.basis, families, ResidualField(uniform=np.asarray(field)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gauss = np.array([0, .0001, .0003, .001, .003, .01, .03, .1, .3, 1.0])
    tesla = gauss*1e-4
    directions = np.eye(3)
    force = np.empty((3, len(tesla), 3)); coherence = np.empty((3, len(tesla)))
    excited = np.empty((3, len(tesla)))
    for axis, direction in enumerate(directions):
        for index, magnitude in enumerate(tesla):
            obe = solver(magnitude*direction)
            rho = obe.steady_state_realization()
            ng = len(obe.basis.ground)
            force[axis, index] = obe.force(np.zeros(3), np.zeros(3), rho=rho)
            off = rho.copy(); np.fill_diagonal(off, 0)
            coherence[axis, index] = np.linalg.norm(off)
            excited[axis, index] = np.trace(rho[ng:, ng:]).real
    line = solver(np.zeros(3)).basis.line
    larmor = abs(line.lande_gf("ground", 2))*MU_B*tesla/hbar/(2*np.pi)
    pumping = line.gamma_rad_s*.08/(2*(1+(2*3)**2))/(2*np.pi)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for axis, label in enumerate(("Bx", "By", "Bz")):
        axes[0].semilogx(np.maximum(gauss*1e3, .03), coherence[axis], "o-", label=label)
        axes[1].semilogx(np.maximum(gauss*1e3, .03), force[axis,:,0]*1e21, "o-", label=label)
    for ax in axes:
        for marker in (1, 10, 100, 500): ax.axvline(marker, color="grey", alpha=.2)
        ax.set_xlabel("|B| (mG; zero displayed at 0.03 mG)"); ax.legend()
    axes[0].set_ylabel("Density-matrix off-diagonal norm")
    axes[0].set_title("Full-vector 24-state OBE, fixed phase realization")
    axes[1].set_ylabel("Instantaneous Fx at v=0 ($10^{-21}$ N)")
    axes[1].set_title("Orientation-dependent force offset (not friction)")
    fig.tight_layout(); fig.savefig(OUT/"vector_residual_obe.png", dpi=220); plt.close(fig)
    metadata = {"simulation_version": __version__, "fidelity": "24-state full-vector OBE; fixed optical-phase realization; v=0 point states",
                "recipe": "rb87_d2_multilevel.yaml powers/detunings; cooling phases [0,0,0,pi/2,0,pi/4]",
                "earth_reference_mg": 500,
                "limitation": "No beta or temperature threshold: moving spatial-period convergence and OBE-consistent force-noise diffusion are not established.",
                "literature_context": "Dalibard and Cohen-Tannoudji, JOSA B 6, 2023 (1989), DOI 10.1364/JOSAB.6.002023; qualitative magnetic disruption only, unmatched geometry."}
    np.savez_compressed(OUT/"vector_residual_obe.npz", field_g=gauss, field_t=tesla,
                        directions=directions, force_n=force, coherence_norm=coherence,
                        excited_population=excited, larmor_hz=larmor,
                        weak_drive_pumping_scale_hz=pumping,
                        metadata_json=json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
