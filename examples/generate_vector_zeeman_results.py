"""Generate exact vector hyperfine-Zeeman reference diagnostics."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cold_atom_mot import __version__
from cold_atom_mot.atomic.species import get_atomic_line
from cold_atom_mot.atomic.zeeman import hyperfine_zeeman_hamiltonian, linear_zeeman_energies

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "atomic_structure"


def save(name):
    plt.tight_layout()
    plt.savefig(OUTPUT/f"{name}.png", dpi=220, bbox_inches="tight")
    plt.savefig(OUTPUT/f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    plt.close()


def spectrum(line, manifold, fields):
    exact = np.array([np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line, manifold, [0,0,b]))
                      for b in fields])/(2*np.pi)
    linear = np.array([np.sort(linear_zeeman_energies(line, manifold, [0,0,b]))
                       for b in fields])/(2*np.pi)
    centre = exact[len(fields)//2].mean()
    return exact-centre, linear-centre


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams["svg.hashsalt"] = "cold-atom-mot-vector-zeeman"
    line = get_atomic_line("87Rb", "D2")
    ground_field = np.linspace(-0.15, 0.15, 241)
    excited_field = np.linspace(-0.03, 0.03, 241)
    ground_exact, ground_linear = spectrum(line, "ground", ground_field)
    excited_exact, excited_linear = spectrum(line, "excited", excited_field)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(ground_field*1e3, ground_exact/1e9, color="tab:blue", lw=.8)
    axes[0].set(xlabel="Magnetic field Bz (mT)", ylabel="Frequency from centroid (GHz)",
                title="87Rb 5S1/2 exact hyperfine–Zeeman spectrum")
    axes[1].plot(excited_field*1e3, excited_exact/1e6, color="tab:orange", lw=.7)
    axes[1].set(xlabel="Magnetic field Bz (mT)", ylabel="Frequency from centroid (MHz)",
                title="87Rb 5P3/2 exact hyperfine–Zeeman spectrum")
    save("exact_zeeman_spectra")

    positive = ground_field >= 0
    deviation = np.max(np.abs(ground_exact-ground_linear), axis=1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for index in range(ground_exact.shape[1]):
        axes[0].plot(ground_field[positive]*1e3, ground_exact[positive,index]/1e9,
                     color="tab:blue", lw=.7)
        axes[0].plot(ground_field[positive]*1e3, ground_linear[positive,index]/1e9,
                     color="black", ls="--", lw=.45, alpha=.55)
    axes[0].set(xlabel="B (mT)", ylabel="Frequency from centroid (GHz)",
                title="Exact (blue) versus linear gF mF (dashed)")
    axes[1].loglog(np.maximum(np.abs(ground_field), 1e-12)*1e3,
                   np.maximum(deviation, 1e-12)/1e3)
    axes[1].set(xlabel="|B| (mT)", ylabel="Largest spectral difference (kHz)",
                title="Nonlinear-Zeeman/hyperfine-mixing correction")
    save("linear_vs_exact_zeeman")

    magnitudes = np.geomspace(1e-10, 3e-2, 120)
    direction_error = []
    for magnitude in magnitudes:
        reference = np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line, "excited", [0,0,magnitude]))
        tilted = np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(
            line, "excited", magnitude*np.array([1,2,-3])/np.sqrt(14)))
        direction_error.append(np.max(np.abs(reference-tilted))/(2*np.pi))
    plt.figure(figsize=(7.2, 4.5))
    plt.loglog(magnitudes*1e3, np.maximum(direction_error, 1e-9))
    plt.xlabel("Field magnitude (mT)"); plt.ylabel("Maximum directional spectral difference (Hz)")
    plt.title("Rotation covariance: z-directed versus tilted field")
    save("zeeman_direction_covariance")

    metadata = {
        "simulation_version": __version__, "isotope": "87Rb", "line": "D2",
        "hamiltonian": "A I.J + electric quadrupole + mu_B(g_J J+g_I I).B",
        "basis": "uncoupled |m_I,m_J>, displayed through coupled |F,m_F> transform",
        "units": {"field": "T", "frequency": "Hz"},
        "limitations": "fine-structure mixing and diamagnetic terms are omitted",
    }
    np.savez_compressed(OUTPUT/"vector_zeeman_reference.npz",
                        ground_field_t=ground_field, ground_exact_hz=ground_exact,
                        ground_linear_hz=ground_linear, excited_field_t=excited_field,
                        excited_exact_hz=excited_exact, excited_linear_hz=excited_linear,
                        direction_field_t=magnitudes,
                        direction_spectral_error_hz=direction_error,
                        metadata_json=json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
