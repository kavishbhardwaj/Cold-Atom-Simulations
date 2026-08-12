# Part V — From framework to digital twin: learning and reproduction

# 21. How a real apparatus would turn this framework into a digital twin

To make the simulation quantitatively predictive for one laboratory, the model would ingest measured quantities such as:

- six beam powers, waists, positions and pointing vectors;
- Jones polarization after the vacuum windows;
- cooling/repump spectra and AOM frequency offsets;
- three-axis magnetic calibration matrix;
- measured residual DC and 50/60-Hz fields;
- MOT-coil switch-off and eddy-current waveform;
- cell/reservoir temperatures and Rb pressure;
- background-gas pressure/composition;
- measured one- and two-body loss coefficients;
- chamber/capture geometry;
- measured sequence timing.

The same code hierarchy could then propagate those measured imperfections into force, capture, loading and eventually temperature predictions. This is the distinction between a **physics framework** and a **calibrated digital twin**.

---
# 22. Suggested learning path for a student

1. Run `configs/rb87_d2_mot.yaml` and understand the effective force map.
2. Derive the two-beam Doppler force by hand and compare with the code.
3. Inspect the generated 24-state basis and selection rules.
4. Run the multilevel rate equation and watch repump transfer population back into F=2.
5. Study the two-level OBE and reproduce the analytical \(\rho_{ee}\).
6. Inspect the vector Zeeman spectrum and compare linear vs exact shifts.
7. Run the reduced PGC example and follow one Sisyphus optical period.
8. Compare zero and finite residual fields, remembering that the current 9.4 mG marker is a timescale comparison, not a temperature threshold.
9. Launch thermal atoms through the capture calculation and inspect why only the slow tail contributes.
10. Add the loading equation and then the optional collective cloud.
11. Only after understanding these pieces, move to the 24-state moving OBE and external validation.

---
# 23. Repository map for this tutorial

- `src/cold_atom_mot/atomic/` — isotope/line data, hyperfine basis, vector Zeeman physics.
- `src/cold_atom_mot/laser/` — Gaussian beams, polarization, Jones optics, coherence groups and apparatus topology.
- `src/cold_atom_mot/magnetic/` — quadrupoles, residual fields, Helmholtz/anti-Helmholtz coils.
- `src/cold_atom_mot/physics/` — effective force, rate equations, OBEs, PGC and collective mean field.
- `src/cold_atom_mot/solvers/` — deterministic and photon-event trajectories.
- `src/cold_atom_mot/simulation/` — capture criteria, ensembles and experimental sequence.
- `src/cold_atom_mot/vacuum.py` — vapour pressure, thermal flux, collision/loss and loading equations.
- `configs/` — explicit physical and numerical assumptions for reproducible runs.
- `results/` — figures, numerical arrays and metadata.
- `docs/validation.md` — what is externally verified and what remains pending.

---
# 24. Reproducing the repository from a clean checkout

This section turns the tutorial into an executable workflow. The repository uses **Python 3.10 or newer** and keeps physical calculations in SI units. The package metadata declares NumPy, SciPy, Matplotlib, PyYAML and SymPy as core dependencies; QuTiP and PyLCP are optional pinned validation dependencies.

## 24.1 Create an isolated environment

From a clean clone of the repository:

```bash
git clone https://github.com/kavishbhardwaj/cold-atom-simulations.git
cd cold-atom-simulations

python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1

python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
```

`requirements-dev.txt` installs the package in editable mode together with the test runner. The optional independent-validation environment is:

```bash
python -m pip install -r requirements-validation.txt
```

At the time this tutorial was written, the validation file pins QuTiP 5.3.1 and PyLCP 1.0.2 so that external comparisons use a known software environment.

## 24.2 Run the five main physics entry points

The following commands reproduce the core solver paths described in this tutorial:

```bash
python -m cold_atom_mot simulate configs/rb87_d2_mot.yaml
python -m cold_atom_mot rate-equation configs/rb87_d2_multilevel.yaml
python -m cold_atom_mot obe configs/rb87_d2_two_level_obe.yaml
python -m cold_atom_mot subdoppler configs/rb87_d2_polarization_gradient.yaml
python -m cold_atom_mot loading configs/rb_vapor_loading.yaml
```

These are not interchangeable commands. They deliberately invoke different fidelity levels:

- `simulate`: effective semiclassical MOT force and trajectories;
- `rate-equation`: 24-state population dynamics without coherence;
- `obe`: the transparent two-level coherent benchmark;
- `subdoppler`: the reduced phase-resolved Sisyphus population model;
- `loading`: thermal vapour flux, trajectory-derived capture and loading/loss.

The 24-state moving OBE is primarily a research-level backend used by diagnostics and validation scripts rather than a large-grid default CLI calculation.

## 24.3 Regenerate representative committed figures

The repository stores generated figures together with numerical arrays/metadata. Representative regeneration commands are:

```bash
python examples/generate_vector_zeeman_results.py
python examples/generate_six_beam_apparatus_results.py
python examples/generate_magnetic_apparatus_results.py
python examples/generate_sequence_results.py
python examples/generate_capture_loading_results.py
python examples/generate_collective_mot_results.py
```

For the independent software checks:

```bash
python examples/generate_external_validation_results.py
```

The corresponding outputs are written below `results/`. The result archives store the arrays and metadata used to construct the figures; the figure itself should never be the only record of a calculation.

## 24.4 Run the validation suite

```bash
python -m pytest -q
python -m compileall src tests examples
git diff --check
```

A reproduction should not be judged only by whether a plot looks similar. The following physics checks are more informative:

1. the two-level OBE should reproduce its analytical steady-state population;
2. the matched QuTiP result should agree at numerical precision for the tested cases;
3. the normalized two-beam force should reproduce the matched PyLCP curve;
4. the generated 87Rb Zeeman spectrum should agree with the matched PyLCP spectrum at the reported sub-Hz scale;
5. probability generators must conserve total population;
6. density matrices must remain trace one, Hermitian and positive within the numerical tolerance;
7. the ideal quadrupole field must satisfy `trace(gradient)=0`;
8. stochastic runs must reproduce when the stored random seed is reused;
9. research calculations must be refined in time step, integration duration, phase sampling and ensemble size rather than trusting one default setting.

## 24.5 Reference configurations used throughout the tutorial

The main effective-MOT reference is stored in `configs/rb87_d2_mot.yaml`:

| Parameter | Reference value | Meaning |
|---|---:|---|
| isotope/line | 87Rb D2 | reference atom |
| power per cooling beam | 10 mW | six independent beams |
| 1/e^2 waist | 8 mm | Gaussian transverse radius |
| cooling detuning | -2 Gamma | effective MOT reference |
| radial quadrupole gradient | 0.10 T/m | 10 G/cm |
| gravity | 9.80665 m/s^2 | along -z |

The multilevel configuration `configs/rb87_d2_multilevel.yaml` uses the same cooling-beam reference and adds a 0.5 mW-per-beam repump tuned near `F=1 -> F'=2`.

The reduced PGC configuration `configs/rb87_d2_polarization_gradient.yaml` uses `Delta=-3 Gamma` and `s=0.08` per beam. Its committed phase/coherence choice is a controlled model configuration, not a claim that every real MOT has the same phase topology.

The experimental timeline in `configs/rb87_d2_reference_sequence.yaml` is explicitly illustrative: its ramps, eddy-current time constant and residual field are reproducible inputs, not a fitted experimental optimum.

The vapour-loading configuration `configs/rb_vapor_loading.yaml` keeps reservoir, vapour and background-gas temperatures distinct, uses natural-isotope fractions, and exposes all collision/loss parameters instead of inventing them.

## 24.6 Numerical reproducibility versus physical reproducibility

These are different concepts.

**Numerical reproducibility** means the same code, configuration, package versions, seed and tolerances reproduce the stored arrays within expected floating-point variation.

**Physical reproducibility** means another apparatus with the same measured beam profiles, magnetic fields, pressures, loss coefficients and timing produces the same physical observable. The repository currently targets the first strongly and provides the interfaces required for the second, but it is not automatically calibrated to a laboratory.

Therefore every serious comparison should record:

- commit SHA;
- configuration file;
- Python/package versions;
- random seed;
- solver mode and tolerances;
- integration duration and step controls;
- phase-sampling settings for incoherent coherent-OBE calculations;
- ensemble size for stochastic/capture calculations;
- which parameters are calculated, literature-sourced, user supplied or experimentally calibrated.

---
# 25. Equation-to-code map

This table is intended to let a student move directly from a derivation in this tutorial to the implementation that evaluates it.

| Physics/equation | Main implementation | Configuration/result entry point |
|---|---|---|
| atomic constants, hyperfine energies, Wigner/CG transition strengths | `src/cold_atom_mot/atomic/species.py` | atomic-structure generators |
| vector hyperfine-Zeeman Hamiltonian | `src/cold_atom_mot/atomic/zeeman.py` | `generate_vector_zeeman_results.py` |
| Gaussian beam intensity and propagation | `src/cold_atom_mot/laser/beam.py` | six-beam apparatus generator |
| Jones optics and spherical polarization fractions | `src/cold_atom_mot/laser/polarization.py` | apparatus/PGC calculations |
| ideal/residual magnetic fields | `src/cold_atom_mot/magnetic/fields.py` | MOT and residual-field configs |
| physical circular/Helmholtz/anti-Helmholtz coils | `src/cold_atom_mot/magnetic/coils.py` | magnetic-apparatus generator |
| effective scattering force | `src/cold_atom_mot/physics/force.py` | `rb87_d2_mot.yaml` |
| multilevel population rate equations | `src/cold_atom_mot/physics/multilevel.py` or corresponding rate-equation module | `rb87_d2_multilevel.yaml` |
| two-level OBE/Lindblad benchmark | `src/cold_atom_mot/physics/optical_bloch.py` | two-level OBE config/validation |
| 24-state moving OBE and Hamiltonian-gradient force | `src/cold_atom_mot/physics/multilevel_obe.py` | vector-field and research diagnostics |
| reduced Sisyphus/PGC model | `src/cold_atom_mot/physics/subdoppler.py` | `rb87_d2_polarization_gradient.yaml` |
| Newton/RK45 trajectories | `src/cold_atom_mot/solvers/` | effective/capture runs |
| thermal vapour flux and loading/loss equations | `src/cold_atom_mot/vacuum.py` | `rb_vapor_loading.yaml` |
| capture classification and rare-event estimator | `src/cold_atom_mot/simulation/capture.py` | capture/loading generator |
| experimental control timeline | `src/cold_atom_mot/simulation/sequence.py` | `rb87_d2_reference_sequence.yaml` |
| collective Gaussian MOT | `src/cold_atom_mot/physics/collective.py` | collective-MOT generator |

If a path changes in a later refactor, the governing principle is more important than the exact filename: configurations specify physical assumptions, `physics/` contains the model equations, `solvers/` propagates dynamics, and `results/` stores generated evidence plus metadata.

---
# 26. A reproducibility checklist before accepting a new result

Before interpreting any new curve as physics, ask:

1. **What observable is being calculated?** Force, population, trajectory, capture probability, loading rate, friction or temperature are not interchangeable.
2. **Which model generated it?** Effective, rate-equation, reduced PGC or full OBE?
3. **Which terms are omitted by that model?** State the omission next to the result.
4. **Are all units and detuning conventions explicit?** In particular distinguish Hz from rad/s and `Gamma` from `Gamma/2pi`.
5. **Has the numerical answer converged?** Refine time step, tolerances, integration window, phase samples or ensemble size as appropriate.
6. **Has the relevant simpler limit been recovered?** Examples: Lorentzian scattering, analytical two-level OBE, Maxwell-consistent gradient, zero-density collective limit.
7. **Is there an independent benchmark?** QuTiP, PyLCP, analytical theory, or a fully specified experiment.
8. **Are experimental parameters actually measured?** If not, label them as scenario inputs.
9. **Is uncertainty represented?** Rare-event zero counts need confidence bounds; one stochastic seed is not an error bar.
10. **Does the claimed conclusion require physics absent from the model?** The present PGC diffusion issue is the canonical example: a friction curve alone is not a defensible sub-Doppler temperature.

This checklist is the practical implementation of the repository's central philosophy: make approximations explicit before making conclusions.

---
# 27. Final perspective

The repository evolved by repeatedly asking a simple question:

> **What is the least complicated model that can answer the present physical question without hiding an approximation?**

For force maps and capture, the effective MOT is often sufficient. For optical pumping, use rate equations. For vector magnetic mixing and coherences, use the OBE. For Sisyphus physics, resolve polarization gradients. For atom number, add thermal flux and explicit loss inputs. For large clouds, add collective effects only after the single-atom model is understood.

The most important scientific decision was therefore not one particular equation. It was to preserve a **model hierarchy with explicit validity boundaries**, and to validate simpler pieces independently before treating a more complicated calculation as predictive.

That makes the repository useful both as a tutorial and as a foundation for a future experiment-specific cold-atom digital twin.
