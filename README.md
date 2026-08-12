# Cold-atom and magneto-optical-trap simulations

A reproducible SI-unit Python framework for **rubidium MOTs, laser cooling, and cold-atom apparatus modelling**. The default reference system is **87Rb D2**, with generated hyperfine/Zeeman structure, six physical laser beams, magnetic coils and residual fields, multilevel light-force models, optical Bloch equations, polarization-gradient cooling, trajectories, vapour loading, experimental sequences, and optional collective-cloud physics.

The goal is to connect atomic physics to laboratory-level questions while keeping every approximation visible. This repository is **not a calibrated digital twin of a particular experiment by default**: measured apparatus parameters and externally validated high-fidelity cooling models are still required for quantitative temperature, atom-number, and tolerance predictions.

> **New to the repository?** Start with the **[continuous equation-to-result tutorial](docs/tutorial/continuous_walkthrough.md)**. It defines symbols before they appear, derives the physics in calculation order, explains why each model/tool was chosen, states the approximation, shows the corresponding result or pedagogical plot, and then explains why the next model is needed. The [equation inventory](docs/tutorial/00_notation_and_equation_inventory.md) is the completeness/glossary companion and the [visual atlas](docs/tutorial/equation_visual_atlas.md) collects the main equation plots.

## What is trustworthy today?

**Implemented** means the model exists and has internal tests. **Externally verified** means a matched result agrees with an independent implementation. **Experimentally calibrated** would require measured apparatus inputs and quantitative agreement with a specified experiment.

| Capability | Current status |
|---|---|
| Two-level OBE, Rabi dynamics, spontaneous decay | **Analytically + QuTiP verified** |
| Normalized 1-D two-beam Doppler force | **PyLCP verified** |
| 87Rb ground hyperfine–Zeeman spectrum | **PyLCP verified** |
| 24-state moving 87Rb D2 OBE | **Implemented + internally tested; external multilevel validation pending** |
| Polarization-gradient cooling | **Mechanism/trends implemented; quantitative OBE-consistent temperature pending** |
| Six-beam optics, coils, sequence timing, capture/loading | **Implemented + internally tested; apparatus calibration dependent** |
| Collective Gaussian MOT / multiple scattering | **Optional mean-field model; literature-trend verified, not experimentally calibrated** |

Independent benchmarks use **no fitted parameters**. In the matched tests, the two-level steady-state population agrees with QuTiP to `5.55e-17`, the normalized two-beam force agrees with PyLCP to a maximum relative difference of about `7.9e-15`, and the 87Rb ground Zeeman spectrum agrees with PyLCP within about `0.57 Hz` over the tested fields. Full conventions and limits are in [validation](docs/validation.md).

## What is modelled?

| Layer | Purpose |
|---|---|
| Atomic structure | 85Rb/87Rb D1/D2 hyperfine and Zeeman bases, Wigner-generated transitions |
| Effective MOT | Fast Doppler/Zeeman force, trajectories, capture studies |
| Multilevel rate equations | Cooling + repump populations and optical pumping |
| 24-state moving OBE | Vector Zeeman Hamiltonian, coherences, beam-resolved Doppler shifts and forces |
| Polarization-gradient model | Phase-resolved light shifts, pumping and Sisyphus force |
| Six-beam apparatus | Independent beams, Gaussian propagation, Jones optics, QWP errors, retroreflection |
| Magnetic apparatus | Anti-Helmholtz/Helmholtz coils, three-axis compensation, gradients, AC fields, eddy currents |
| Experimental sequence | MOT load → CMOT → field switch-off → settling → PGC/molasses → TOF |
| Capture/loading | Thermal surface flux, trajectory-derived capture, loading/loss equations |
| Collective MOT | Optional Gaussian density, two-body loss, shadowing, multiple scattering, radiation-trapping proxy |

For 87Rb D2 the atomic basis contains the full **8 ground + 16 excited = 24 hyperfine-Zeeman states**. The moving OBE is implemented, but the repository deliberately does not label the complete multilevel force or PGC temperature as externally validated yet.

D1 structure is generated, but D1 gray molasses is **not** faked: Raman dark-state physics still needs a dedicated coherent implementation and validation.

## Selected results

| Independent software validation | Exact vector Zeeman structure |
|---|---|
| ![QuTiP and PyLCP validation](results/validation/independent_software_comparison.png) | ![Exact Zeeman spectra](results/atomic_structure/exact_zeeman_spectra.png) |

| Six physical MOT beams | Three-axis magnetic apparatus |
|---|---|
| ![Six-beam apparatus](results/laser_apparatus/six_beam_apparatus.png) | ![Magnetic field maps](results/magnetic_apparatus/compensated_field_maps.png) |

| Polarization-gradient force | Experimental sequence |
|---|---|
| ![Sub-Doppler force](results/polarization_gradient/subdoppler_force_velocity.png) | ![Sequence timeline](results/sequence/sequence_timeline.png) |

| Vapour capture/loading | Collective MOT mean field |
|---|---|
| ![Vapour loading](results/capture_loading/vapor_capture_loading.png) | ![Collective MOT diagnostics](results/collective_mot/collective_mot_diagnostics.png) |

See the **[scientific results gallery](results/README.md)** for captions, numerical data, provenance, held-fixed parameters, and fidelity limits.

## Repository guide

Start with the **[continuous tutorial](docs/tutorial/continuous_walkthrough.md)**; use the [tutorial index](docs/tutorial/README.md) and [documentation map](docs/README.md) for deeper reference pages. The main technical references are:

- [Notation/equation inventory](docs/tutorial/00_notation_and_equation_inventory.md) — symbols and governing equations.
- [Equation visual atlas](docs/tutorial/equation_visual_atlas.md) — equation plots and corresponding simulation results.
- [Model hierarchy](docs/model_hierarchy.md) — what each solver includes and neglects.
- [Validation](docs/validation.md) — independent checks, error metrics, and remaining validation gates.
- [Cooling physics](docs/cooling_physics.md) — force, optical pumping, recoil, and PGC assumptions.
- [Six-beam apparatus](docs/six_beam_apparatus.md) — real beam geometry, Jones optics, and imperfections.
- [Magnetic apparatus](docs/magnetic_apparatus.md) — compensation coils, stray fields, switching and eddy currents.
- [Experimental sequences](docs/experimental_sequences.md) — time-dependent laboratory cycle.
- [Collective MOT physics](docs/collective_mot.md) — optional density-dependent mean-field extension.

## Reproduce

```bash
python -m pip install -r requirements-dev.txt

python -m cold_atom_mot simulate configs/rb87_d2_mot.yaml
python -m cold_atom_mot rate-equation configs/rb87_d2_multilevel.yaml
python -m cold_atom_mot obe configs/rb87_d2_two_level_obe.yaml
python -m cold_atom_mot subdoppler configs/rb87_d2_polarization_gradient.yaml
python -m cold_atom_mot loading configs/rb_vapor_loading.yaml
```

Generate representative results:

```bash
python examples/generate_vector_zeeman_results.py
python examples/generate_six_beam_apparatus_results.py
python examples/generate_magnetic_apparatus_results.py
python examples/generate_sequence_results.py
python examples/generate_capture_loading_results.py
python examples/generate_collective_mot_results.py
```

Independent validation requires the pinned optional packages:

```bash
python -m pip install -r requirements-validation.txt
python examples/generate_external_validation_results.py
```

Run the test suite with:

```bash
python -m pytest -q
```

Expensive research grids are intentionally not part of CI.

## Package map

```text
atomic/       isotopes, D lines, hyperfine/Zeeman bases, Wigner transitions
laser/        physical Gaussian beams, polarization, Jones optics, apparatus
magnetic/     quadrupoles, residual fields, Helmholtz/anti-Helmholtz coils
physics/      effective force, rate equations, OBEs, PGC, collective mean field
solvers/      deterministic and photon-event trajectories
simulation/   capture criteria, ensembles, experimental sequences
vacuum.py     Rb vapour pressure/density and configurable loading/loss
foundations.py ballistic, thermal, Gaussian-beam and dipole-trap benchmarks
```
