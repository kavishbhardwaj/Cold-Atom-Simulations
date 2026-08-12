# Cold-atom and magneto-optical-trap simulations

A reproducible SI-unit Python framework for rubidium spectroscopy, realistic
six-beam MOT apparatus, semiclassical and multilevel light forces, optical Bloch
benchmarks, polarization-gradient cooling, stochastic recoil, capture, and
configurable vapour loading/loss. **87Rb D2 is the default reference MOT**, but
the atomic engine generates the full hyperfine/Zeeman structure of 85Rb and
87Rb on both D lines.

This is inspectable scientific software, not a calibrated prediction of a
particular experiment. Every solver documents what it neglects.

## Atomic and solver support

| System | Atomic basis/data | Effective MOT | Rate-equation MOT | Multilevel OBE | Polarization gradient / gray molasses |
|---|---:|---:|---:|---:|---:|
| **87Rb D2** | full | yes | **F=2→F′=3**, repump **F=1→F′=2** | sparse framework; research-scale validation pending | adiabatic F=2→F′=3 population model |
| **85Rb D2** | full | yes | **F=3→F′=4**, repump **F=2→F′=3** | operator framework | not validated |
| 87Rb D1 | full | line benchmark only | no conventional closed-cycle MOT | operator framework | gray molasses not implemented |
| 85Rb D1 | full | line benchmark only | no conventional closed-cycle MOT | operator framework | gray molasses not implemented |

D2 provides the convenient stretched cycling transition used by conventional
MOTs. D1 is not inferior: its Λ systems, Raman coherence and dark states are
central to gray molasses, but those require a coherently validated solver and
are not faked here. See [atomic systems](docs/atomic_systems.md).

## Physical model hierarchy

| Capability | Purpose and approximation |
|---|---|
| Foundations | Ballistics, thermal distributions, Gaussian propagation and dipole traps |
| Effective semiclassical MOT | Fast shared-saturation Doppler/Zeeman force and trajectories |
| Multilevel rate equations | Generated hyperfine/Zeeman populations, cooling and repump; no coherences |
| Two-level OBE benchmark | Exact analytical/numerical coherence benchmark |
| Sparse multilevel OBE | Exact vector hyperfine–Zeeman Hamiltonian and collapse operators; full six-beam validation remains |
| Polarization-gradient model | Phase/coherence-group-resolved light shifts, pumping and Sisyphus force; populations only |
| Photon-event Monte Carlo | Exact ≥1 Poisson-event probability per step, absorption and spontaneous recoil |
| Capture and loading/loss | Thermal surface-flux trajectories, isotope-resolved loading, and calibrated loss inputs |

Detailed equations and boundaries: [model hierarchy](docs/model_hierarchy.md),
[cooling physics](docs/cooling_physics.md), [numerics](docs/numerical_methods.md),
and [validation](docs/validation.md).
Time-dependent laboratory controls are described in [experimental sequences](docs/experimental_sequences.md).
Independent physical beams and Jones optics are described in [six-beam apparatus](docs/six_beam_apparatus.md).
Bias compensation, coil imperfections, and switching fields are described in [magnetic apparatus](docs/magnetic_apparatus.md).

## Validation status

| Benchmark | Status |
|---|---|
| Two-level OBE vs analytic formulas and QuTiP | **ANALYTICALLY / INDEPENDENT-SOFTWARE VERIFIED** |
| Normalized 1-D two-beam force vs PyLCP | **INDEPENDENT-SOFTWARE VERIFIED** |
| 87Rb ground vector-Zeeman spectrum vs PyLCP | **INDEPENDENT-SOFTWARE VERIFIED** |
| Full 24-state OBE forces/populations vs external software | **NOT YET VALIDATED** |
| PGC mechanism vs primary theory | **LITERATURE-TREND VERIFIED** |
| Quantitative 87Rb MOT/PGC experiment | **NOT YET VALIDATED** |

Exact conventions, numerical residuals, unmatched assumptions, and citations
are in [validation](docs/validation.md).

## Selected results

| Multilevel force and populations | Polarization-gradient force |
|---|---|
| ![Effective and multilevel force](results/multilevel/effective_vs_multilevel_force.png) | ![Sub-Doppler and Doppler forces](results/polarization_gradient/subdoppler_force_velocity.png) |

| Damping versus power | Beam-waist capture tradeoff |
|---|---|
| ![Damping and scattering](results/parameter_studies/damping_power_physics.png) | ![Waist and capture](results/parameter_studies/beam_waist_capture.png) |

| Physical anti-Helmholtz field | Deterministic trajectories |
|---|---|
| ![Coil field](results/effective_mot/antihelmholtz_field.png) | ![Trajectories](results/effective_mot/deterministic_trajectories.png) |

| Vapour capture and loading | Loading-loss sensitivity |
|---|---|
| ![Vapour pressure, capture, and loading](results/capture_loading/vapor_capture_loading.png) | ![Loading curves under calibrated loss rates](results/capture_loading/loading_loss_sensitivity.png) |

All captions, held-fixed parameters, SVG alternatives, NPZ data and limitations
are in the **[results gallery](results/README.md)**.

## Reproduce

```bash
python -m pip install -r requirements-dev.txt
python -m cold_atom_mot simulate configs/rb87_d2_mot.yaml
python -m cold_atom_mot rate-equation configs/rb87_d2_multilevel.yaml
python -m cold_atom_mot obe configs/rb87_d2_two_level_obe.yaml
python -m cold_atom_mot subdoppler configs/rb87_d2_polarization_gradient.yaml
python -m cold_atom_mot loading configs/rb_vapor_loading.yaml

python examples/generate_foundations.py
python examples/generate_effective_mot_results.py
python examples/generate_multilevel_results.py
python examples/generate_optical_bloch_results.py
python examples/generate_polarization_gradient_results.py
python examples/generate_parameter_studies.py
python examples/generate_capture_loading_results.py
```

## Package map

```text
atomic/       isotope, line, hyperfine basis and Wigner-generated transitions
laser/        Gaussian beams, polarization and coherence groups
magnetic/     quadrupoles, residual fields and segmented Biot–Savart coils
physics/      effective force, rate equations, OBEs and polarization gradients
solvers/      deterministic and photon-event trajectories
simulation/   explicit capture criteria and ensemble metrics
vacuum.py     Rb vapour pressure/density and configurable loading/loss
foundations.py ballistic, thermal, Gaussian and dipole-trap benchmarks
```

CI deliberately uses `python -m pytest -q`; expensive research grids are not run
in CI.
