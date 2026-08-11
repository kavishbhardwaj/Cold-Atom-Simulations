# Cold Atom MOT Simulations

A physically documented, reproducible 3D simulation framework for a six-beam
rubidium-87 vapour-cell magneto-optical trap (MOT).  The current **Phase 1**
implements ideal/coil magnetic fields, an effective two-level semiclassical
radiation-pressure model, adaptive deterministic trajectories, and discrete
photon-event Monte Carlo trajectories with absorption and spontaneous recoil.

This is scientific software under active development—not experiment-control
software and not a quantitatively calibrated prediction of a particular MOT.
The implemented approximation level and missing physics are stated explicitly.

## Phase-1 capabilities

- sourced 87Rb D2 mass, wavelength, lifetime, linewidth, saturation convention,
  recoil scales, hyperfine intervals, Landé factors and normalized hyperfine
  strengths;
- six independent 3D Gaussian travelling beams with power, waist, origin,
  direction, detuning, frequency offset, propagation-relative helicity and
  polarization purity;
- local σ−/π/σ+ polarization decomposition about an arbitrary quantization axis;
- ideal rotated quadrupole fields and composable uniform/gradient/AC stray fields;
- geometrical anti-Helmholtz coils evaluated by converged segmented-wire
  Biot–Savart integration, including tilt, displacement and current imbalance;
- Level-A shared-saturation force with Gaussian intensity, Doppler shift, signed
  Zeeman shift, per-beam momentum transfer and an arbitrary gravity vector;
- adaptive RK45 mean-force trajectories;
- seeded photon-event trajectories with responsible-beam selection, absorption
  momentum and isotropic spontaneous-emission recoil;
- validated YAML configuration, CLI execution and NPZ result data with metadata.

## Fidelity and limitations

| Level | Physical model | Repository status |
|---|---|---|
| A | Effective two-level Doppler/scattering force | **Implemented** |
| B | Multilevel hyperfine/Zeeman rate equations + repump | Phase 2 roadmap |
| C | Multilevel optical Bloch equations | Phase 3 roadmap |
| D | Phase-resolved polarization-gradient/sub-Doppler model | Phase 4 roadmap |
| E | Vapour loading, collision loss and experiment calibration | Phase 5 roadmap |

Level A does **not** model optical pumping, dark states, coherences, repump
dynamics, polarization-gradient cooling, reabsorption, multiple scattering,
loading, collision loss or experimental imperfections not explicitly configured.
It must not be used to claim a sub-Doppler temperature or calibrated atom number.
See [theory](docs/theory.md), [numerics](docs/numerical_methods.md),
[validation](docs/validation.md), and [references](docs/references.md).

## Reproducible quick run

Install the package and execute the committed SI-unit configuration:

```bash
python -m pip install -r requirements-dev.txt
python -m cold_atom_mot simulate configs/rb87_standard_mot.yaml
```

The command writes `results/phase1/phase1_run.npz`, including deterministic and
Monte Carlo arrays plus JSON metadata recording configuration, units, model
fidelity, solver, seed, atom number, time step and package version.

Regenerate the documented reference dataset and every PNG/SVG pair with:

```bash
python examples/generate_phase1_results.py
```

## Selected Phase-1 results

Every displayed PNG has a vector SVG alternative and underlying numerical arrays
in [`phase1_reference.npz`](results/phase1/phase1_reference.npz).

| Apparatus geometry | Physical anti-Helmholtz field |
|---|---|
| ![Six beam and coil geometry](results/phase1/apparatus_geometry.png) | ![Anti-Helmholtz magnetic field](results/phase1/antihelmholtz_field.png) |
| [SVG](results/phase1/apparatus_geometry.svg) | [SVG](results/phase1/antihelmholtz_field.svg) |

| Deterministic force map | Adaptive 3D trajectories |
|---|---|
| ![Force as a function of position and velocity](results/phase1/force_map_x_vx.png) | ![Three dimensional deterministic trajectories](results/phase1/deterministic_trajectories.png) |
| [SVG](results/phase1/force_map_x_vx.svg) | [SVG](results/phase1/deterministic_trajectories.svg) |

| Monte Carlo convergence | Photon-recoil distribution |
|---|---|
| ![Monte Carlo atom number convergence](results/phase1/monte_carlo_convergence.png) | ![Photon recoil velocity distribution](results/phase1/recoil_velocity_distribution.png) |
| [SVG](results/phase1/monte_carlo_convergence.svg) | [SVG](results/phase1/recoil_velocity_distribution.svg) |

![Beam and magnetic parameter sensitivities](results/phase1/parameter_sensitivities.png)

[Vector SVG](results/phase1/parameter_sensitivities.svg)

These are model diagnostics, not fitted experimental data.  The force map uses
the ideal Level-A quadrupole configuration; the field map separately validates
the geometrical coil backend.

## Package map

```text
src/cold_atom_mot/
  atomic/rb87.py             sourced constants and derived recoil scales
  laser/beam.py              independent Gaussian beams and six-beam geometry
  laser/polarization.py      propagation/local spherical polarization bases
  magnetic/fields.py         ideal quadrupole and residual fields
  magnetic/coils.py          segmented-wire circular coils
  physics/force.py           documented Level-A radiation-pressure force
  solvers/deterministic.py   adaptive mean-force trajectories
  solvers/monte_carlo.py     discrete photon events and recoil
  io/config.py               validated YAML configuration
  cli.py                     reproducible command-line runs
```

No empty OBE, rate-equation or sub-Doppler modules are included merely to imply
capabilities that do not yet exist.

## Validation and tests

```bash
python -m pytest -q
python -m compileall src tests
git diff --check
```

The tests cover power normalization, polarization conventions, opposite photon
momenta, Maxwell-consistent gradients, anti-Helmholtz symmetry and segmentation
convergence, coil-error field-zero motion, restoring/damping force symmetry,
3D integration, recoil statistics, fixed-seed reproducibility and timestep
refinement.  CI uses `python -m pytest -q`; expensive research scans are excluded.

## Foundational models

The original `cold_atom.py` remains as an educational collection of ballistic,
thermal-velocity and Gaussian optical-dipole-trap formulas.  It is intentionally
kept separate from the flagship MOT package because those calculations describe
different physical systems and approximation levels.

## Roadmap

1. **Phase 2:** explicit F=1,2 and F'=0,1,2,3 Zeeman bases, cooling/repump optical
   pumping, branching, and matched PyLCP rate-equation comparisons.
2. **Phase 3:** documented Hamiltonian/Lindblad OBE backend with tolerance studies
   and selected PyLCP cross-validation.
3. **Phase 4:** phase-resolved six-beam polarization gradients, light shifts,
   Sisyphus dynamics, diffusion and residual-field suppression studies.
4. **Phase 5:** sourced Rb vapour pressure, independently specified background
   pressure, trajectory-derived capture, loading/loss and experiment calibration.
