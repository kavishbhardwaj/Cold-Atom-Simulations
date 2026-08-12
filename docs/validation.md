# Validation status

This document describes the **current** scientific status of the repository. Historical implementation milestones are intentionally omitted so that “implemented”, “tested”, and “externally validated” are not confused.

## Status language

- **INTERNAL TESTED** — analytical limits, symmetries, conservation laws, numerical convergence, or regression tests pass inside this repository.
- **ANALYTICALLY VERIFIED** — the implementation reproduces a closed-form result under matched assumptions.
- **INDEPENDENT-SOFTWARE VERIFIED** — a matched calculation agrees with a separate package through its public API.
- **LITERATURE-TREND VERIFIED** — the qualitative mechanism/trend agrees with primary literature, but assumptions are not matched well enough for a numerical residual claim.
- **EXPERIMENTALLY COMPARED** — a published experiment is discussed, but the present apparatus is not sufficiently specified to claim quantitative reproduction.
- **NOT YET VALIDATED** — the feature may be implemented and internally tested, but the stated high-fidelity claim has not passed the relevant external gate.

A successful test of one layer does **not** validate a more complicated layer built on top of it.

## Independent validation matrix

| Claim | Analytical / internal | QuTiP 5.3.1 | PyLCP 1.0.2 | Literature / experiment | Current status |
|---|---|---|---|---|---|
| Two-level OBE steady state, Liouvillian, Rabi dynamics, decay | exact formulas + convergence tests | matched public-API comparison | not needed | standard OBE physics | **ANALYTICALLY + INDEPENDENT-SOFTWARE VERIFIED** |
| Normalized 1-D two-beam Doppler force | exact low-level force expression | not needed | matched `heuristiceq` comparison | standard Doppler theory | **INDEPENDENT-SOFTWARE VERIFIED** |
| 87Rb ground vector-Zeeman spectrum | weak-field slopes, rotation covariance, exact hyperfine construction | not compared | matched `hyperfine_coupled` comparison | Breit-Rabi/Zeeman physics | **INDEPENDENT-SOFTWARE VERIFIED** |
| 24-state moving 87Rb D2 OBE | trace, Hermiticity, positivity, Doppler/frame/force/convergence tests | two-level limit only | full multilevel case not yet compared | not quantitatively matched | **IMPLEMENTED + INTERNAL TESTED; EXTERNAL VALIDATION PENDING** |
| Adiabatic population PGC | pumping, light-shift, periodicity and force checks | not comparable | not yet matched | Dalibard–Cohen-Tannoudji mechanism | **LITERATURE-TREND VERIFIED** |
| Quantitative PGC temperature / residual-field tolerance | incomplete matched diffusion | not compared | not compared | qualitative magnetic-disruption context | **NOT YET VALIDATED** |
| MOT capture/loading | analytical flux/loading limits + convergence tests | not applicable | simple force only | apparatus parameters incomplete | **INTERNAL TESTED; NOT EXPERIMENTALLY CALIBRATED** |
| Gaussian collective MOT / multiple scattering | Gaussian integrals + low-density recovery | not applicable | not compared | Walker/Sesko/Wieman trend | **LITERATURE-TREND VERIFIED** |

## Quantitative independent-software results

### QuTiP: two-level optical Bloch equations

The benchmark constructs exactly the internal convention

`H/hbar = [[0, Omega*/2], [Omega/2, -Delta]]`

with collapse operator

`sqrt(Gamma) |g><e|`

and saturation convention

`s = 2 |Omega|^2 / Gamma^2`.

The comparison spans multiple detunings and saturations, includes damped Rabi dynamics and pure spontaneous decay, and explicitly maps the different Liouvillian vectorization conventions.

Representative maximum residuals from the committed benchmark:

- steady-state excited-population absolute difference: **5.55e-17**;
- maximum relative excited-population difference: **4.84e-16**;
- Liouvillian entries: within **2e-14** test tolerance;
- Rabi/decay density matrices: within **2e-8**;
- pure spontaneous decay: within **3e-8**.

No fitted parameter is used.

### PyLCP: normalized two-beam force

The public PyLCP `heuristiceq` comparison uses matched dimensionless units

`Gamma = k = hbar = 1`, `s = 0.05` per beam, `Delta = -2 Gamma`,

with opposite propagation directions and the same shared saturation denominator.
Across 101 velocities from `-0.5 Gamma/k` to `+0.5 Gamma/k`:

- force RMS difference: **2.15e-19 hbar k Gamma**;
- maximum relative force difference: **7.93e-15**.

This verifies the normalized two-level/two-beam force. It does **not** validate the full 87Rb multilevel MOT or PGC solver.

### PyLCP: 87Rb vector-Zeeman spectrum

The independent comparison uses PyLCP's public `hyperfine_coupled` helper with the same 87Rb nuclear spin, electronic angular momentum, hyperfine constant, electronic `g` factor, and Bohr magneton. The nuclear-`g` sign-convention difference is mapped explicitly.

At magnetic fields of `0`, `1 uT`, and `100 uT`, the maximum absolute ground-manifold spectral difference is approximately **0.57 Hz**.

This independently verifies the implemented 87Rb ground hyperfine-Zeeman spectrum over the tested field range. It does not by itself validate light coupling, optical pumping, or force.

## Implemented high-fidelity 87Rb D2 OBE

The repository now contains a **24-state moving-atom OBE** for the full 87Rb D2 hyperfine-Zeeman basis: 8 ground states and 16 excited states. It includes full vector Zeeman structure within the adopted rotating-frame approximations, cooling/repump frequencies, beam-resolved Doppler shifts, off-resonant excited-hyperfine coupling, spontaneous branching, phase/coherence-group handling, and force operators from the interaction-Hamiltonian gradient.

Internal tests cover:

- complete 24-state basis and transition graph;
- trace, Hermiticity and positivity tolerances;
- beam-specific Doppler sign and velocity-dependent force;
- stationary versus time-dependent Hamiltonian detection;
- cooling + repump block-rotating frame;
- controlled cross-ground-manifold RWA diagnostics;
- incoherent phase-ensemble refinement;
- analytic interaction-gradient validation;
- travelling-wave force limit;
- time-average convergence and deliberate rejection of unconverged RESEARCH cases.

This is **implemented and internally tested**, but the complete multilevel force and populations are **not yet independently verified against a matched PyLCP 87Rb D2 calculation**. That remains the principal external validation gate.

## Residual magnetic fields and sub-Doppler cooling

The reduced five-population PGC model supports only a magnetic field parallel to its fixed quantization axis. It now rejects unsupported transverse fields instead of silently projecting them.

A companion full-vector 24-state OBE diagnostic scans residual fields independently along x, y, and z from `0` to `1000 mG`. For the configured `87Rb D2`, `Delta=-3 Gamma`, `s=0.08/beam` recipe:

- F=2 Larmor scale: approximately **699.6 Hz/mG**;
- weak-drive optical-pumping scale: approximately **6.56 kHz**;
- the two simple rates cross near **9.4 mG**.

The 9.4 mG value is a **competition-of-timescales marker**, not a temperature threshold and not a “10% PGC degradation” criterion.

The vector OBE shows orientation-dependent coherence changes, but a defensible `beta(B)` requires converged moving-lattice force calculations at positive and negative velocities. A defensible `T(B)` additionally requires momentum diffusion/force noise at compatible fidelity, including internal-state and dipole-force fluctuations. The repository therefore does **not** manufacture a quantitative sub-Doppler temperature from an inconsistent friction/diffusion pair.

Dalibard and Cohen-Tannoudji, *JOSA B* **6**, 2023 (1989), DOI `10.1364/JOSAB.6.002023`, provides the primary trend-level context for polarization-gradient cooling and magnetic-precession disruption. Geometry, dimensionality and diffusion assumptions are not sufficiently matched for a numerical reproduction claim.

## Effective, rate-equation, and apparatus checks

The lower-cost models retain extensive internal validation:

- generated transition strengths and selection rules;
- normalized spontaneous branching;
- rate-generator probability conservation;
- repump population transfer;
- restoring/damping symmetries;
- Gaussian-beam power normalization;
- ideal quadrupole `trace(gradient)=0`;
- segmented anti-Helmholtz convergence and field-zero behavior;
- deterministic and photon-event trajectory reproducibility;
- recoil statistics and timestep refinement;
- six-beam polarization/Jones normalization;
- three-axis compensation and coil-response checks;
- experimental-sequence ramp and switching regressions.

These are implementation and numerical checks, not experimental calibration.

## Vapour capture, loading, and loss

The loading framework is internally checked against:

- the one-sided thermal flux `n<v>/4`;
- analytical surface-flux speed statistics and cosine incidence;
- isotope-resolved flux weighting;
- stratified rare-event sampling and Wilson intervals;
- capture-classification solver refinement;
- arbitrary-`N0` one-body loading solutions;
- nonlinear steady-state balance;
- Gaussian two-body effective volume;
- separated background/hot-Rb loss bookkeeping and double-count protection.

The acceptance sphere is an **acceptance boundary**, not a complete chamber model. Collision-loss coefficients, wall recycling, dispenser geometry and apparatus-specific pressure composition remain measured/literature inputs. Atom-number curves using assumed inputs are scenarios, not calibrated predictions.

## Collective MOT extension

The optional Gaussian collective model adds density-dependent two-body loss, Beer-Lambert shadowing, a Walker/Sesko/Wieman-style Coulomb multiple-scattering approximation, optical-depth diagnostics and a radiation-trapping recoil proxy. It exactly recovers the independent-atom loading limit when collective terms are removed.

Its cross sections, temperature, restoring coefficient and two-body coefficient must retain explicit provenance. The model is **LITERATURE-TREND VERIFIED**, not exact radiative transfer and not experimentally validated.

Primary trend reference: Walker, Sesko, and Wieman, *Phys. Rev. Lett.* **64**, 408 (1990), DOI `10.1103/PhysRevLett.64.408`.

## Published experimental context

Lett et al., *JOSA B* **6**, 2084 (1989), DOI `10.1364/JOSAB.6.002084`, is retained as foundational sub-Doppler experimental context but uses sodium rather than 87Rb.

Townsend et al., *Phys. Rev. A* **52**, 1423 (1995), DOI `10.1103/PhysRevA.52.1423`, provides established MOT context. The present repository does not claim quantitative reproduction because beam profiles, coherence topology, residual vector field, coil transfer function, chamber/capture geometry and loss parameters are not all matched.

No experimental parameter is silently tuned to reduce a discrepancy.

## Next validation gates

The most valuable next steps are deliberately narrow:

1. matched PyLCP **87Rb D2 multilevel** populations and force versus velocity/position;
2. matched 1-D `F=2 -> F'=3` polarization-gradient benchmark using the coherent OBE;
3. OBE-consistent force-noise/diffusion before reporting quantitative sub-Doppler `T(B)`;
4. one fully specified 87Rb MOT/PGC experiment with measured or published apparatus inputs;
5. capture/loading comparison against a fully specified vapour-cell apparatus.

The reproducible external-validation data and software versions are stored under `results/validation/`.
