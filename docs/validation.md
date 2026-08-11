# Validation status

## Completed Phase-1 checks

- atomic transition-strength tables normalize within their documented manifolds;
- Gaussian intensity integrates to configured power;
- propagation-relative circular polarization reverses correctly and local
  spherical fractions normalize;
- opposing beam momenta cancel for an ideal six-beam geometry;
- the ideal quadrupole has `trace(gradient)=0`;
- a symmetric segmented anti-Helmholtz pair has a central zero and the expected
  near-centre gradient ratio `(1,1,-2)`;
- 128-to-256 coil segmentation convergence is tested away from the wire;
- tilt/offset moves the numerical field zero;
- Level-A force is restoring and damping near the origin and respects inversion
  symmetry when gravity and stray fields vanish;
- isotropic recoil has zero mean and component variance 1/3;
- a fixed seed exactly reproduces Monte Carlo trajectories;
- Monte Carlo mean velocity is stable under a factor-two timestep refinement.

## Meaning of validation

These are internal analytic, symmetry and convergence checks.  They do **not**
constitute experimental validation.  No measured capture fraction, temperature,
cloud size, loading rate or atom number is predicted.

## PyLCP status

PyLCP is cited and planned as an independent Level-B/Level-C reference.  It is
not a Phase-1 dependency and no PyLCP comparison is claimed.  A scientifically
useful comparison must match basis, beam convention, magnetic Hamiltonian,
linewidth and saturation convention; a superficial plot overlay would not be
validation.

## Next validation gates

Phase 2 must compare Clebsch–Gordan-resolved populations and force profiles with
PyLCP rate equations.  Phase 3 must compare steady-state OBEs and tolerance
refinement.  Phase 4 must reproduce known 1D polarization-gradient limits before
attempting phase-coherent 3D claims.  Phase 5 requires published or measured
collision inputs and experiment-specific calibration.
