# Validation status

## Completed effective-model checks

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
- effective force is restoring and damping near the origin and respects inversion
  symmetry when gravity and stray fields vanish;
- isotropic recoil has zero mean and component variance 1/3;
- a fixed seed exactly reproduces Monte Carlo trajectories;
- Monte Carlo mean velocity is stable under a factor-two timestep refinement.

## Meaning of validation

These are internal analytic, symmetry and convergence checks.  They do **not**
constitute experimental validation.  No measured capture fraction, temperature,
cloud size, loading rate or atom number is predicted.

## PyLCP status

PyLCP is cited and planned as an independent rate-equation/coherent reference.  It is
not a effective-model dependency and no PyLCP comparison is claimed.  A scientifically
useful comparison must match basis, beam convention, magnetic Hamiltonian,
linewidth and saturation convention; a superficial plot overlay would not be
validation.

## Next validation gates

the multilevel model must compare Clebsch–Gordan-resolved populations and force profiles with
PyLCP rate equations.  the OBE benchmark must compare steady-state OBEs and tolerance
refinement.  the polarization-gradient model must reproduce known 1D polarization-gradient limits before
attempting phase-coherent 3D claims.  the loading/loss model requires published or measured
collision inputs and experiment-specific calibration.

## multilevel checks

- the Wigner/CG generation reproduces known Clebsch–Gordan values and selection
  rule zeros;
- the basis contains exactly 8 ground and 16 excited Zeeman states;
- every generated transition satisfies `m'=m+q` and the stretched cycling
  transition has unit strength under the saturation convention;
- spontaneous branching from every excited Zeeman state sums to one;
- every rate-generator column sums to zero;
- stationary populations are non-negative and normalized;
- adding repump light transfers steady-state population back into ground F=2;
- the ideal rate-equation force is restoring in position and damping in velocity.

These checks validate implementation consistency and analytic limits. They are
not an independent atomic-physics benchmark. A matched PyLCP comparison remains
required before claiming external rate-equation validation; the comparison must use
the same hyperfine basis, saturation convention, beam polarizations, detunings
and magnetic Hamiltonian.

## two-level OBE reduced-OBE checks

- the stationary density matrix matches the exact two-level excited-population
  expression over four decades of saturation and red/blue detunings;
- density matrices are normalized, Hermitian and positive semidefinite;
- zero drive relaxes to the ground state and produces zero radiation force;
- time evolution preserves trace and approaches the stationary state;
- tighter integration tolerance and smaller maximum step converge;
- a pure-dephasing Lindblad channel agrees with its generalized analytical
  stationary population;
- the atomic registry distinguishes 85Rb/87Rb and D1/D2 while exposing that only
  87Rb D2 has high-fidelity solver support.

This validates the reduced solver, not a full multilevel MOT OBE. External PyLCP
comparison and a 24-state coherent Hamiltonian remain validation gates before
claiming full coherent MOT support.

## polarization-gradient checks and boundary

Tests require a periodic transverse field; normalized local σ−/π/σ+ fractions;
a probability-conserving pumping generator with non-negative transfer rates;
physical stationary populations; finite state-resolved shifts and forces; the
expected axial Zeeman spacing; and force convergence under grid refinement.
These are internal checks, not PyLCP or experimental validation. A multilevel
OBE retaining ground coherences is the next external-validation gate.

## Vapour, capture, and loading checks

- the surface-flux speed sampler reproduces its analytical mean;
- all sampled spherical trajectories point inward and fixed seeds reproduce;
- stratified sample weights equal their analytical flux probabilities;
- loading rate equals incident isotope flux times weighted capture probability;
- rubidium partial pressure and non-Rb background pressure remain independent;
- the one-sided flux equals `n<v>/4`;
- zero-loss and one-body loading equations reproduce their analytical limits;
- two-body loss is rejected without a positive effective volume.
- zero-capture strata retain a nonzero Wilson upper bound;
- capture classification is stable under RK45 step/tolerance refinement;
- arbitrary initial populations reproduce analytical one-body solutions;
- nonlinear steady state satisfies the positive-root balance equation;
- Gaussian-cloud effective volume matches the analytical density integral;
- background component losses add and aggregate/component double counting is rejected.

These tests validate sampling and bookkeeping, not collision cross sections or
an experimental loading rate. The acceptance surface and effective loss inputs
must be matched to an apparatus before quantitative comparison.

## Vector hyperfine–Zeeman checks

- exact zero-field Hamiltonians reproduce tabulated hyperfine energies;
- weak-field slopes recover generated Landé `g_F` factors;
- equal-magnitude x, y, z, and tilted fields have identical spectra;
- the vector Hamiltonian is Hermitian and continuous through `B=0`;
- larger fields resolve a nonlinear departure from `g_F m_F mu_B B`;
- the transverse-field multilevel Liouvillian remains trace preserving.

## Residual fields and polarization-gradient cooling

The former axial residual-field scan is not a vector-field validation.  Its
five-population adiabatic model projects B onto a fixed quantization axis and
now rejects transverse B.  It reports a matched recoil-event diffusion tensor
(directional absorption plus isotropic emission), but omits internal-state and
dipole-force fluctuations; consequently no Einstein temperature is claimed.

The companion full-vector 24-state OBE scan covers 0, 0.1, 0.3, 1, 3, 10, 30,
100, 300 and 1000 mG independently along x, y and z.  These equal 0--100 µT or
0--1e-4 T; an Earth-scale 500 mG reference is marked.  For F=2, the Larmor scale
is 699.6 Hz/mG, while the weak-drive optical-pumping estimate for the explicitly
configured -3 Gamma, s=0.08 recipe is 6.56 kHz.  Thus the scales cross near
9.4 mG, but this is a competition-of-rates marker, **not** a measured 10%
friction threshold.

The static OBE scan retains vector Zeeman mixing, optical coherences, full
excited hyperfine structure and branching.  It finds orientation-dependent
density-matrix coherences, but symmetry makes Fx(v=0) vanish.  A defensible
beta(B) requires long-time moving-lattice calculations at both velocity signs;
a defensible temperature additionally requires internal/dipole-force noise from
the same OBE.  Neither is inferred from the present point scan.  Dalibard and
Cohen-Tannoudji, JOSA B **6**, 2023 (1989), DOI 10.1364/JOSAB.6.002023, provides
the primary theoretical context that magnetic precession disrupts Sisyphus
optical pumping. Geometry and atomic assumptions are not matched, so the
comparison is qualitative and parameters were not fitted.
