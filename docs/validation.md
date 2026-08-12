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

## Independent validation matrix

Labels are deliberately narrow: verification of one limit does not validate a
different solver or apparatus model.

| Claim | Analytical | Independent internal | QuTiP 5.3.1 | PyLCP 1.0.2 | Published theory | Published experiment |
|---|---|---|---|---|---|---|
| Two-level OBE steady state, Liouvillian, Rabi dynamics, decay | **ANALYTICALLY VERIFIED** | **INTERNAL TESTED** | **INDEPENDENT-SOFTWARE VERIFIED** | not needed | standard OBE | **NOT YET VALIDATED** |
| Normalized 1-D two-beam molasses force | **ANALYTICALLY VERIFIED** | **INTERNAL TESTED** | not needed | **INDEPENDENT-SOFTWARE VERIFIED** | standard Doppler theory | **NOT YET VALIDATED** |
| 87Rb ground vector-Zeeman spectrum | weak-field and rotation checks | independent uncoupled/coupled construction | not compared | **INDEPENDENT-SOFTWARE VERIFIED** | Breit-Rabi trend | spectroscopy not compared |
| 24-state moving OBE | two-level/trace limits only | **INTERNAL TESTED** | full case **NOT YET VALIDATED** | full case **NOT YET VALIDATED** | **NOT YET VALIDATED** | **NOT YET VALIDATED** |
| Adiabatic population PGC | limiting formulas | **INTERNAL TESTED** | not comparable | not compared | **LITERATURE-TREND VERIFIED** | **NOT YET VALIDATED** |
| MOT/loading predictions | component formulas | **INTERNAL TESTED** | not compared | simple normalized force only | qualitative | **EXPERIMENTALLY COMPARED; not quantitatively reproduced** |
| Gaussian collective MOT / multiple scattering | Gaussian integrals and low-density limit | **INTERNAL TESTED** | not applicable | not compared | **LITERATURE-TREND VERIFIED** | **NOT YET VALIDATED** |

### QuTiP methodology and errors

The benchmark constructs through QuTiP's public API exactly the internal
Hamiltonian `[[0,Omega*/2],[Omega/2,-Delta]]`, collapse operator
`sqrt(Gamma)|g><e|`, and saturation convention
`s=2|Omega|^2/Gamma^2`. Row-stacked internal and column-stacked QuTiP
Liouvillians are related by an explicit permutation. Tests span detuning and
saturation, damped Rabi oscillations, and pure spontaneous decay. The generated
NPZ records maximum population errors; no fitted parameter is used.

### PyLCP methodology and scope

PyLCP is called only through its public API. Its `heuristiceq` two-beam plane-
wave case is matched in dimensionless units (`Gamma=k=hbar=1`), including each
beam's `s=0.05`, detuning `-2 Gamma`, propagation direction, and shared
saturation denominator. Force is compared at 101 velocities by RMS and maximum
relative error. This verifies the normalized heuristic force, **not** the full
87Rb multilevel implementation.

The independent 87Rb comparison uses PyLCP's public `hyperfine_coupled` helper
with the same `I=3/2`, `J=1/2`, hyperfine A, electronic g factor and Bohr
magneton. PyLCP's documented nuclear-g convention has the opposite sign to the
Steck convention used here, so the mapping is explicit. Ground-manifold spectra
at 0, 1 uT and 100 uT agree within 0.7 Hz. A matched PyLCP 87Rb D2 force,
population, 1-D MOT, and 3-D MOT benchmark remains **NOT YET VALIDATED**; the
repository does not disguise the normalized two-level comparison as that work.

### Literature and experiment

Dalibard and Cohen-Tannoudji, *JOSA B* **6**, 2023 (1989), DOI
`10.1364/JOSAB.6.002023`, is the primary theoretical benchmark. The repository
reproduces the red-detuned polarization-gradient mechanism—spatial light shifts,
state-dependent pumping, sub-Doppler odd force and magnetic-precession
suppression—as a **LITERATURE-TREND VERIFIED** result. Geometry, reduced
manifolds, dimensionality and diffusion are not sufficiently matched for a
numerical curve claim.

Lett et al., *JOSA B* **6**, 2084 (1989), DOI `10.1364/JOSAB.6.002084`, is an
experimental observation of sub-Doppler cooling, but it uses sodium rather than
87Rb and therefore is not an 87Rb validation. Townsend et al., *Phys. Rev. A*
**52**, 1423 (1995), DOI `10.1103/PhysRevA.52.1423`, supplies an established MOT
experimental context, but beam profiles, phase topology, residual vector field,
coil transfer function and capture boundary are insufficiently matched to the
present reference apparatus. It is labeled **EXPERIMENTALLY COMPARED; not
quantitatively reproduced**, with measured/assumed/unknown inputs kept separate.
No experimental parameter was tuned to reduce a residual.
