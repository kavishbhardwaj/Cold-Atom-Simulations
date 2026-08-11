# Numerical methods and reproducibility

## Units and configuration

Calculations use SI units internally.  `configs/rb87_standard_mot.yaml` records
the Level-A quick-run parameters and fixed seed.  YAML validation rejects
negative power, non-positive waist, missing sections and non-positive time
scales.  `configs/rb87_misaligned_coils.yaml` records a separate physical coil
geometry study.  A production coil configuration will be integrated into the
CLI after Phase-1 convergence and performance profiling.

## Deterministic integration

`scipy.integrate.solve_ivp` with adaptive RK45 is used with configurable
`rtol`, `atol` and maximum step.  The maximum step resolves the mechanical
motion; optical excitation is adiabatically eliminated in Level A, so RK45 does
not resolve the 26 ns lifetime.

## Biot–Savart quadrature

Circular conductors are split into equal straight segments.  Each contribution
is evaluated at its midpoint.  The default is 256 segments per loop.  The test
suite requires the 128-to-256 segment change at a representative off-wire point
to be below 0.1%.  Field Jacobians use centred finite differences.  The wire is
idealized as having zero thickness, so evaluation on it is rejected.

## Photon-event integration

The quick configuration uses `dt=5 ns`, selected from the calculated scattering
probability rather than image smoothness.  The solver rejects a step whenever
`max(R_total dt)>0.1`.  Absorption-beam selection is categorical using the six
rates; spontaneous directions are isotropic.  A `numpy.random.Generator` with an
explicit seed makes every run reproducible.

## Stored results

`results/phase1/phase1_reference.npz` stores arrays used by the figures plus a
JSON metadata record containing configuration, units, solver, model fidelity,
seed, step, atom count and package version.  PNG and SVG files are generated
from the same Matplotlib figure object.  PNG is the GitHub display format; SVG
is the document/vector alternative.

## Quick and research modes

Only `quick` is operational in Phase 1 and is suitable for CI-scale regression.
A future `research` profile will increase ensemble size, coil segmentation and
scan resolution only after convergence is demonstrated.  Naming an expensive
profile without implementing its error controls would be misleading, so no
inactive research config is shipped yet.

## Level-B stationary rate equations

The 24×24 population generator uses a column convention, `dp/dt=A p`. Every
optical and decay process is inserted as equal source/sink terms, and tests
require each generator column to sum to zero within floating-point roundoff. A
normalization row replaces one redundant stationary equation and the resulting
linear system is solved with `numpy.linalg.solve`. Solutions are checked for
normalization and negative population.

The reference position/velocity scans solve this system independently at every
point; no previous solution is used as an initial guess. This is slower than an
adiabatic table interpolation but makes the committed data straightforward to
reproduce. Phase 2 remains a stationary-internal-state approximation: it assumes
internal populations relax faster than external motion. A future time-dependent
population integrator must test that separation where it is not valid.

## Reduced OBE numerics

The stationary two-state OBE is solved as a 4×4 complex Liouvillian null problem
with one row replaced by `Tr(ρ)=1`. Time-dependent density matrices are packed
into eight real variables and integrated with adaptive RK45. Configurable
`rtol`, `atol` and maximum step are expressed relative to the lifetime in the
Phase-3 YAML file. Tests compare two tolerance/step settings and require their
final density matrices to converge.

The committed power–detuning map varies both parameters instead of presenting a
single unexplained power trace. The waist figure presents fixed-power and
fixed-peak-intensity conditions at the same x=2 mm evaluation point. All arrays
and the held-fixed choices are stored in `results/phase3/phase3_reference.npz`.
