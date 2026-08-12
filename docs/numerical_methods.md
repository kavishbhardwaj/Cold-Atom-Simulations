# Numerical methods and reproducibility

## Units and configuration

Calculations use SI units internally.  `configs/rb87_d2_mot.yaml` records
the effective quick-run parameters and fixed seed.  YAML validation rejects
negative power, non-positive waist, missing sections and non-positive time
scales.  `configs/rb87_misaligned_coils.yaml` records a separate physical coil
geometry study.  A production coil configuration will be integrated into the
CLI after effective-model convergence and performance profiling.

## Deterministic integration

`scipy.integrate.solve_ivp` with adaptive RK45 is used with configurable
`rtol`, `atol` and maximum step.  The maximum step resolves the mechanical
motion; optical excitation is adiabatically eliminated in effective, so RK45 does
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

`results/effective_mot/effective_mot_reference.npz` stores arrays used by the figures plus a
JSON metadata record containing configuration, units, solver, model fidelity,
seed, step, atom count and package version.  PNG and SVG files are generated
from the same Matplotlib figure object.  PNG is the GitHub display format; SVG
is the document/vector alternative.

## Quick and research modes

Only `quick` is operational in the effective model and is suitable for CI-scale regression.
A future `research` profile will increase ensemble size, coil segmentation and
scan resolution only after convergence is demonstrated.  Naming an expensive
profile without implementing its error controls would be misleading, so no
inactive research config is shipped yet.

## rate-equation stationary rate equations

The 24×24 population generator uses a column convention, `dp/dt=A p`. Every
optical and decay process is inserted as equal source/sink terms, and tests
require each generator column to sum to zero within floating-point roundoff. A
normalization row replaces one redundant stationary equation and the resulting
linear system is solved with `numpy.linalg.solve`. Solutions are checked for
normalization and negative population.

The reference position/velocity scans solve this system independently at every
point; no previous solution is used as an initial guess. This is slower than an
adiabatic table interpolation but makes the committed data straightforward to
reproduce. the multilevel model remains a stationary-internal-state approximation: it assumes
internal populations relax faster than external motion. A future time-dependent
population integrator must test that separation where it is not valid.

## Reduced OBE numerics

The stationary two-state OBE is solved as a 4×4 complex Liouvillian null problem
with one row replaced by `Tr(ρ)=1`. Time-dependent density matrices are packed
into eight real variables and integrated with adaptive RK45. Configurable
`rtol`, `atol` and maximum step are expressed relative to the lifetime in the
two-level OBE YAML file. Tests compare two tolerance/step settings and require their
final density matrices to converge.

The committed power–detuning map varies both parameters instead of presenting a
single unexplained power trace. The waist figure presents fixed-power and
fixed-peak-intensity conditions at the same x=2 mm evaluation point. All arrays
and the held-fixed choices are stored in `results/optical_bloch/optical_bloch_reference.npz`.

## Phase-resolved population dynamics

polarization-gradient integrates five ground populations with adaptive RK45 along `x=v t`.
Interference with transverse beams makes the general six-beam field periodic
over one wavelength, so averages use complete λ periods and discard the first
half. The committed convergence scan refines 20–80 samples per period. Numerical
light-shift gradients use `2e-5/k`. Phase, detuning, saturation and uniform field
are stored in YAML and NPZ metadata.

## Vapour loading and loss

Rubidium partial pressure is either supplied directly or obtained from the cited
Alcock fit at the **reservoir/cold-spot temperature**. The vapour kinetic
temperature separately sets `P_Rb/(k_B T_vapour)` and incident velocities; the
background temperature sets background density and relative collision speed.
These can differ when a reservoir/cold spot, cell body, and vacuum environment
are not isothermal. This is distinct from non-Rb background pressure. The loading integrator accepts, but never invents, a
trajectory/geometry-derived loading rate, a calibrated one-body loss rate, and
an optional phenomenological two-body coefficient/effective volume. No committed
curve is presented as a calibrated atom-number prediction.

Atoms incident on the spherical acceptance boundary are sampled from the
one-sided equilibrium flux. The total flux per unit area is
`n sqrt(k_B T/(2 pi m)) = n <v>/4`; its speed distribution is proportional to
`v^3 exp[-m v^2/(2 k_B T)]`, not the bulk Maxwell `v^2` law. Incidence angles
follow the cosine law. Configurable speed strata resolve the extremely rare
slow tail without pretending a small unstratified sample can measure it. Each
stratum carries its exact flux-CDF probability. Geometric tail strata extend to
a configured maximum; a Wilson 95% upper bound from the terminal zero-capture
stratum bounds—not declares zero—the unresolved loading contribution. The result
stores the last speed, omitted flux, upper bound, and whether tolerance was met.

The capture criterion requires an atom to remain inside the specified radius
and below the speed threshold for the configured dwell time. Multiplying its
weighted probability by isotope-specific incident flux produces `R`. The loading
ODE then uses independently configured one-body background/hot-Rb losses and an
optional phenomenological two-body coefficient. A kinetic background-loss mode
is available only when the user supplies gas mass and an effective ejection
cross section.

The temperature scan integrates a temperature-independent numerical response
`P_capture(v,b)` against the surface-flux distribution at each vapour
temperature. It therefore changes density, mean flux, and within-bin thermal
weight consistently rather than reusing conditional samples drawn at one T.

For a Gaussian cloud with rms widths `(sigma_x,sigma_y,sigma_z)`,
`integral(n²dV)=N²/[8 pi^(3/2) sigma_x sigma_y sigma_z]`; the corresponding
two-body effective volume may replace the explicit volume. Widths are inputs,
not inferred from the approximate polarization-gradient temperature.

## Vector magnetic Hamiltonian in coherent models

The sparse multilevel OBE constructs hyperfine and vector-Zeeman blocks in a
stable uncoupled `|m_I,m_J>` basis and transforms them to the coupled basis used
by the dipole graph. It evaluates the complete `B(r,t)` from the configured
quadrupole, residual, composite, or coil field. Transverse magnetic mixing is
retained and the Hamiltonian remains continuous at field zero. The effective
force, rate equations, and population-only polarization-gradient model retain
documented reduced Zeeman approximations for speed; they are not equivalent to
this coherent Hamiltonian.

## Moving-atom multilevel OBE

`MultilevelOBE` uses the complete generated D2 dipole graph, rather than
discarding couplings outside a laser's nominal tuning manifold.  For 87Rb this
is the 24-state basis (8 ground and 16 excited states), a 576-component complex
density matrix, and spontaneous collapse operators made from the calculated
branching matrix.  Cooling and repump target labels determine laser frequency;
they do not filter the allowed excited hyperfine states.

One global excited-manifold optical reference is used.  Every beam then carries
its own phase

`k_i . (r0 + v t) - delta_omega_i t + phi_i`,

where `delta_omega_i` contains the selected hyperfine tuning, detuning, and AOM
offset.  Consequently the trajectory-frame frequency is independently
`delta_omega_i-k_i.v`.  Equal trajectory-frame frequencies admit a stationary
rotating-frame sparse null solve.  Unequal frequencies are not folded into an
effective detuning: the sparse, explicitly time-dependent Liouvillian is
integrated and its late-time force is averaged.

The per-beam force is evaluated as `-hbar Tr(rho gradient(H_i/hbar))`.  A spatial
central derivative includes both the travelling-wave radiation-pressure term
and Gaussian-envelope dipole force, including stimulated/coherent effects.
Spontaneous emission has zero mean recoil in the isotropic approximation and
therefore adds diffusion but no mean force.

`quick` mode uses 12 lifetimes and 49 output samples for a time-dependent point;
`research` uses 40 lifetimes and 241 samples with tighter integration tolerance.
These modes are intended for representative phase-space points and short 1-D
scans.  A full multi-frequency 24-state three-dimensional grid is deliberately
not a CI workload; rate equations remain the appropriate weak-coherence grid
model.
