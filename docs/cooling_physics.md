# Cooling physics

The simulator exposes an effective Doppler/scattering model, multilevel
population rate equations, two-level and sparse multilevel OBE machinery, and
an adiabatic polarization-gradient population model. These are complementary
approximations, not chronological stages; see `model_hierarchy.md`.

The deterministic effective calculation is semiclassical: position and momentum
are classical, while photon momentum and the Lorentzian scattering rate retain
`hbar`.  The Monte Carlo solver uses the same effective rates but samples photon
events.  Neither mode is an experimentally calibrated model.

## Gaussian beams

For power `P` and 1/e² radius `w`, diffraction is negligible across the small
effective-model MOT volume and each travelling beam uses

```text
I(r_perp) = 2 P/(pi w²) exp(-2 |r_perp|²/w²),   s = I/I_sat.
```

The transverse integral of `I` is tested to equal `P`.  Six separate beam
objects are used; symmetry is a configuration, not an implementation shortcut.
Polarization helicity is defined relative to each beam's own propagation vector.
The local spherical components are projections onto q=-1,0,+1 about a supplied
quantization axis.  effective records this geometry but does not yet evolve Zeeman
populations from those components.

## Radiation pressure

For beam `i`, the effective detuning is

```text
delta_i = delta_L + delta_frequency - k_i dot v + delta_Z,i,
delta_Z,i = -sign(k_axis) sign(G_axis) mu_eff B_axis / hbar.
```

The shared-saturation scattering rate and mean optical force are

```text
R_i = (Gamma/2) s_i / [1 + sum_j s_j + (2 delta_i/Gamma)²],
F_opt = sum_i hbar k_i R_i,
F_total = F_opt + m g.
```

Shared saturation represents competition for one effective excited state.  The
signed Zeeman construction reproduces the restoring branch of an ideal MOT but
is **not** a multilevel Clebsch–Gordan calculation.  The configurable full vector
`g` permits arbitrary apparatus orientation.  Linear coefficients are numerical
central derivatives, `beta=-dF/dv` and `kappa=-dF/dx`; they are meaningful only
near a stable fixed point.

## Magnetic fields

The ideal field is

```text
B(r) = R diag(b', b', -2b') R^T (r-r0),
```

whose gradient has zero trace, as required by `div B=0`.  The physical coil
backend evaluates

```text
dB = mu0 I N/(4 pi) dl cross (r-r') / |r-r'|³
```

on midpoint straight-wire segments around independently positioned and oriented
circular loops.  Opposed currents form an anti-Helmholtz pair.  Radius,
separation, turns, currents, relative current imbalance, lateral displacement,
and relative tilt are explicit.  Segment refinement is tested away from the
wire.  Uniform, linear-gradient, and sinusoidal residual fields are composable.

## Deterministic and stochastic dynamics

The adaptive solver integrates `dr/dt=v`, `m dv/dt=F(r,v,t)` with RK45.  In the
photon-event solver an absorption from beam `i` adds `+hbar k_i`; the subsequent
spontaneous photon adds a recoil of fixed magnitude in an isotropically sampled
direction.  Its mean is zero and component variance is one third, both tested.
The step is rejected if the total scattering probability exceeds 0.1, because
the present Bernoulli scheme permits at most one event per atom per step.

Momentum diffusion from these random absorption/emission events is therefore
present in Monte Carlo runs.  No equilibrium temperature is asserted from short
effective-model trajectories.  In particular, the Doppler formula is not a sub-Doppler
prediction.

## Models deliberately not implemented in the effective model

A genuine F=1,2 / F'=0,1,2,3 calculation requires basis-resolved dipole matrix
elements, repump coupling, optical pumping and spontaneous branching.  OBEs add
a Hamiltonian, coherences and Lindblad decay.  Polarization-gradient cooling
additionally requires phase-resolved fields and spatial light shifts.  Empty
APIs for these models are intentionally absent.  Vapour pressure, capture and
loading/loss require loading/loss collision and calibrated geometry inputs and are
also not approximated by arbitrary curves.

## rate-equation: multilevel population rate equations

the multilevel model uses the complete 87Rb D2 hyperfine/Zeeman population basis

```text
5S1/2: F=1 (3 states), F=2 (5 states)
5P3/2: F'=0,1,2,3 (1+3+5+7 states), total=24.
```

Allowed electric-dipole couplings obey `m'=m+q`, `q∈{-1,0,+1}` and `ΔF=0,±1`
(excluding 0→0). Matrix elements are generated from

```text
|<F,m;1,q|F',m'>|²
```

using the SymPy Wigner 6-j and Clebsch–Gordan APIs, scaled to the documented hyperfine line
strengths and then to unity for the closed stretched F=2,m=2 → F'=3,m'=3
transition. This last normalization matches the circular cycling-transition
saturation-intensity convention. Spontaneous branching is generated from the
same squared dipole matrix elements and normalized independently for every
excited Zeeman state.

For travelling beam `b` and transition `g→e`, the stimulated rate is

```text
W_b,ge = (Γ/2) s_b C_ge² P_b(q) /
         [1 + (2 δ_b,ge/Γ)²],
δ_b,ge = δ_b + δ_hfs - k_b·v
         - (μB/ℏ)(g_e m_e - g_g m_g)|B|.
```

`P_b(q)` is the local spherical polarization fraction about `B/|B|`. At a field
zero the code uses a fixed z quantization axis; for a fully symmetric six-beam
configuration observables are basis invariant, but asymmetric zero-field cases
should be interpreted with care. Polarization impurity mixes the calculated
fractions with an isotropic one-third background.

Unlike effective, this expression has no inserted shared-saturation denominator:
stimulated absorption/emission and the finite ground/excited populations produce
saturation dynamically. Adding both mechanisms would double-count saturation.
Cooling and repump use the common D2 reference wavevector; their 6.8 GHz
frequency difference changes photon momentum by less than 2×10⁻⁵ and is
neglected at this fidelity.

The population vector obeys

```text
dp_e/dt = Σ_g,b W_b,ge (p_g-p_e) - Γ p_e,
dp_g/dt = Σ_e,b W_b,ge (p_e-p_g) + Γ Σ_e b_e→g p_e.
```

The stationary solution is the normalized null vector of this conservative
linear generator. The beam force retains net stimulated momentum transfer,

```text
F_b = ℏ k_b Σ_ge W_b,ge (p_g-p_e),
F = Σ_b F_b + m g.
```

Cooling beams address F=2→F'=3 and repump beams address F=1→F'=2, while their
off-resonant couplings to the other allowed excited hyperfine levels are kept.
Couplings to the other ground hyperfine manifold are neglected because their
6.835 GHz separation is far outside the configured linewidths.

### What rate-equation still cannot describe

The state vector contains populations, not a density matrix. Consequently the
model excludes ground/excited coherences, coherent population trapping, dark
states, Raman resonances, spatially varying light shifts, stimulated-force
interference and Sisyphus/polarization-gradient cooling. Those require coherent
OBEs and polarization-gradient phase-resolved fields. A rate-equation temperature is therefore not
reported from the steady-state force alone.

## coherent: reduced optical Bloch equations

the OBE benchmark introduces a coherence-resolving backend for one explicitly selected
transition, currently the effective 87Rb D2 stretched transition. In the
rotating frame and basis `{|g>,|e>}`,

```text
H/ℏ = [[0, Ω*/2],
       [Ω/2, -δ]].
```

Spontaneous decay is represented by `C=sqrt(Γ)|g><e|`, and the density matrix
obeys the Lindblad master equation

```text
dρ/dt = -i[H/ℏ,ρ] + CρC† - 1/2(C†Cρ + ρC†C).
```

The saturation convention is `s=2|Ω|²/Γ²`. The stationary excited population is

```text
ρee = (s/2) / [1+s+(2δ/Γ)²],
```

and a single travelling wave produces `F=ℏk Γρee`. The numerical Liouvillian
solution is tested against this analytical expression over detuning and
saturation, while adaptive time evolution tests coherent transients, trace,
Hermiticity, positivity and tolerance refinement.

This is a genuine OBE, but it is a **reduced two-state OBE**, not the eventual
24-state cooling+repump OBE. It cannot form Raman dark states or represent six
coherent beam phases. Its purpose is to establish a transparent Hamiltonian,
collapse-operator, observable and validation infrastructure before enlarging the
basis.

## Why damping can decrease as beam power increases

The damping coefficient is a derivative, `β=-∂F/∂v`, not the total scattering
force. At low saturation, more power increases scattering and steepens the
velocity response. At higher saturation, the transition population saturates
and power broadening flattens the Lorentzian slope at a fixed detuning. Thus β
has an optimum and can fall even while the maximum available radiation pressure
continues toward its saturation limit. Changing detuning moves that optimum;
a single power curve without its held-fixed detuning is incomplete.
The two-level OBE result therefore includes both a two-dimensional `(power, detuning)`
map and fixed-detuning cuts. The descending branch is not evidence that “more
light cools less” in every sense: the maximum scattering force can increase
while its small-velocity derivative decreases.

## Why force can decrease as Gaussian waist increases

For a Gaussian travelling beam,

```text
I0 = 2P/(πw²),   I(r_perp)=I0 exp(-2r_perp²/w²).
```

At fixed power, increasing `w` improves spatial overlap but reduces peak
intensity as `1/w²`; at small `w`, an off-axis atom can instead lie in the weak
Gaussian wing. Their competition naturally produces a non-monotonic optimum.
At fixed peak intensity, power must scale as `w²`, removing the direct intensity
penalty, although transverse beams can still change shared saturation at an
off-axis point. Therefore “force versus waist” has no unique meaning until
power/intensity, evaluation position, detuning and other beams are specified.

## Pure dephasing in the reduced OBE

Laser linewidth, collisions and technical phase noise can destroy optical
coherence without directly removing excited population. the OBE benchmark represents a
Markovian approximation with

```text
C_phi = sqrt(gamma_phi/2) sigma_z,
gamma_2 = Gamma/2 + gamma_phi.
```

The optical coherence then decays at `gamma_2`. The generalized stationary
excited population is

```text
rho_ee = |Omega|² gamma_2 /
         [2 Gamma (delta²+gamma_2²) + 2 |Omega|² gamma_2].
```

This is a controlled homogeneous-dephasing model, not a replacement for
explicit stochastic laser-frequency noise or a measured laser spectrum.

## polarization-gradient: phase-resolved polarization-gradient model

the polarization-gradient model selects the closed 87Rb D2 `F=2 -> F'=3` manifold. D2 is retained for
continuity with the MOT cycling transition; D1 gray molasses and 85Rb require
different hyperfine graphs and are not represented by changed constants. Each
beam contributes `E(r)=sum_b sqrt(s_b) epsilon_b exp[i(k_b.r+phi_b)]`.
Projection about the fixed quantization axis gives local `s_q(r)`. For
`m -> m'=m+q`, the squared coupling is `C²_mq s_q`. Low-saturation adiabatic
elimination gives

```text
U_m = sum_q hbar delta_mq Gamma² C²_mq s_q / [8(delta_mq²+Gamma²/4)],
R_mq = Gamma³ C²_mq s_q / [8(delta_mq²+Gamma²/4)].
```

The transition detuning includes the projected linear Zeeman shift. Excitation
and normalized spontaneous branching form a conservative five-state generator
`dp/dt=A(r)p`. Along `x=vt`, initial transients are discarded and
`F_x=-average[sum_m p_m partial_x U_m]` is evaluated over complete optical
periods. This is the Sisyphus correlation between climbing a dressed potential
and pumping into another Zeeman state.

The approximation neglects ground coherences, excited populations, nonadiabatic
corrections, recoil during the force solve, repump coupling, Gaussian envelopes
and uncontrolled phase noise. Transverse magnetic suppression requires
coherences; only axial Zeeman sensitivity is claimed. The recoil-diffusion
Einstein estimate is stored only for positive friction and is not an
experimental or complete equilibrium-temperature prediction.
