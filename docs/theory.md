# Phase-1 physical model

## Fidelity labels

The repository uses explicit labels rather than treating “MOT simulation” as a
single approximation.

| Level | Meaning | Status |
|---|---|---|
| A | Effective two-level Doppler/scattering model | Implemented |
| B | Multilevel rate equations with cooling and repump | Phase 2 |
| C | Multilevel optical Bloch equations (OBEs) | Phase 3 |
| D | Phase-resolved polarization-gradient/sub-Doppler dynamics | Phase 4 |
| E | Experimentally calibrated vapour loading/loss | Phase 5 |

The deterministic Level-A calculation is semiclassical: position and momentum
are classical, while photon momentum and the Lorentzian scattering rate retain
`hbar`.  The Monte Carlo solver uses the same Level-A rates but samples photon
events.  Neither mode is an experimentally calibrated model.

## Gaussian beams

For power `P` and 1/e² radius `w`, diffraction is negligible across the small
Phase-1 MOT volume and each travelling beam uses

```text
I(r_perp) = 2 P/(pi w²) exp(-2 |r_perp|²/w²),   s = I/I_sat.
```

The transverse integral of `I` is tested to equal `P`.  Six separate beam
objects are used; symmetry is a configuration, not an implementation shortcut.
Polarization helicity is defined relative to each beam's own propagation vector.
The local spherical components are projections onto q=-1,0,+1 about a supplied
quantization axis.  Level A records this geometry but does not yet evolve Zeeman
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
Phase-1 trajectories.  In particular, the Doppler formula is not a sub-Doppler
prediction.

## Models deliberately not implemented in Phase 1

A genuine F=1,2 / F'=0,1,2,3 calculation requires basis-resolved dipole matrix
elements, repump coupling, optical pumping and spontaneous branching.  OBEs add
a Hamiltonian, coherences and Lindblad decay.  Polarization-gradient cooling
additionally requires phase-resolved fields and spatial light shifts.  Empty
APIs for these models are intentionally absent.  Vapour pressure, capture and
loading/loss require Phase-5 collision and calibrated geometry inputs and are
also not approximated by arbitrary curves.

## Level B: multilevel population rate equations

Phase 2 uses the complete 87Rb D2 hyperfine/Zeeman population basis

```text
5S1/2: F=1 (3 states), F=2 (5 states)
5P3/2: F'=0,1,2,3 (1+3+5+7 states), total=24.
```

Allowed electric-dipole couplings obey `m'=m+q`, `q∈{-1,0,+1}` and `ΔF=0,±1`
(excluding 0→0). Matrix elements are generated from

```text
|<F,m;1,q|F',m'>|²
```

using the Racah factorial expression, scaled to the documented hyperfine line
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

Unlike Level A, this expression has no inserted shared-saturation denominator:
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

### What Level B still cannot describe

The state vector contains populations, not a density matrix. Consequently the
model excludes ground/excited coherences, coherent population trapping, dark
states, Raman resonances, spatially varying light shifts, stimulated-force
interference and Sisyphus/polarization-gradient cooling. Those require Level C
OBEs and Level D phase-resolved fields. A Level-B temperature is therefore not
reported from the steady-state force alone.
