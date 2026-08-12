# Part II — MOT force models, OBEs, and sub-Doppler physics

# 6. First MOT model: effective semiclassical radiation pressure

The first force model intentionally sacrifices detailed hyperfine structure for speed and interpretability.

For beam \(i\), the effective detuning is

\[
\delta_i=\delta_L+\delta_{\mathrm{offset}}-\mathbf k_i\cdot\mathbf v+\delta_{Z,i}.
\]

The \(-\mathbf k_i\cdot\mathbf v\) term is the Doppler shift. The approximate Zeeman term changes the resonance oppositely for the two counterpropagating beams so that an ideal MOT becomes restoring.

The shared-saturation scattering rate is

\[
R_i=\frac{\Gamma}{2}
\frac{s_i}
{1+\sum_j s_j+\left(2\delta_i/\Gamma\right)^2}.
\]

The mean optical force is

\[
\mathbf F_\mathrm{opt}=\sum_i\hbar\mathbf k_iR_i,
\]

and the total force includes gravity,

\[
\mathbf F=\mathbf F_\mathrm{opt}+m\mathbf g.
\]

Near a stable point, the force can be linearized:

\[
F_x\approx-\kappa x-\beta v_x,
\]

where

\[
\kappa=-\left.\frac{\partial F_x}{\partial x}\right|_0,
\qquad
\beta=-\left.\frac{\partial F_x}{\partial v_x}\right|_0.
\]

### Why this model exists

It is fast enough for trajectories, capture scans, beam-waist studies and apparatus sensitivity calculations. It is also transparent: if a force changes, one can usually trace the change to detuning, intensity, Doppler shift or Zeeman shift.

### Main approximations

- one effective cooling transition;
- no explicit repump population dynamics;
- no Zeeman-state optical pumping;
- no coherences or dark states;
- no true polarization-gradient cooling;
- shared saturation is inserted phenomenologically for the effective cycling transition.

Therefore this model is suitable for **MOT-scale force and capture studies**, but not for claiming a quantitative sub-Doppler temperature.

---
# 7. Multilevel rate equations: adding the real hyperfine population graph

The next model retains all 24 87Rb D2 hyperfine-Zeeman populations but still neglects coherence.

For beam \(b\) and transition \(g\rightarrow e\), the stimulated rate has the structure

\[
W_{b,ge}=\frac{\Gamma}{2}
\frac{s_b\,C_{ge}^2\,P_b(q)}
{1+\left(2\delta_{b,ge}/\Gamma\right)^2}.
\]

Here:

- \(C_{ge}^2\) is the generated dipole strength;
- \(P_b(q)\) is the local \(\sigma^-\), \(\pi\) or \(\sigma^+\) fraction;
- \(\delta_{b,ge}\) includes laser detuning, hyperfine offset, Doppler shift and the transition Zeeman shift.

The population equations are

\[
\frac{dp_e}{dt}
=\sum_{g,b}W_{b,ge}(p_g-p_e)-\Gamma p_e,
\]

\[
\frac{dp_g}{dt}
=\sum_{e,b}W_{b,ge}(p_e-p_g)
+\Gamma\sum_e b_{e\rightarrow g}p_e.
\]

The stationary population is the normalized null vector of the rate generator.

The beam force is

\[
\mathbf F_b=\hbar\mathbf k_b
\sum_{ge}W_{b,ge}(p_g-p_e),
\]

and

\[
\mathbf F=\sum_b\mathbf F_b+m\mathbf g.
\]

## Why there is no shared-saturation denominator here

Unlike the effective model, saturation arises through the finite ground/excited populations and stimulated absorption/emission. Adding the effective-model shared denominator would double-count saturation.

## Cooling and repump

The reference configuration uses six cooling beams near \(F=2\rightarrow F'=3\) and six repump beams near \(F=1\rightarrow F'=2\). Off-resonant coupling to other allowed excited hyperfine levels is retained.

Coupling of a cooling carrier to the other ground hyperfine manifold is neglected because the ground-state splitting is several GHz, far beyond the optical linewidth in the reference regime. The higher-fidelity OBE exposes a quantitative bound on this rotating-wave approximation.

## What the rate equations still cannot do

A population vector cannot represent phase coherence. Therefore this model excludes:

- coherent population trapping;
- Raman dark states;
- ground-state coherences;
- standing-wave interference force;
- genuine Sisyphus cooling.

It is the correct intermediate model when optical pumping matters but coherence can be neglected.

---
# 8. Optical Bloch equations: the coherence-resolving layer

## 8.1 Two-level OBE benchmark

Before building the 24-state OBE, the repository validates the machinery on a two-level atom. In a rotating frame,

\[
\frac{H}{\hbar}=
\begin{pmatrix}
0 & \Omega^*/2\\
\Omega/2 & -\delta
\end{pmatrix}.
\]

Spontaneous decay uses

\[
C=\sqrt{\Gamma}\,|g\rangle\langle e|,
\]

and the Lindblad master equation is

\[
\dot\rho=-\frac{i}{\hbar}[H,\rho]
+C\rho C^\dagger
-\frac12\left(C^\dagger C\rho+\rho C^\dagger C\right).
\]

With

\[
s=\frac{2|\Omega|^2}{\Gamma^2},
\]

the stationary excited population is

\[
\rho_{ee}=\frac{s/2}{1+s+(2\delta/\Gamma)^2}.
\]

A travelling wave then gives

\[
F=\hbar k\Gamma\rho_{ee}.
\]

### Validation result

The two-level steady-state population agrees with QuTiP to a maximum absolute difference of approximately

\[
5.55\times10^{-17},
\]

and the Liouvillian/dynamical tests also agree within the committed numerical tolerances. The normalized two-beam Doppler force agrees with PyLCP to a maximum relative difference of about

\[
7.9\times10^{-15}.
\]

These results are important because they validate conventions for \(\Gamma\), detuning, Rabi frequency, saturation and force before moving to a much larger basis.

## 8.2 The 24-state moving-atom OBE

The high-fidelity 87Rb D2 solver uses the 24-state density matrix and the full generated spontaneous branching graph.

For a moving atom,

\[
\mathbf r(t)=\mathbf r_0+\mathbf vt.
\]

Each beam has a phase

\[
\phi_i(t)=\mathbf k_i\cdot(\mathbf r_0+\mathbf vt)
-\delta\omega_i t+\phi_{i,0}.
\]

Thus every physical beam carries its own Doppler shift \(\mathbf k_i\cdot\mathbf v\).

### Why a block-rotating frame was introduced

Cooling and repump optical frequencies differ by several GHz. Explicitly integrating the corresponding sub-nanosecond beat would be numerically wasteful because the experimentally relevant internal/mechanical dynamics are much slower.

The solver therefore uses one excited-manifold reference and independent carrier rotations for the two ground hyperfine manifolds. The large cooling-repump carrier separation is removed analytically. Only physically relevant residual same-manifold AOM/frequency/Doppler beats remain time dependent.

### Controlled rotating-wave approximations

Two approximations remain explicit:

1. one optical carrier addresses its selected ground hyperfine manifold, while the far-off-resonant other-ground-manifold drive is discarded;
2. magnetic matrix elements connecting different ground-F blocks oscillate at the ground hyperfine splitting in the block frame and are secularly discarded.

The code evaluates quantitative diagnostics for the optical approximation using the smallest actual discarded transition detuning. In the reference MOT, the reported population-scale bound is of order a few parts in \(10^6\), making the approximation controlled in the weak-field MOT regime.

### Force from the Hamiltonian gradient

Instead of estimating force only from total excited-state decay, the OBE uses

\[
\mathbf F_i
=-\hbar\,\mathrm{Tr}\left(\rho\,\nabla\frac{H_i}{\hbar}\right).
\]

The analytic gradient contains:

- \(i\mathbf k_i\) from the travelling-wave phase;
- the Gaussian amplitude gradient \(-2\mathbf r_\perp/w^2\).

This retains radiation pressure and coherent/dipole-force contributions allowed by the model.

### Incoherent beams

When beam groups are mutually incoherent, one fixed arbitrary phase realization is not physical. The solver therefore averages observables over deterministic phase ensembles and refines the ensemble until the result converges. Coherent beams retain their relative phase.

### What remains unvalidated

The full 24-state moving OBE is internally tested but its complete 87Rb multilevel populations and force have not yet been independently matched against a full PyLCP 87Rb calculation. This is one of the most important remaining validation steps.

---
# 9. Polarization-gradient cooling: why sub-Doppler physics needs a different model

A standard Doppler force can cool below the initial thermal speed, but it cannot explain the lowest temperatures observed in multilevel alkali molasses. Polarization-gradient cooling depends on spatially varying light shifts and optical pumping among magnetic substates.

The repository therefore includes a reduced, phase-resolved 87Rb D2

\[
F=2\rightarrow F'=3
\]

population model.

Each beam contributes a complex field amplitude schematically as

\[
\mathbf E(\mathbf r)\propto
\sum_b\sqrt{s_b}\,\boldsymbol\epsilon_b
\exp[i(\mathbf k_b\cdot\mathbf r+\phi_b)].
\]

Projection onto local spherical components produces \(s_q(\mathbf r)\). For a ground state \(m\) and polarization \(q\), low-saturation adiabatic elimination gives a state-dependent light shift of the form

\[
U_m(\mathbf r)
=\sum_q
\frac{\hbar\,\delta_{mq}\,\Gamma^2
C_{mq}^2s_q(\mathbf r)}
{8\left(\delta_{mq}^2+\Gamma^2/4\right)},
\]

and an optical-pumping rate

\[
R_{mq}(\mathbf r)
=\frac{\Gamma^3 C_{mq}^2s_q(\mathbf r)}
{8\left(\delta_{mq}^2+\Gamma^2/4\right)}.
\]

The five \(F=2\) ground populations obey a position-dependent rate equation as the atom moves through the optical lattice. The Sisyphus force is evaluated from

\[
F_x=-\left\langle
\sum_m p_m\frac{\partial U_m}{\partial x}
\right\rangle.
\]

The physical picture is: atoms climb a state-dependent optical potential, are optically pumped near the top into another sublevel, and lose kinetic energy repeatedly.

## Reference PGC parameters

The committed reduced-model reference uses

\[
\Delta=-3\Gamma,
\qquad s=0.08\ \text{per beam}.
\]

It integrates over 24 optical periods, discarding the first 12 as transient.

## Important approximation boundary

This population PGC model omits ground-state coherence. It therefore cannot describe the full magnetic-field suppression of PGC, Raman dark states, or a complete force-noise diffusion tensor. It is a pedagogical/mechanistic Sisyphus model, not a final quantitative temperature solver.

---
# 10. Stray magnetic fields and sub-Doppler cooling

The repository explicitly avoids converting an incomplete diffusion model into a false temperature prediction.

A full-vector 24-state OBE diagnostic scans residual fields along x, y and z. For the configured

\[
\Delta=-3\Gamma,\qquad s=0.08/\text{beam}
\]

reference, two useful internal timescales are:

\[
\frac{\omega_L}{2\pi}\approx699.6\ \mathrm{Hz/mG}
\]

for the \(F=2\) Larmor scale, and a weak-drive optical-pumping scale of about

\[
6.56\ \mathrm{kHz}.
\]

These simple rates become comparable around

\[
B\approx9.4\ \mathrm{mG}.
\]

### How to interpret 9.4 mG

It means magnetic precession has reached the same order of timescale as the optical pumping that establishes the sub-Doppler internal-state distribution. It is therefore a physically meaningful warning scale.

It is **not**:

- a predicted 10% temperature threshold;
- a measured compensation requirement;
- a universal MOT number.

A defensible temperature curve \(T(B)\) requires both:

\[
\beta(B)=-\left.\frac{\partial F}{\partial v}\right|_{v=0}
\]

from converged moving-lattice OBE calculations, and a momentum-diffusion tensor at compatible fidelity. Internal-state switching and dipole-force fluctuations are not yet fully included in the diffusion calculation, so the repository does not manufacture \(T(B)\) from mismatched models.

This is an example of a deliberate scientific decision: **it is better to report a validated timescale marker than an attractive but unjustified temperature number.**

---
