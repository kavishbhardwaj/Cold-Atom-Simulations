# Part I — Physical system and apparatus

This chapter builds the physical inputs before any MOT force is calculated. For each layer, the tutorial gives the equation, the modelling decision, the implementation choice, the approximation, and the result that checks or illustrates that layer.

# 1. The physical problem

A magneto-optical trap cools and confines neutral atoms by combining:

1. near-resonant laser radiation pressure;
2. Doppler shifts, which make the force velocity dependent;
3. Zeeman shifts in a magnetic-field gradient, which make the force position dependent;
4. polarization selection rules and optical pumping among magnetic sublevels;
5. spontaneous-emission recoil, which produces diffusion and ultimately limits cooling;
6. in a real vapour-cell experiment, atomic loading, collisions, finite beam geometry, coil imperfections and time-dependent experimental controls.

A single equation is not enough to describe all of these effects at every useful fidelity. The central design decision of the repository was therefore **not to build one monolithic “exact MOT model.”** Instead, it implements a hierarchy of models. Fast models answer engineering questions cheaply; higher-fidelity models retain atomic structure and coherence but are computationally expensive.

The hierarchy is:

**atomic constants and basis**
→ **six laser beams + magnetic field**
→ **effective scattering-force MOT**
→ **multilevel population rate equations**
→ **optical Bloch equations (OBEs)**
→ **polarization-gradient / sub-Doppler models**
→ **classical trajectories and photon recoil**
→ **vapour capture and loading**
→ **experimental timing sequence**
→ **optional collective-cloud effects**.

The code deliberately keeps these layers separate so that an approximation used for speed is not silently presented as a high-fidelity prediction.

---
# 2. Why the reference system is 87Rb D2

The main reference system is

$$
^{87}\mathrm{Rb}:\quad 5S_{1/2}\rightarrow 5P_{3/2}\quad (D_2),
$$

with wavelength

$$
\lambda = 780.241209686\ \mathrm{nm}.
$$

The cooling laser is referenced to

$$
F=2\rightarrow F'=3,
$$

and the repump laser to

$$
F=1\rightarrow F'=2.
$$

This was chosen because the D2 line provides the familiar stretched cycling transition used in conventional rubidium MOTs, the isotope is widely used in cold-atom experiments, and authoritative atomic constants are available. More importantly, choosing one complete isotope/line makes the simulation auditable before attempting several partially supported systems.

The repository also stores atomic data for 85Rb and for both D1 lines. However, a D1 gray-molasses calculation is **not** produced by merely changing the wavelength. D1 gray molasses requires coherent Lambda-system/Raman dark-state physics, so it remains explicitly outside the currently validated model.

## Useful scales derived from the repository constants

For the 87Rb D2 lifetime

$$
\tau = 26.2348\ \mathrm{ns},
$$

the decay rate is

$$
\Gamma=\frac{1}{\tau}\approx 3.8117\times10^7\ \mathrm{s^{-1}},
$$

or

$$
\frac{\Gamma}{2\pi}\approx 6.0666\ \mathrm{MHz}.
$$

The two-level Doppler temperature scale is

$$
T_D=\frac{\hbar\Gamma}{2k_B}\approx 145.6\ \mu\mathrm K.
$$

The recoil velocity is

$$
v_r=\frac{\hbar k}{m}\approx 5.88\ \mathrm{mm/s},
$$

and the single-photon recoil temperature is approximately

$$
T_r=\frac{(\hbar k)^2}{2mk_B}\approx0.181\ \mu\mathrm K.
$$

These numbers are useful sanity checks throughout the code: a normal Doppler MOT temperature is expected around the 100-microkelvin scale, while genuine polarization-gradient cooling can go below the Doppler scale but remains well above a single recoil temperature in ordinary D2 molasses.

**Relevant code:** `src/cold_atom_mot/atomic/species.py`, `docs/atomic_systems.md`.

---
# 3. Atomic structure: from constants to a 24-state basis

## 3.1 Hyperfine energies

For a hyperfine level with nuclear spin $I$, electronic angular momentum $J$ and total angular momentum $F$, the code uses

$$
K=F(F+1)-I(I+1)-J(J+1).
$$

The magnetic-dipole hyperfine contribution is

$$
E_A=\frac{A}{2}K.
$$

When the electric-quadrupole constant $B$ is applicable, the code adds the standard rank-two hyperfine term. These energies generate the ground and excited hyperfine offsets rather than hard-coding a table of every Zeeman state.

For 87Rb D2 the coupled basis contains:

- ground $5S_{1/2}$: $F=1$ with 3 Zeeman states and $F=2$ with 5 states;
- excited $5P_{3/2}$: $F'=0,1,2,3$ with $1+3+5+7=16$ states.

Therefore

$$
N_\mathrm{states}=8+16=24.
$$

A full density matrix has

$$
24\times24=576
$$

complex entries before using Hermiticity and trace constraints. This is why the full OBE is used for representative points and short scans rather than huge three-dimensional grids.

## 3.2 Dipole selection rules and transition strengths

Electric-dipole couplings obey

$$
\Delta F=0,\pm1,\qquad F=0\not\leftrightarrow F'=0,
$$

and

$$
m_F'=m_F+q,\qquad q\in\{-1,0,+1\},
$$

where $q=-1,0,+1$ correspond to $\sigma^-$, $\pi$, and $\sigma^+$ components relative to the chosen quantization axis.

The code generates relative hyperfine line strengths using Wigner 6-j symbols and then the Zeeman-resolved strengths using Clebsch-Gordan coefficients. In schematic form,

$$
|\langle F,m;1,q|F',m'\rangle|^2
$$

is multiplied by the reduced hyperfine strength. The stretched cycling transition

$$
F=2,m_F=2\rightarrow F'=3,m_F'=3
$$

is normalized to unit strength under the cycling-transition saturation convention.

Spontaneous branching ratios are generated from the same dipole graph and normalized separately for each excited state. Thus the stimulated and spontaneous parts of the model share one consistent atomic basis.

### Why generate the basis instead of hard-coding it?

Because this makes isotope/line changes physical rather than cosmetic. 85Rb has a different nuclear spin and different hyperfine manifolds; D1 has a different excited-state $J$. Generating the graph from angular-momentum algebra prevents a D2 table from being accidentally reused as D1 data.

**Tool decision — use SymPy angular-momentum algebra.** The code uses SymPy Wigner 6-j and Clebsch–Gordan functions rather than a hand-entered transition table. This makes selection rules auditable, makes isotope changes reproducible, and removes a large class of transcription errors. NumPy arrays are then used for the numerical basis and branching matrices.

**Result of this step.** For 87Rb D2, the generated basis contains exactly 8 ground and 16 excited Zeeman states, giving the 24-state system used by the rate-equation and multilevel-OBE layers. The same generated dipole graph supplies both laser coupling strengths and spontaneous branching, so the two parts cannot silently use inconsistent transition tables.

---
# 4. Magnetic physics: from a MOT gradient to a full vector Zeeman Hamiltonian

## 4.1 Fast quadrupole field

The ideal MOT field is represented as

$$
\mathbf B(\mathbf r)
=R\,\mathrm{diag}(b',b',-2b')\,R^T(\mathbf r-\mathbf r_0).
$$

The trace of the gradient is zero:

$$
\nabla\cdot\mathbf B=0,
$$

as required by Maxwell's equations in the trapping region.

The reference effective configuration uses

$$
b'=0.10\ \mathrm{T/m}=10\ \mathrm{G/cm}.
$$

Therefore an axial displacement of 1 mm corresponds to a field scale of order 1 G on a radial axis in the ideal linear model.

## 4.2 Physical coils

For apparatus modelling, the code can replace the ideal gradient with segmented circular coils using the Biot-Savart law,

$$
d\mathbf B=\frac{\mu_0 I N}{4\pi}\frac{d\boldsymbol\ell\times(\mathbf r-\mathbf r')}{|\mathbf r-\mathbf r'|^3}.
$$

Anti-Helmholtz pairs generate the MOT gradient. Helmholtz-like pairs provide three-axis bias compensation. Radius error, separation error, tilt, displacement, turn mismatch and current imbalance can be introduced explicitly.

A measured calibration can be represented as

$$
\mathbf B=M\mathbf I+\mathbf B_\mathrm{offset},
$$

where $M$ is a measured 3x3 coil-response matrix. The compensation currents are obtained by least squares, allowing cross-axis coupling instead of assuming three perfectly orthogonal coils.

## 4.3 Full vector hyperfine-Zeeman Hamiltonian

For coherent calculations, a scalar $g_Fm_FB$ shift is insufficient when the field is transverse or when the quantization axis changes. The repository therefore constructs

$$
H=H_\mathrm{hfs}+H_Z,
$$

with

$$
H_Z=\mu_B\left(g_J\mathbf J+g_I\mathbf I\right)\cdot\mathbf B.
$$

The Hamiltonian is constructed in a field-independent uncoupled basis $|m_I,m_J\rangle$, then transformed to the coupled $|F,m_F\rangle$ ordering used by the optical dipole graph. This avoids a discontinuous “follow the local field” basis when $\mathbf B\rightarrow0$.

At weak field it recovers

$$
\Delta E\approx g_F\mu_Bm_FB,
$$

while retaining transverse mixing and nonlinear Zeeman shifts within the selected fine-structure manifold.

### Approximation

Coupling to other fine-structure/electronic manifolds and diamagnetic terms is omitted. The model therefore describes hyperfine decoupling within the selected D line, not an unrestricted high-field Paschen-Back calculation.

### Validation result

The 87Rb ground hyperfine-Zeeman spectrum was compared with PyLCP under matched conventions. The maximum reported spectral difference over the tested fields was approximately

$$
0.57\ \mathrm{Hz}.
$$

This is an independent validation of the atomic magnetic Hamiltonian, not yet a validation of the complete light-force calculation.

**Why this result matters.** The vector Zeeman Hamiltonian is one of the foundations used later for transverse stray fields and the 24-state OBE. Agreement at the sub-hertz level with an independent PyLCP construction means later disagreement is much less likely to come from the basic hyperfine/Zeeman spectrum itself.

![Exact 87Rb vector-Zeeman spectra](../../results/atomic_structure/exact_zeeman_spectra.png)

![Linear approximation compared with the exact vector-Zeeman calculation](../../results/atomic_structure/linear_vs_exact_zeeman.png)

---
# 5. The six laser beams: why the code does not use one “MOT intensity”

A real three-dimensional MOT contains three counterpropagating beam pairs. The repository therefore represents **six physical beams**, each of which can have its own:

- power;
- direction and origin;
- waist and ellipticity;
- focus;
- detuning and AOM offset;
- optical phase;
- linewidth;
- coherence group;
- Jones vector and optical-element train;
- pointing error and retroreflection loss.

## 5.1 Gaussian intensity

For a collimated circular Gaussian beam,

$$
I(r_\perp)=\frac{2P}{\pi w^2}\exp\left(-\frac{2r_\perp^2}{w^2}\right).
$$

The on-axis saturation parameter is

$$
s=\frac{I}{I_\mathrm{sat}}.
$$

The reference effective MOT uses

$$
P=10\ \mathrm{mW},\qquad w=8\ \mathrm{mm}.
$$

Hence

$$
I_0=\frac{2(0.010)}{\pi(0.008)^2}\approx99.47\ \mathrm{W/m^2}.
$$

With the stored 87Rb D2 saturation intensity $I_\mathrm{sat}=16.69\ \mathrm{W/m^2}$, the peak cycling-transition saturation per beam is approximately

$$
s_0\approx5.96.
$$

This is a fairly strong illustrative MOT configuration, which is why saturation and power broadening cannot be ignored in the effective force.

## 5.2 Polarization: why “three sigma+ and three sigma- beams” is incomplete

$\sigma^+$, $\pi$ and $\sigma^-$ are defined relative to a quantization axis. A circularly polarized beam viewed along its own propagation direction need not remain pure $\sigma^+$ relative to a different laboratory or local magnetic-field axis.

The code therefore stores the complex polarization/Jones vector and computes the local spherical fractions

$$
P_{-1}+P_0+P_{+1}=1.
$$

This is essential for multilevel optical pumping and for understanding polarization gradients.

## 5.3 Coherence groups

Not every beam in a real MOT is mutually phase coherent. The OBE therefore distinguishes:

- six independent beams;
- three coherent counterpropagating pairs;
- a fully coherent six-beam field.

Beams in one coherence group retain their relative phases. Mutually incoherent groups are averaged by deterministic phase cycling with convergence refinement. This prevents unrelated lasers from being silently treated as one phase-stable standing-wave field.

### Apparatus sensitivity result

For the specified effective-MOT reference recipe, changing the power of one x beam by +/-10% moved the solved MOT centre by approximately

$$
-0.74\ \mathrm{mm}\quad\mathrm{to}\quad +0.69\ \mathrm{mm}.
$$

A 5 mrad pointing perturbation produced a displacement of about

$$
36\ \mu\mathrm m.
$$

These are recipe-specific model sensitivities, not universal laboratory tolerances.

**Result and interpretation.** The six-beam layer is therefore not decorative geometry: it converts specific optical imperfections into a calculable force imbalance and MOT-centre displacement. The result also explains why the repository keeps six physical beam objects instead of one global “MOT intensity” parameter.

![Six physical MOT beams](../../results/laser_apparatus/six_beam_apparatus.png)

![Controlled beam-imperfection study](../../results/laser_apparatus/apparatus_imperfections.png)

---
