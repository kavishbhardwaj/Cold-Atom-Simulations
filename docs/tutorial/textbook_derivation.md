# A Worked Textbook Derivation of the `cold-atom-simulations` MOT Model

This chapter is written as a worked derivation. The goal is not merely to tell you **what the code calls a solver**, but to show how the physical model is assembled mathematically from one rubidium atom all the way to a MOT force, a 24-state density matrix, trajectories, loading, and the present sub-Doppler model.

The repeated structure is:

> **define the physical quantity → derive the equation → state the approximation → show the exact numerical object solved by the code → show the result → explain why the next equation is needed.**

The reference system is \(^{87}\mathrm{Rb}\) on the D2 line. SI units are used internally. Spectroscopic constants tabulated in Hz are explicitly converted to angular frequency when they enter Hamiltonians.

---

# 1. What a MOT calculation must ultimately produce

A magneto-optical trap must both cool and confine. Near the centre, a useful local target is

$$
F_x(x,v_x)\approx-\kappa x-\beta_v v_x,
$$

where

- \(x\) is displacement from the magnetic-field zero;
- \(v_x\) is atomic velocity;
- \(\kappa>0\) is the local spring constant;
- \(\beta_v>0\) is the velocity-damping coefficient.

The code does **not** assume this equation globally. It calculates the full force and uses finite differences near the origin to obtain

$$
\kappa=-\left.\frac{\partial F_x}{\partial x}\right|_{0},
\qquad
\beta_v=-\left.\frac{\partial F_x}{\partial v_x}\right|_{0}.
$$

To calculate \(F_x\), however, we first need the atomic energy levels, transition strengths, laser fields, and magnetic field.

---

# 2. Atomic constants and natural scales

For the reference \(^{87}\mathrm{Rb}\) D2 line,

$$
5S_{1/2}\rightarrow5P_{3/2},
\qquad
\lambda=780.241209686\ \mathrm{nm}.
$$

Define

$$
k=\frac{2\pi}{\lambda},
$$

where \(k\) is the optical wave number. The excited-state lifetime is

$$
\tau=26.2348\ \mathrm{ns},
$$

so the spontaneous population-decay rate is

$$
\Gamma=\frac{1}{\tau}
\approx3.8117\times10^7\ \mathrm{s^{-1}},
$$

or

$$
\frac{\Gamma}{2\pi}\approx6.0666\ \mathrm{MHz}.
$$

Two useful temperature scales are

$$
T_D=\frac{\hbar\Gamma}{2k_B}
\approx145.6\ \mu\mathrm K,
$$

and

$$
T_r=\frac{(\hbar k)^2}{2mk_B}
\approx0.181\ \mu\mathrm K.
$$

The recoil velocity is

$$
v_r=\frac{\hbar k}{m}\approx5.88\ \mathrm{mm/s}.
$$

These scales are used as sanity checks. A temperature of \(50\ \mu\mathrm K\) is sub-Doppler relative to \(T_D\), but is still hundreds of recoil temperatures.

---

# 3. Hyperfine structure from angular momentum

## 3.1 Define the quantum numbers

For one fine-structure manifold:

- \(I\) = nuclear spin;
- \(J\) = total electronic angular momentum;
- \(F\) = total hyperfine angular momentum, with

$$
\mathbf F=\mathbf I+\mathbf J;
$$

- \(m_F\) = projection of \(F\) on the chosen laboratory quantization axis.

For \(^{87}\mathrm{Rb}\),

$$
I=\frac32.
$$

The allowed \(F\) values are

$$
F=|I-J|,\ |I-J|+1,\ldots,I+J.
$$

Therefore:

- \(5S_{1/2}\), \(J=1/2\Rightarrow F=1,2\);
- \(5P_{3/2}\), \(J=3/2\Rightarrow F'=0,1,2,3\).

## 3.2 Magnetic-dipole hyperfine interaction

Use

$$
\mathbf I\cdot\mathbf J
=\frac12\left[F(F+1)-I(I+1)-J(J+1)\right].
$$

Define the shorthand

$$
K\equiv F(F+1)-I(I+1)-J(J+1).
$$

The symbol \(A_{\rm hfs}\) is the **magnetic-dipole hyperfine constant**. It is a measured atomic constant tabulated in Hz. The magnetic-dipole energy is

$$
\frac{E_A}{h}=\frac{A_{\rm hfs}}{2}K.
$$

For \(^{87}\mathrm{Rb}\) \(5S_{1/2}\), the stored value is

$$
A_{\rm hfs}=3.417341305452145\ \mathrm{GHz}.
$$

For \(F=2\), \(K=1.5\). For \(F=1\), \(K=-2.5\). Hence

$$
\frac{E_{F=2}-E_{F=1}}{h}
=\frac{A_{\rm hfs}}{2}(1.5-(-2.5))
=2A_{\rm hfs}
\approx6.83468\ \mathrm{GHz}.
$$

That number later explains why cooling and repump carriers are separated by several GHz.

## 3.3 Electric-quadrupole term

The symbol \(B_{\rm hfs}\) is the **electric-quadrupole hyperfine constant**. It is not the magnetic field. The magnetic field is always written as the vector \(\mathbf B\).

When both \(I\ge1\) and \(J\ge1\), the code uses

$$
\frac{E_{\rm hfs}}{h}
=
\frac{A_{\rm hfs}}{2}K
+B_{\rm hfs}
\frac{\frac34K(K+1)-I(I+1)J(J+1)}
{2I(2I-1)J(2J-1)}.
$$

For the \(5S_{1/2}\) ground state, \(J=1/2\), so the quadrupole term vanishes. For \(5P_{3/2}\), it is retained.

![Hyperfine levels generated from the stored constants](figures/hyperfine_energy_levels.svg)

---

# 4. Construct the 24-state basis explicitly

The code orders states by increasing \(F\), then increasing \(m_F\).

## 4.1 Ground basis: 8 states

$$
\begin{aligned}
|g_1\rangle&=|F=1,m_F=-1\rangle,\\
|g_2\rangle&=|1,0\rangle,\\
|g_3\rangle&=|1,+1\rangle,\\
|g_4\rangle&=|2,-2\rangle,\\
|g_5\rangle&=|2,-1\rangle,\\
|g_6\rangle&=|2,0\rangle,\\
|g_7\rangle&=|2,+1\rangle,\\
|g_8\rangle&=|2,+2\rangle.
\end{aligned}
$$

## 4.2 Excited basis: 16 states

$$
\begin{aligned}
|e_1\rangle&=|F'=0,m_F'=0\rangle,\\
|e_2\rangle&=|1,-1\rangle,
&|e_3\rangle&=|1,0\rangle,
&|e_4\rangle&=|1,+1\rangle,\\
|e_5\rangle&=|2,-2\rangle,
&|e_6\rangle&=|2,-1\rangle,
&|e_7\rangle&=|2,0\rangle,\\
|e_8\rangle&=|2,+1\rangle,
&|e_9\rangle&=|2,+2\rangle,\\
|e_{10}\rangle&=|3,-3\rangle,
&|e_{11}\rangle&=|3,-2\rangle,
&|e_{12}\rangle&=|3,-1\rangle,\\
|e_{13}\rangle&=|3,0\rangle,
&|e_{14}\rangle&=|3,+1\rangle,
&|e_{15}\rangle&=|3,+2\rangle,
&|e_{16}\rangle&=|3,+3\rangle.
\end{aligned}
$$

The complete basis is

$$
\{|1\rangle,\ldots,|24\rangle\}
=
\{|g_1\rangle,\ldots,|g_8\rangle,
|e_1\rangle,\ldots,|e_{16}\rangle\}.
$$

This gives

$$
N=8+16=24.
$$

---

# 5. Dipole-allowed transitions and their strengths

Electric-dipole selection rules are

$$
\Delta F=0,\pm1,
\qquad
F=0\not\leftrightarrow F'=0,
$$

and

$$
\Delta m_F=q\in\{-1,0,+1\}.
$$

Here

- \(q=-1\Rightarrow\sigma^-\);
- \(q=0\Rightarrow\pi\);
- \(q=+1\Rightarrow\sigma^+\).

The relative hyperfine reduced strength is generated from a Wigner 6-j symbol:

$$
S_{F\rightarrow F'}
\propto
(2F'+1)(2J_g+1)
\begin{Bmatrix}
J_e & F' & I\\
F & J_g & 1
\end{Bmatrix}^{2}.
$$

For one Zeeman-resolved transition,

$$
S_{Fm\rightarrow F'm'}
\propto
S_{F\rightarrow F'}
\left|
\langle F,m;1,q|F',m'\rangle
\right|^2.
$$

The code generates every allowed transition using SymPy Wigner and Clebsch-Gordan functions, then normalizes the strongest transition to unit strength.

Spontaneous branching from one excited state is obtained from the same strength graph:

$$
b_{e\rightarrow g}
=
\frac{S_{ge}}
{\sum_{g'}S_{g'e}}.
$$

Thus

$$
\sum_g b_{e\rightarrow g}=1
$$

for every excited state.

**Decision.** The transition graph is generated rather than hard-coded, so changing isotope or D line changes the angular-momentum algebra rather than merely changing a wavelength.

---

# 6. Magnetic Hamiltonian

## 6.1 Weak-field Landé factor

For a coupled hyperfine state,

$$
g_F=
\frac{
 g_J[F(F+1)+J(J+1)-I(I+1)]
+g_I[F(F+1)+I(I+1)-J(J+1)]
}{2F(F+1)}.
$$

The familiar weak-field shift is

$$
\Delta E\approx g_F\mu_Bm_FB.
$$

This is useful for intuition but insufficient for arbitrary vector fields.

## 6.2 Exact operator used by the code

The code constructs the hyperfine Hamiltonian first in the uncoupled basis

$$
|m_I,m_J\rangle.
$$

In operator form,

$$
\frac{H_{\rm hfs}}{h}
=A_{\rm hfs}\,\mathbf I\cdot\mathbf J
+B_{\rm hfs}
\frac{
3(\mathbf I\cdot\mathbf J)^2
+\frac32\mathbf I\cdot\mathbf J
-I(I+1)J(J+1)\mathbf 1
}
{2I(2I-1)J(2J-1)}.
$$

The exact vector Zeeman interaction is

$$
H_Z
=\mu_B
\left(g_J\mathbf J+g_I\mathbf I\right)\cdot\mathbf B.
$$

Therefore

$$
H_{\rm atom}=H_{\rm hfs}+H_Z.
$$

The coupled and uncoupled bases are related by Clebsch-Gordan coefficients:

$$
|F,m_F\rangle
=\sum_{m_I,m_J}
\langle I,m_I;J,m_J|F,m_F\rangle
|m_I,m_J\rangle.
$$

The matrix is transformed into the fixed laboratory \( |F,m_F\rangle \) ordering used by the optical dipole graph.

**Decision.** A fixed basis is used rather than choosing the quantization axis to follow the local \(\mathbf B\). This avoids a basis singularity when \(\mathbf B\rightarrow0\) and allows transverse magnetic mixing to appear as off-diagonal Hamiltonian elements.

![Exact vector Zeeman spectra](../../results/atomic_structure/exact_zeeman_spectra.png)

The calculated ground-state spectrum agrees with PyLCP within about \(0.57\) Hz over the tested field range.

---

# 7. The six physical laser beams

For a collimated elliptical Gaussian beam,

$$
I(x,y)=
\frac{2P}{\pi w_xw_y}
\exp\left[-2\left(\frac{x^2}{w_x^2}+\frac{y^2}{w_y^2}\right)\right].
$$

For a propagating Gaussian beam,

$$
w_x(z)=w_{x0}\sqrt{1+\left(\frac{z}{z_{Rx}}\right)^2},
$$

$$
w_y(z)=w_{y0}\sqrt{1+\left(\frac{z}{z_{Ry}}\right)^2},
$$

with Rayleigh ranges

$$
z_{Rx}=\frac{\pi w_{x0}^2}{\lambda},
\qquad
z_{Ry}=\frac{\pi w_{y0}^2}{\lambda}.
$$

The saturation parameter is

$$
s_i(\mathbf r)=\frac{I_i(\mathbf r)}{I_{\rm sat}}.
$$

The complex optical field is represented schematically as

$$
\mathbf E_i(\mathbf r)
\propto
\sqrt{I_i(\mathbf r)}\,
\boldsymbol\epsilon_i
\exp[i\Phi_i(\mathbf r)],
$$

where \(\boldsymbol\epsilon_i\) is the complex polarization vector.

The wave vector is

$$
\mathbf k_i=\frac{2\pi}{\lambda}\hat{\mathbf k}_i.
$$

A circular polarization relative to the propagation direction is constructed from a transverse basis \(\mathbf e_1,\mathbf e_2\):

$$
\boldsymbol\epsilon_{\pm}
=\frac{\mathbf e_1\pm i\mathbf e_2}{\sqrt2}.
$$

Relative to a fixed quantization axis, any polarization is decomposed into spherical components:

$$
P_q=|\boldsymbol\epsilon_q^\dagger\boldsymbol\epsilon|^2,
\qquad
\sum_{q=-1}^{+1}P_q=1.
$$

This is why the code does not treat “\(\sigma^+\)” as a global beam label independent of coordinate system.

![Six physical MOT beams](../../results/laser_apparatus/six_beam_apparatus.png)

---

# 8. First reduction: effective two-level MOT force

Before solving 24 internal states, the code uses a fast effective force for large trajectories and capture ensembles.

For beam \(i\), define

$$
\delta_i
=\delta_{L,i}+\delta_{{\rm AOM},i}
-\mathbf k_i\cdot\mathbf v
+\delta_{Z,i}.
$$

The terms are, respectively:

1. laser detuning;
2. configured AOM/frequency offset;
3. Doppler shift;
4. effective Zeeman shift.

If the beam has angular linewidth \(\gamma_{L,i}\), define

$$
\Gamma_i^{\rm eff}=\Gamma+\gamma_{L,i}.
$$

The code uses the shared-saturation scattering rate

$$
R_i=
\frac{\Gamma}{2}
\frac{s_i\Gamma/\Gamma_i^{\rm eff}}
{1+\sum_js_j+(2\delta_i/\Gamma_i^{\rm eff})^2}.
$$

Each absorption gives mean momentum \(\hbar\mathbf k_i\), hence

$$
\mathbf F_{\rm opt}
=\sum_i\hbar\mathbf k_iR_i.
$$

Gravity is added:

$$
\mathbf F=\mathbf F_{\rm opt}+m\mathbf g.
$$

![Two-beam Doppler damping](figures/doppler_force_vs_velocity.svg)

![MOT restoring-force illustration](figures/mot_restoring_force.svg)

![Full effective-MOT force map](../../results/effective_mot/force_map_x_vx.png)

**Approximation.** This model has no explicit \(F,m_F\) populations and no quantum coherence. It exists because it is fast enough for large mechanical ensembles.

---

# 9. Second reduction: 24-population rate equations

Keep all 24 populations but set all coherences to zero.

For beam \(b\) and transition \(g\leftrightarrow e\), define

$$
s_{b,ge}^{\rm eff}
=s_bS_{ge}P_b(q).
$$

The stimulated transition rate is

$$
W_{b,ge}
=\frac{\Gamma}{2}
\frac{s_{b,ge}^{\rm eff}\Gamma/\Gamma_b^{\rm eff}}
{1+(2\delta_{b,ge}/\Gamma_b^{\rm eff})^2}.
$$

Ground and excited populations obey

$$
\dot p_e
=\sum_{g,b}W_{b,ge}(p_g-p_e)-\Gamma p_e,
$$

$$
\dot p_g
=\sum_{e,b}W_{b,ge}(p_e-p_g)
+\Gamma\sum_e b_{e\rightarrow g}p_e.
$$

This can be written

$$
\dot{\mathbf p}=A_{\rm rate}\mathbf p,
$$

where \(A_{\rm rate}\) is the population generator matrix. It is not the hyperfine constant \(A_{\rm hfs}\).

The stationary state satisfies

$$
A_{\rm rate}\mathbf p_{\rm ss}=0,
$$

with

$$
\sum_i p_i=1.
$$

The force from beam \(b\) is

$$
\mathbf F_b
=\hbar\mathbf k_b
\sum_{ge}W_{b,ge}(p_g-p_e).
$$

![Effective versus multilevel force](../../results/multilevel/effective_vs_multilevel_force.png)

![Manifold populations](../../results/multilevel/manifold_populations.png)

**Why this is not enough.** A population model only stores diagonal density-matrix elements. It cannot describe phase coherence, coherent population trapping, Raman dark states, standing-wave interference, or the full Sisyphus mechanism.

---

# 10. The density matrix: from two states to 24 states

## 10.1 Two-level example first

For states \(|g\rangle,|e\rangle\),

$$
\rho=
\begin{pmatrix}
\rho_{gg}&\rho_{ge}\\
\rho_{eg}&\rho_{ee}
\end{pmatrix}.
$$

The diagonal elements are populations. The off-diagonal elements are coherences.

A physical density matrix satisfies

$$
\rho=\rho^\dagger,
$$

$$
\mathrm{Tr}\rho=1,
$$

and

$$
\rho\ge0.
$$

## 10.2 The actual 24-state matrix

For the basis

$$
\{|g_1\rangle,\ldots,|g_8\rangle,
|e_1\rangle,\ldots,|e_{16}\rangle\},
$$

the density matrix is exactly

$$
\rho=
\begin{pmatrix}
\rho_{gg}&\rho_{ge}\\
\rho_{eg}&\rho_{ee}
\end{pmatrix},
$$

with dimensions

$$
\rho_{gg}:8\times8,
\quad
\rho_{ge}:8\times16,
$$

$$
\rho_{eg}:16\times8,
\quad
\rho_{ee}:16\times16.
$$

![Block structure of the complete 24-state density matrix](figures/density_matrix_24x24.svg)

Written by indices,

$$
\rho_{ij}=\langle i|\rho|j\rangle,
\qquad i,j=1,\ldots,24.
$$

Examples:

$$
\rho_{44}
=P(F=2,m_F=-2),
$$

$$
\rho_{8,24}
=\langle F=2,m_F=+2|\rho|F'=3,m_F'=+3\rangle,
$$

which is the optical coherence of the stretched cycling transition.

The matrix contains

$$
24^2=576
$$

complex storage locations in the numerical representation. Hermiticity means only one triangular half is independent. A Hermitian \(24\times24\) matrix has \(24^2=576\) independent real parameters; the trace-one constraint reduces the physical parameter count to

$$
576-1=575
$$

real degrees of freedom before the nonlinear positivity restriction is considered.

Printing all 576 symbols individually would hide the physics. The block form above is the **complete matrix definition**: every one of the 576 entries is generated by \(\rho_{ij}\) with the explicit state ordering listed in Sec. 4.

---

# 11. The two-level optical Bloch equation

The rotating-frame Hamiltonian used for validation is

$$
\frac{H}{\hbar}
=
\begin{pmatrix}
0&\Omega^*/2\\
\Omega/2&-\delta
\end{pmatrix}.
$$

Here \(\Omega\) is the complex Rabi frequency and \(\delta\) is detuning.

Closed-system dynamics would give

$$
\dot\rho=-\frac{i}{\hbar}[H,\rho].
$$

Spontaneous decay is represented by

$$
C=\sqrt\Gamma|g\rangle\langle e|.
$$

Define the Lindblad dissipator

$$
\mathcal D[C]\rho
=C\rho C^\dagger
-\frac12\left(C^\dagger C\rho+\rho C^\dagger C\right).
$$

Then

$$
\boxed{
\dot\rho
=-\frac{i}{\hbar}[H,\rho]
+\mathcal D[C]\rho
}.
$$

With no laser drive,

$$
\rho_{ee}(t)=e^{-\Gamma t}.
$$

![Lindblad spontaneous decay](figures/lindblad_spontaneous_decay.svg)

With coherent drive, damped Rabi oscillations appear:

![Damped Rabi oscillations](figures/obe_rabi_oscillations.svg)

Using

$$
s=\frac{2|\Omega|^2}{\Gamma^2},
$$

the zero-pure-dephasing stationary population is

$$
\rho_{ee}
=\frac{s/2}{1+s+(2\delta/\Gamma)^2}.
$$

![OBE steady-state lineshape](figures/obe_steady_state_lorentzian.svg)

A travelling-wave force becomes

$$
\mathbf F
=\hbar\mathbf k\,\Gamma\rho_{ee}.
$$

The repository agrees with QuTiP at numerical precision in the matched two-level tests and reproduces the corresponding normalized PyLCP force.

---

# 12. Build the full 24-state Hamiltonian

It is convenient to work with

$$
h\equiv\frac{H}{\hbar},
$$

which has units rad/s.

In the ground/excited block basis,

$$
h(t)=
\begin{pmatrix}
h_g(t)&V^\dagger(t)\\
V(t)&h_e(t)
\end{pmatrix}.
$$

The blocks are:

- \(h_g\): \(8\times8\) ground hyperfine/Zeeman block plus rotating-frame carrier shifts;
- \(h_e\): \(16\times16\) excited hyperfine/Zeeman block;
- \(V\): \(16\times8\) optical coupling matrix.

## 12.1 Bare atomic blocks

The code evaluates

$$
h_g^{\rm bare}=\frac{H_{\rm hfs,g}+H_{Z,g}}{\hbar}-\omega_{g,\rm ref}\mathbf 1,
$$

$$
h_e^{\rm bare}=\frac{H_{\rm hfs,e}+H_{Z,e}}{\hbar}-\omega_{e,\rm ref}\mathbf 1.
$$

A global scalar energy origin is physically irrelevant, so subtracting reference energies changes no observable.

## 12.2 Optical coupling element

For beam \(i\) and allowed transition \(g\rightarrow e\), define

$$
\Omega_{i,ge}(\mathbf r)
=\Gamma
\sqrt{\frac{s_i(\mathbf r)S_{ge}P_i(q)}{2}}.
$$

The code inserts

$$
V_{eg}^{(i)}(t)
=\frac{\Omega_{i,ge}}{2}
e^{i\phi_i(t)}.
$$

The Hermitian-conjugate element is

$$
V_{ge}^{(i)}(t)
=\left[V_{eg}^{(i)}(t)\right]^*.
$$

For a moving atom,

$$
\mathbf r(t)=\mathbf r_0+\mathbf vt,
$$

and the phase is

$$
\phi_i(t)
=\mathbf k_i\cdot[\mathbf r(t)-\mathbf r_{i,0}]
-\Delta\omega_i t
+\phi_{i,0}.
$$

Thus the Doppler shift arises directly from the time derivative of the optical phase.

---

# 13. Why the code uses a block rotating frame

The cooling and repump frequencies differ by approximately the ground hyperfine splitting, several GHz. Explicit integration of that carrier beat would force a time step of order tenths of a nanosecond even though optical pumping and mechanical motion are much slower.

For each beam family, the reference laser angular-frequency offset is

$$
\omega_{L,i}^{\rm off}
=2\pi
\left[
(\nu^e_{F'_t}-\nu^e_{F'_{\max}})
-(\nu^g_F-\nu^g_{F_{\max}})
\right]
+\delta_i+\omega_{{\rm AOM},i}.
$$

For each addressed ground \(F\) block, one carrier \(\omega_F^{\rm car}\) is chosen. The retained beat for beam \(i\) is then

$$
\Delta\omega_i^{\rm retained}
=\omega_{L,i}^{\rm off}
-\omega_F^{\rm car}
-\mathbf k_i\cdot\mathbf v.
$$

The several-GHz cooling-repump separation is therefore removed analytically; same-manifold AOM and Doppler differences remain explicit.

## Controlled cross-ground RWA

A laser assigned to one ground hyperfine manifold does not explicitly drive the other ground manifold. The code estimates the largest omitted amplitude by

$$
\epsilon_{\rm RWA}
=\frac{\Omega_{\max}}{\Delta_{\min}},
$$

with a population-scale estimate

$$
P_{\rm omitted}\sim\epsilon_{\rm RWA}^2.
$$

For the reference MOT,

$$
\epsilon_{\rm RWA}\approx1.65\times10^{-3},
$$

$$
P_{\rm omitted}\lesssim2.72\times10^{-6}.
$$

Cross-\(F\) magnetic elements are also secularly removed in the independently rotating ground blocks. This is a weak-field approximation and would need reconsideration near the hyperfine Paschen-Back regime.

---

# 14. Spontaneous emission in the 24-state model

For every allowed spontaneous decay \(e\rightarrow g\), the code creates one collapse operator

$$
C_{ge}
=\sqrt{\Gamma b_{e\rightarrow g}}
|g\rangle\langle e|.
$$

The complete master equation is

$$
\boxed{
\dot\rho
=-i[h(t),\rho]
+\sum_{g,e}\mathcal D[C_{ge}]\rho
}.
$$

If optional explicit ground-state mixing is requested, extra collapse operators are added:

$$
C_{ts}^{\rm mix}
=\sqrt{\frac{\gamma_{\rm mix}}{N_g}}
|t\rangle\langle s|,
$$

where \(N_g=8\). The default mixing rate is zero; it is not silently introduced as a numerical fix.

---

# 15. Convert the \(24\times24\) equation into a 576-dimensional linear system

To solve the master equation numerically, the matrix is vectorized column by column:

$$
|\rho\rangle\rangle
\equiv\operatorname{vec}(\rho).
$$

Since \(24\times24=576\),

$$
|\rho\rangle\rangle
\in\mathbb C^{576}.
$$

Use the identity

$$
\operatorname{vec}(A\rho B)
=(B^T\otimes A)\operatorname{vec}(\rho).
$$

The coherent term becomes

$$
\mathcal L_H
=-i\left(\mathbf1\otimes h-h^T\otimes\mathbf1\right).
$$

For one collapse operator,

$$
\mathcal L_C
=C^*\otimes C
-\frac12\mathbf1\otimes C^\dagger C
-\frac12(C^\dagger C)^T\otimes\mathbf1.
$$

Therefore the full sparse Liouvillian is

$$
\boxed{
\mathcal L
=-i(\mathbf1\otimes h-h^T\otimes\mathbf1)
+\sum_c\mathcal L_{C_c}
}.
$$

and

$$
\boxed{
\frac{d}{dt}|\rho\rangle\rangle
=\mathcal L(t)|\rho\rangle\rangle.
}
$$

![From 24x24 density matrix to 576x576 Liouvillian](figures/liouvillian_576.svg)

The matrix \(\mathcal L\) is formally \(576\times576\), but is stored sparsely because most atomic states are not directly connected by one dipole transition.

## 15.1 Stationary solution

For a time-independent rotating frame,

$$
\mathcal L|\rho_{\rm ss}\rangle\rangle=0.
$$

The null equation alone has an arbitrary normalization. One row is therefore replaced with

$$
\mathrm{Tr}\rho=1.
$$

The resulting sparse linear system is solved directly.

## 15.2 Time-dependent solution

When residual optical beats or the magnetic field are time dependent,

$$
\dot{\rho}=\mathcal L(t)\rho
$$

is integrated with adaptive `solve_ivp`. The code determines time dependence from known physical beat frequencies and from the magnetic-field object's declared time dependence; it does not infer stationarity from two coincident time samples.

---

# 16. Force from the full quantum Hamiltonian

The force is not guessed from total fluorescence. For beam \(i\),

$$
\boxed{
\mathbf F_i
=-\left\langle\nabla H_i\right\rangle
=-\hbar\operatorname{Tr}
\left[
\rho\nabla\left(\frac{H_i}{\hbar}\right)
\right].
}
$$

For a collimated Gaussian travelling wave, the complex field amplitude contains

$$
\sqrt I\,e^{i\mathbf k\cdot\mathbf r}.
$$

Its gradient gives

$$
\nabla\ln E
=i\mathbf k
-2\frac{\mathbf r_\perp}{w^2}.
$$

The first term produces travelling-wave radiation pressure; the second gives the Gaussian-envelope dipole-force contribution.

The code independently checked this analytical gradient against finite differences. In the travelling-wave limit it reproduces

$$
F=\hbar k\Gamma\rho_{ee}
$$

to a reported relative error of approximately

$$
2.43\times10^{-15}.
$$

---

# 17. Coherent and incoherent laser groups

If two beams are phase coherent, their relative phase is a physical control parameter and must be retained.

If two groups are mutually incoherent, an arbitrary fixed relative phase is not physical. The code therefore computes

$$
\langle O\rangle_{\phi}
=\frac{1}{N_\phi}
\sum_{n=1}^{N_\phi}O(\phi_n)
$$

for the requested observable \(O\), then refines \(N_\phi\) until the result converges.

A useful convergence metric is

$$
\epsilon_\phi
=\frac{\|O_{2N}-O_N\|}
{\max(\|O_{2N}\|,\|O_N\|,O_{\rm floor})}.
$$

Research mode requires the phase-averaged observable to meet the configured tolerance.

---

# 18. Polarization-gradient cooling as a controlled reduced model

The full OBE is expensive, so the repository also contains a transparent adiabatic population model for the closed

$$
F=2\rightarrow F'=3
$$

manifold.

The coherent field is

$$
\mathbf E(\mathbf r)
\propto\sum_b
\sqrt{s_b}\,\boldsymbol\epsilon_b
\exp[i(\mathbf k_b\cdot\mathbf r+\phi_b)].
$$

Projection onto spherical components gives local \(s_q(\mathbf r)\).

For ground state \(m\) and excited state \(m'\), the transition detuning including the projected linear Zeeman shift is

$$
\delta_{mm'}
=\Delta
-\frac{\mu_BB_\parallel}{\hbar}
(g_e m'-g_gm).
$$

The adiabatically eliminated light shift is

$$
U_m(\mathbf r)
=\sum_{m',q}
\frac{\hbar\delta_{mm'}\Gamma^2
S_{mm'q}s_q(\mathbf r)}
{8[\delta_{mm'}^2+(\Gamma/2)^2]}.
$$

The excitation/pumping rate is

$$
R_{mm'}(\mathbf r)
=\frac{\Gamma^3S_{mm'q}s_q(\mathbf r)}
{8[\delta_{mm'}^2+(\Gamma/2)^2]}.
$$

The five ground populations obey

$$
\dot{\mathbf p}=A_{\rm pump}(\mathbf r)\mathbf p.
$$

The conservative Sisyphus force is

$$
\mathbf F_{\rm Sis}
=-\sum_mp_m\nabla U_m.
$$

For motion through the lattice, the code integrates \(\mathbf p(t)\), discards the transient optical periods, and averages the force.

![PGC light shifts and pumping](../../results/polarization_gradient/light_shifts_pumping.png)

![Sub-Doppler force versus velocity](../../results/polarization_gradient/subdoppler_force_velocity.png)

A local friction coefficient is

$$
\beta_v
=-\left.\frac{\partial F}{\partial v}\right|_{v=0}.
$$

The present recoil-only diffusion tensor is

$$
D_{pp}^{\rm recoil}
=\frac{(\hbar k)^2R_{\rm sc}}{2}
\left[
\sum_b w_b\hat{\mathbf k}_b\hat{\mathbf k}_b^T
+\frac13\mathbf1
\right].
$$

It includes absorption shot noise and isotropic spontaneous recoil but omits internal-state and dipole-force fluctuations.

Therefore the repository deliberately does **not** combine the full coherent force with this reduced diffusion tensor to claim a quantitative temperature.

---

# 19. Stray magnetic fields and PGC

The \(F=2\) Larmor frequency scale is approximately

$$
\frac{\omega_L}{2\pi}
\approx699.6\ \mathrm{Hz/mG}.
$$

For the reference \(\Delta=-3\Gamma,\ s=0.08\) recipe, a simple weak-drive optical-pumping scale is about

$$
R_{\rm pump}\approx6.56\ \mathrm{kHz}.
$$

Equating these simple rates gives

$$
B\sim9.4\ \mathrm{mG}.
$$

This is a **timescale marker**, not a predicted temperature threshold.

A defensible temperature curve would require, at matching fidelity,

$$
\beta_v(B)
=-\left.\frac{\partial F(v,B)}{\partial v}\right|_0
$$

and a full momentum-diffusion tensor \(D_{pp}(B)\). Only in the appropriate linear Fokker-Planck regime could one then use an Einstein-type relation such as

$$
k_BT\sim\frac{D_{pp}}{\beta_v}
$$

with carefully matched definitions.

![Vector residual-field OBE diagnostic](../../results/polarization_gradient/vector_residual_obe.png)

---

# 20. External atomic motion

Once an internal-state model supplies a mean force,

$$
\frac{d\mathbf r}{dt}=\mathbf v,
$$

$$
m\frac{d\mathbf v}{dt}=\mathbf F(\mathbf r,\mathbf v,t).
$$

The deterministic code uses adaptive RK45.

Absorption from beam \(i\) transfers

$$
\Delta\mathbf p_{\rm abs}=+\hbar\mathbf k_i.
$$

Spontaneous emission gives a recoil of magnitude \(\hbar k\) in a sampled direction. For isotropic recoil,

$$
\langle\Delta\mathbf p_{\rm em}\rangle=0,
$$

and for one Cartesian component

$$
\langle(\Delta p_x)^2\rangle
=\frac{(\hbar k)^2}{3}.
$$

![Deterministic MOT trajectories](../../results/effective_mot/deterministic_trajectories.png)

---

# 21. Time-dependent experimental sequence

A laboratory sequence is represented as successive control stages:

$$
\text{MOT load}
\rightarrow\text{CMOT}
\rightarrow\text{field off}
\rightarrow\text{settling}
\rightarrow\text{PGC/molasses}
\rightarrow\text{TOF}.
$$

A smoothstep ramp uses

$$
f(u)=3u^2-2u^3,
\qquad0\le u\le1,
$$

so a parameter \(X\) becomes

$$
X(u)=X_0+[X_1-X_0]f(u).
$$

An example current/gradient decay is

$$
G(t)=G_0e^{-(t-t_0)/\tau_{\rm coil}}.
$$

Residual magnetic field can include

$$
\mathbf B(t)
=\mathbf B_{\rm DC}
+\mathbf B_{\rm eddy}e^{-(t-t_0)/\tau_{\rm eddy}}
+\mathbf B_{\rm AC}\sin(2\pi ft+\phi).
$$

![Experimental sequence](../../results/sequence/sequence_timeline.png)

---

# 22. Vapour loading: from thermal gas to captured atoms

The ideal-gas number density is

$$
n=\frac{P}{k_BT}.
$$

The one-sided thermal flux through a surface is

$$
\frac{\Phi}{A}
=n\sqrt{\frac{k_BT}{2\pi m}}
=\frac{n\langle v\rangle}{4}.
$$

Atoms crossing a surface do **not** follow the ordinary Maxwell speed distribution. Flux weighting gives

$$
p_{\rm flux}(v)
\propto
v^3\exp\left(-\frac{mv^2}{2k_BT}\right).
$$

![Bulk Maxwell versus surface-flux speed distribution](figures/thermal_flux_distribution.svg)

The trajectory-derived loading rate is

$$
R_{\rm load}
=\Phi_{\rm incident}P_{\rm capture}.
$$

![Capture response map](../../results/capture_loading/capture_response_map.png)

---

# 23. Loading and losses

With one-body loss,

$$
\dot N=R_{\rm load}-\gamma N.
$$

The solution is

$$
N(t)=N_\infty+[N(0)-N_\infty]e^{-\gamma t},
$$

with

$$
N_\infty=\frac{R_{\rm load}}{\gamma}.
$$

With two-body loss,

$$
\dot N
=R_{\rm load}-\gamma N
-\beta_2\int n^2(\mathbf r)d^3r.
$$

For a Gaussian cloud,

$$
n(\mathbf r)
=\frac{N}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z}
\exp\left[-\frac12\left(
\frac{x^2}{\sigma_x^2}
+\frac{y^2}{\sigma_y^2}
+\frac{z^2}{\sigma_z^2}
\right)\right].
$$

Then

$$
\int n^2d^3r
=\frac{N^2}{8\pi^{3/2}\sigma_x\sigma_y\sigma_z}.
$$

Define

$$
V_{2,\rm eff}
=8\pi^{3/2}\sigma_x\sigma_y\sigma_z.
$$

Then

$$
\boxed{
\dot N
=R_{\rm load}-\gamma N
-\frac{\beta_2}{V_{2,\rm eff}}N^2.
}
$$

![Loading/loss dynamics](figures/loading_loss_dynamics.svg)

The code does not invent \(\beta_2\) or collision cross sections. They must be literature-, user-, or experiment-supplied quantities.

---

# 24. Optional collective MOT model

For a Gaussian cloud, the peak density is

$$
n_0
=\frac{N}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z}.
$$

The central column density along axis \(i\) is

$$
\mathcal N_i
=\frac{N}{2\pi\sigma_j\sigma_k}.
$$

Optical depth is

$$
\mathrm{OD}_i=\sigma_{\rm opt}\mathcal N_i.
$$

A simple reabsorption probability proxy is

$$
P_{\rm reabs}
=1-e^{-\langle\mathrm{OD}\rangle}.
$$

The mean-field multiple-scattering coefficient is

$$
Q
=\frac{\sigma_L\sigma_R I_{\rm tot}}{4\pi c}.
$$

The corresponding radial repulsive force is modeled as

$$
F_{\rm rep}(r)
=Q\frac{N_{\rm enc}(r)}{r^2}\operatorname{sgn}(r).
$$

For a spherical Gaussian cloud,

$$
\frac{N_{\rm enc}(r)}{N}
=
\operatorname{erf}\left(\frac{r}{\sqrt2\sigma}\right)
-\sqrt{\frac2\pi}\frac{r}{\sigma}
\exp\left(-\frac{r^2}{2\sigma^2}\right).
$$

This is an optional mean-field approximation, not a full radiative-transfer calculation.

![Collective MOT diagnostics](../../results/collective_mot/collective_mot_diagnostics.png)

---

# 25. The calculation as an algorithm

A complete reference calculation can now be written as a sequence.

## Step 1 — load atomic constants

Read \(m,\lambda,\tau,A_{\rm hfs},B_{\rm hfs},g_J,g_I,I,I_{\rm sat}\).

## Step 2 — generate hyperfine levels

For every allowed \(F\), calculate \(K\) and \(E_{\rm hfs}(F)\).

## Step 3 — build all 24 states

Create the 8 ground and 16 excited \(|F,m_F\rangle\) states in the fixed ordering of Sec. 4.

## Step 4 — generate the dipole graph

For every ground/excited pair, test \(\Delta F\) and \(\Delta m_F\), then calculate the Wigner-6j and Clebsch-Gordan strength.

## Step 5 — calculate spontaneous branching

Normalize the allowed decay strengths from each excited state.

## Step 6 — construct six beam objects

For each beam calculate \(I(\mathbf r)\), \(s(\mathbf r)\), \(\mathbf k\), \(\boldsymbol\epsilon\), spherical \(P_q\), phase, linewidth, and coherence group.

## Step 7 — construct the magnetic field

Use either the ideal quadrupole field or physical coils plus measured/configured residual fields.

## Step 8 — choose fidelity

- effective force for broad trajectory/capture scans;
- 24-population rate equation when optical pumping matters but coherence does not;
- 24-state OBE when coherences/vector Zeeman physics matter;
- reduced PGC model for transparent Sisyphus mechanism studies.

## Step 9 — if using the 24-state OBE, build \(h(t)\)

Construct \(h_g,h_e,V\), apply rotating-frame carrier shifts, and assemble

$$
h(t)=
\begin{pmatrix}
h_g&V^\dagger\\V&h_e\end{pmatrix}.
$$

## Step 10 — construct all collapse operators

For every allowed spontaneous branch, add

$$
C_{ge}=\sqrt{\Gamma b_{e\to g}}|g\rangle\langle e|.
$$

## Step 11 — build the Liouvillian

Construct the sparse \(576\times576\) matrix \(\mathcal L\).

## Step 12 — solve the internal state

Stationary case:

$$
\mathcal L|\rho_{\rm ss}\rangle\rangle=0,
\quad \mathrm{Tr}\rho=1.
$$

Time-dependent case:

$$
\frac{d}{dt}|\rho\rangle\rangle
=\mathcal L(t)|\rho\rangle\rangle.
$$

## Step 13 — calculate force

$$
\mathbf F_i
=-\hbar\mathrm{Tr}\left[
\rho\nabla(H_i/\hbar)
\right].
$$

## Step 14 — propagate motion

$$
\dot{\mathbf r}=\mathbf v,
\qquad
m\dot{\mathbf v}=\mathbf F.
$$

## Step 15 — if studying loading, sample the correct incident thermal flux

Use \(p_{\rm flux}(v)\propto v^3e^{-mv^2/(2k_BT)}\), integrate trajectories, apply the explicit capture criterion, and obtain \(P_{\rm capture}\).

## Step 16 — calculate loading

$$
R_{\rm load}=\Phi P_{\rm capture},
$$

then integrate the chosen loss equation.

## Step 17 — validate the level appropriate to the claim

The repository currently has independent validation for the two-level OBE, normalized two-beam force, and \(^{87}\mathrm{Rb}\) ground vector-Zeeman spectrum. Full 24-state \(^{87}\mathrm{Rb}\) force/population validation against a matched PyLCP calculation remains an important pending gate.

---

# 26. Why several models remain instead of one enormous calculation

The 24-state OBE is the most complete internal-state model currently implemented, but using it for every thermal atom in a large three-dimensional capture Monte Carlo would be prohibitively expensive.

The hierarchy is therefore deliberate:

$$
\text{effective force}
\subset
\text{rate equations}
\subset
\text{coherent OBE}
$$

in the sense of increasing retained internal physics, not strict mathematical set inclusion.

The scientific rule is:

> **Use the least expensive model that still contains the physics required by the observable being claimed.**

This is why an effective force is acceptable for broad capture studies but not for a quantitative sub-Doppler temperature, and why the code refuses to manufacture a precise \(T(B)\) from a high-fidelity force and a lower-fidelity diffusion model.

---

# 27. What is fully established and what is still pending

Externally checked at present:

- two-level analytical and QuTiP steady-state OBE;
- two-level Liouvillian dynamics;
- normalized two-beam Doppler force versus PyLCP;
- \(^{87}\mathrm{Rb}\) vector-Zeeman spectrum versus PyLCP.

Implemented but not yet fully externally matched:

- complete 24-state moving-\(^{87}\mathrm{Rb}\) OBE force/populations;
- quantitatively matched PGC force/diffusion and final temperature;
- apparatus-specific atom number and temperature without measured calibration inputs.

That distinction is intentional. The textbook shows the whole implemented calculation while preserving the boundary between **implemented mathematics** and **experimentally validated prediction**.
