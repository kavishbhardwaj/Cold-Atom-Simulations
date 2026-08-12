# Before the tutorial — notation, symbols, and complete equation map

This page exists to prevent the most common problem in technical tutorials: a symbol suddenly appearing in an equation without the reader knowing what it means.

The rule used throughout the tutorial is:

> **Every symbol is defined before or immediately after it first appears.**

A second rule is equally important:

> **The tutorial includes every governing physical equation that directly determines a core result in this repository.**

It does not list trivial array reshaping, plotting transforms, or every line of numerical bookkeeping as a separate “equation.” Numerical definitions that materially affect a scientific result — for example the capture criterion, Wilson confidence interval, phase averaging, or ODE tolerances — are explained where they enter the calculation.

## 0.1 Frequency and unit convention

The code uses SI units internally. Laser detunings and frequency offsets used inside Hamiltonians are angular frequencies, so they are measured in rad/s. The spontaneous decay rate

$$
\Gamma=\frac{1}{\tau}
$$

has units s$^{-1}$ and is compared directly with those angular-frequency detunings. For spectroscopy it is often easier to quote

$$
\frac{\Gamma}{2\pi}
$$

in MHz. Whenever the tutorial writes a detuning such as $\Delta=-3\Gamma$, it means the angular detuning is minus three natural linewidths in the convention used by the code.

## 0.2 Core symbols used repeatedly

| Symbol | Meaning |
|---|---|
| $I$ | nuclear spin quantum number; for $^{87}$Rb, $I=3/2$ |
| $J$ | total electronic angular momentum in one fine-structure manifold |
| $F$ | total hyperfine angular momentum, $\mathbf F=\mathbf I+\mathbf J$ |
| $m_F$ | projection of $F$ on the chosen quantization axis |
| $A_{\rm hfs}$ | magnetic-dipole hyperfine constant, stored in Hz |
| $B_{\rm hfs}$ | electric-quadrupole hyperfine constant, stored in Hz; **not** the magnetic field |
| $\mathbf B$ | magnetic-field vector in tesla |
| $g_J,g_I,g_F$ | electronic, nuclear and hyperfine Landé factors |
| $\mu_B$ | Bohr magneton |
| $\lambda$ | optical wavelength |
| $k=2\pi/\lambda$ | optical wave number |
| $\tau$ | excited-state lifetime |
| $\Gamma=1/\tau$ | spontaneous population-decay rate |
| $P$ | laser power |
| $I_{\rm opt}$ | optical intensity; the subscript avoids confusion with nuclear spin $I$ |
| $I_{\rm sat}$ | saturation intensity |
| $s=I_{\rm opt}/I_{\rm sat}$ | saturation parameter |
| $\delta,\Delta$ | laser detuning from resonance |
| $\Omega$ | complex Rabi frequency |
| $\rho$ | density matrix |
| $p_i$ | population of state $i$ in a rate-equation model |
| $\mathbf r,\mathbf v$ | atomic position and velocity |
| $\mathbf F$ | mechanical force |
| $\kappa$ | small-displacement MOT spring constant |
| $\beta_v$ | velocity-damping coefficient; the subscript distinguishes it from two-body loss $\beta_2$ |
| $D_{pp}$ | momentum-diffusion coefficient/tensor |
| $R_{\rm load}$ | MOT loading rate in atoms/s |
| $\gamma$ | one-body loss rate |
| $\beta_2$ | two-body loss coefficient in m$^3$/s |
| $N$ | trapped atom number |
| $n(\mathbf r)$ | atomic number density |

---

# 0.3 Atomic-structure equations

## Hyperfine energy

The hyperfine interaction couples the nuclear angular momentum $\mathbf I$ to the electronic angular momentum $\mathbf J$. In a coupled state $|F,m_F\rangle$, define

$$
K=F(F+1)-I(I+1)-J(J+1).
$$

The symbol $A_{\rm hfs}$ is the **magnetic-dipole hyperfine constant**. It sets the energy scale of the $\mathbf I\cdot\mathbf J$ interaction. The magnetic-dipole contribution is

$$
\frac{E_A}{h}=\frac{A_{\rm hfs}}{2}K.
$$

The symbol $B_{\rm hfs}$ is the **electric-quadrupole hyperfine constant**. It exists only when both $I\ge1$ and $J\ge1$. The repository uses

$$
\frac{E_{\rm hfs}}{h}
=
\frac{A_{\rm hfs}}{2}K
+
B_{\rm hfs}
\frac{\frac34K(K+1)-I(I+1)J(J+1)}
{2I(2I-1)J(2J-1)}.
$$

For the $5S_{1/2}$ ground state, $J=1/2$, so the quadrupole term is absent. For the $5P_{3/2}$ excited state, it is retained.

## Hyperfine Landé factor

For weak magnetic fields the coupled-state Landé factor is

$$
g_F=
\frac{
 g_J\left[F(F+1)+J(J+1)-I(I+1)\right]
+g_I\left[F(F+1)+I(I+1)-J(J+1)\right]
}{2F(F+1)}.
$$

For $F=0$, the code returns $g_F=0$.

## Hyperfine and Zeeman-resolved dipole strengths

The relative hyperfine reduced strength used by the code is

$$
S_{F\rightarrow F'}
\propto
(2F'+1)(2J_g+1)
\begin{Bmatrix}
J_e & F' & I\\
F & J_g & 1
\end{Bmatrix}^{2}.
$$

The Zeeman-resolved strength then contains the Clebsch-Gordan factor

$$
S_{F m\rightarrow F' m'}
\propto
S_{F\rightarrow F'}
\left|
\langle F,m;1,q|F',m'\rangle
\right|^2,
$$

where

$$
q=m'-m\in\{-1,0,+1\}
$$

labels $\sigma^-$, $\pi$, and $\sigma^+$ light. The code normalizes the strongest generated transition to unit relative strength and obtains spontaneous branching by normalizing all allowed decay strengths from each excited state.

---

# 0.4 Magnetic-field equations

## Ideal MOT quadrupole

The fast field model is

$$
\mathbf B(\mathbf r)=
R\,\mathrm{diag}(b',b',-2b')R^T(\mathbf r-\mathbf r_0),
$$

where $b'$ is the radial gradient, $R$ is an optional rotation matrix and $\mathbf r_0$ is the field-zero position. Its gradient is traceless,

$$
\nabla\cdot\mathbf B=0.
$$

## Residual and AC field

A configurable background can contain

$$
\mathbf B_{\rm res}(\mathbf r,t)
=G_{\rm stray}\mathbf r+\mathbf B_0
+\mathbf B_{\rm ac}\sin(2\pi f t+\phi).
$$

## Biot-Savart coil field

Physical circular coils are calculated from

$$
d\mathbf B=
\frac{\mu_0 I_c N_t}{4\pi}
\frac{d\boldsymbol\ell\times(\mathbf r-\mathbf r')}{|\mathbf r-\mathbf r'|^3},
$$

where $I_c$ is coil current and $N_t$ is turn count.

## Full hyperfine-Zeeman Hamiltonian

The coherent model uses

$$
H=H_{\rm hfs}+H_Z,
$$

with

$$
H_Z=\mu_B\left(g_J\mathbf J+g_I\mathbf I\right)\cdot\mathbf B.
$$

Within one fine-structure manifold the hyperfine Hamiltonian is equivalently implemented as

$$
\frac{H_{\rm hfs}}{h}
=A_{\rm hfs}\,\mathbf I\cdot\mathbf J
+B_{\rm hfs}
\frac{3(\mathbf I\cdot\mathbf J)^2+\frac32\mathbf I\cdot\mathbf J-I(I+1)J(J+1)}
{2I(2I-1)J(2J-1)}.
$$

At weak field this approaches

$$
\Delta E\approx g_F\mu_Bm_FB.
$$

---

# 0.5 Laser-beam equations

The optical wave vector is

$$
\mathbf k=\frac{2\pi}{\lambda}\hat{\mathbf k}.
$$

For an elliptical Gaussian beam the repository uses

$$
I_{\rm opt}(x,y,z)
=
\frac{2P}{\pi w_x(z)w_y(z)}
\exp\left[-2\left(\frac{x^2}{w_x^2(z)}+\frac{y^2}{w_y^2(z)}\right)\right].
$$

The Rayleigh ranges and propagated waists are

$$
z_{R,x}=\frac{\pi w_{0x}^2}{\lambda},
\qquad
z_{R,y}=\frac{\pi w_{0y}^2}{\lambda},
$$

$$
w_x(z)=w_{0x}\sqrt{1+(z/z_{R,x})^2},
\qquad
w_y(z)=w_{0y}\sqrt{1+(z/z_{R,y})^2}.
$$

The saturation parameter is

$$
s(\mathbf r)=\frac{I_{\rm opt}(\mathbf r)}{I_{\rm sat}}.
$$

For mutually coherent beams the electric fields add before squaring. For mutually incoherent groups, intensities are added after coherent summation within each group:

$$
I_{\rm total}
=
\sum_g\left|\sum_{i\in g}\mathbf E_i\right|^2
+\sum_{j\in\mathrm{independent}}I_j.
$$

---

# 0.6 Effective MOT force equations

The polarization spin vector used in the fast Zeeman approximation is

$$
\mathbf s_{\rm pol}=\mathrm{Re}\left(i\,\boldsymbol\epsilon\times\boldsymbol\epsilon^*\right),
$$

and its projection on the local magnetic axis produces an effective photon angular-momentum expectation $\bar q$. The corresponding effective Zeeman detuning is

$$
\delta_Z=-\frac{\mu_{\rm eff}}{\hbar}\bar q\,|\mathbf B|.
$$

For beam $i$,

$$
\delta_i=\delta_{L,i}+\delta_{{\rm AOM},i}-\mathbf k_i\cdot\mathbf v+\delta_{Z,i}.
$$

If the laser has angular FWHM linewidth $\gamma_{L,i}$, the width used by the code is

$$
\Gamma_i^{\rm eff}=\Gamma+\gamma_{L,i}.
$$

The effective scattering rate is

$$
R_i=\frac{\Gamma}{2}
\frac{s_i\,\Gamma/\Gamma_i^{\rm eff}}
{1+\sum_j s_j+\left(2\delta_i/\Gamma_i^{\rm eff}\right)^2}.
$$

The mean force is

$$
\mathbf F_{\rm opt}=\sum_i\hbar\mathbf k_iR_i,
\qquad
\mathbf F=\mathbf F_{\rm opt}+m\mathbf g.
$$

Near the trap centre,

$$
F_x\approx-\kappa x-\beta_v v_x,
$$

with

$$
\kappa=-\left.\frac{\partial F_x}{\partial x}\right|_0,
\qquad
\beta_v=-\left.\frac{\partial F_x}{\partial v_x}\right|_0.
$$

---

# 0.7 Multilevel population-rate equations

For a beam family $b$ and transition $g\leftrightarrow e$, define

$$
s_{b,ge}^{\rm eff}=s_bS_{ge}P_b(q),
$$

where $S_{ge}$ is the generated transition strength and $P_b(q)$ is the local spherical-polarization fraction. The code uses

$$
W_{b,ge}=\frac{\Gamma}{2}
\frac{s_{b,ge}^{\rm eff}\,\Gamma/\Gamma_b^{\rm eff}}
{1+\left(2\delta_{b,ge}/\Gamma_b^{\rm eff}\right)^2}.
$$

Stimulated absorption and stimulated emission enter symmetrically. The complete population vector obeys

$$
\dot{\mathbf p}=A_{\rm rate}\mathbf p,
$$

where $A_{\rm rate}$ is the column-conservative rate generator assembled from all stimulated and spontaneous processes. Written state by state,

$$
\dot p_e=\sum_{g,b}W_{b,ge}(p_g-p_e)-\Gamma p_e,
$$

$$
\dot p_g=\sum_{e,b}W_{b,ge}(p_e-p_g)+\Gamma\sum_e b_{e\rightarrow g}p_e.
$$

The per-beam radiation-pressure force is

$$
\mathbf F_b=\hbar\mathbf k_b\sum_{ge}W_{b,ge}(p_g-p_e).
$$

---

# 0.8 Two-level OBE and Lindblad equations

In the rotating frame the two-level Hamiltonian is

$$
\frac{H}{\hbar}=\begin{pmatrix}0&\Omega^*/2\\\Omega/2&-\delta\end{pmatrix}.
$$

$\Omega$ is the complex Rabi frequency: it measures the coherent laser coupling between $|g\rangle$ and $|e\rangle$. The spontaneous-emission collapse operator is

$$
C=\sqrt{\Gamma}\,|g\rangle\langle e|.
$$

For any collapse operator $C$, define the Lindblad dissipator

$$
\mathcal D[C]\rho=C\rho C^\dagger-\frac12\left(C^\dagger C\rho+\rho C^\dagger C\right).
$$

The master equation is

$$
\dot\rho=-\frac{i}{\hbar}[H,\rho]+\mathcal D[C]\rho.
$$

Optional pure dephasing is represented by

$$
C_\phi=\sqrt{\gamma_\phi/2}\,\sigma_z.
$$

Define the transverse coherence-decay rate

$$
\gamma_\perp=\frac{\Gamma}{2}+\gamma_\phi.
$$

The code's analytical steady-state population is

$$
\rho_{ee}=\frac{|\Omega|^2\gamma_\perp}{2\Gamma(\delta^2+\gamma_\perp^2)+2|\Omega|^2\gamma_\perp}.
$$

For zero extra dephasing and

$$
s=\frac{2|\Omega|^2}{\Gamma^2},
$$

this reduces to

$$
\rho_{ee}=\frac{s/2}{1+s+(2\delta/\Gamma)^2}.
$$

A single travelling wave has mean scattering force

$$
\mathbf F=\hbar\mathbf k\,\Gamma\rho_{ee}.
$$

---

# 0.9 Full moving multilevel OBE

For the 24-state basis, the density matrix still obeys the same structural master equation,

$$
\dot\rho=-i[h(t),\rho]+\sum_c\mathcal D[C_c]\rho,
$$

where $h=H/\hbar$ is stored in angular-frequency units and the collapse operators $C_c$ are generated from the spontaneous-branching matrix.

For an atom moving with constant velocity during one internal-state solve,

$$
\mathbf r(t)=\mathbf r_0+\mathbf vt.
$$

Beam $i$ contributes the phase

$$
\phi_i(t)=\mathbf k_i\cdot(\mathbf r_0+\mathbf vt)-\delta\omega_i t+\phi_{i,0}.
$$

For transition $g\leftrightarrow e$, the code uses the coupling scale

$$
\Omega_{i,ge}=\Gamma\sqrt{\frac{s_iS_{ge}P_i(q)}{2}}.
$$

The force from beam $i$ is evaluated from the Hamiltonian gradient,

$$
\mathbf F_i=-\hbar\,\mathrm{Tr}\left[\rho\,\nabla\left(\frac{H_i}{\hbar}\right)\right].
$$

For the Gaussian travelling-wave amplitude the analytic gradient contains

$$
\nabla E_i\supset E_i\left(i\mathbf k_i-\frac{2\mathbf r_\perp}{w^2}\right).
$$

The discarded cross-ground optical coupling is monitored with

$$
\epsilon_{\rm RWA}=\frac{\Omega_{\max}}{\Delta_{\min}},
\qquad
P_{\rm discarded}\sim\epsilon_{\rm RWA}^2.
$$

---

# 0.10 Polarization-gradient / Sisyphus equations

A coherent beam contributes

$$
\mathbf E_b(\mathbf r)\propto\sqrt{s_b}\,\boldsymbol\epsilon_b\exp[i(\mathbf k_b\cdot\mathbf r+\phi_b)].
$$

For ground sublevel $m$ and polarization component $q$, the adiabatically eliminated light shift is

$$
U_m(\mathbf r)=\sum_q\frac{\hbar\,\delta_{mq}\Gamma^2C_{mq}^2s_q(\mathbf r)}{8[\delta_{mq}^2+(\Gamma/2)^2]}.
$$

The corresponding optical-pumping rate is

$$
R_{mq}(\mathbf r)=\frac{\Gamma^3C_{mq}^2s_q(\mathbf r)}{8[\delta_{mq}^2+(\Gamma/2)^2]}.
$$

The state-resolved conservative force is

$$
\mathbf F_m=-\nabla U_m,
$$

and the population-weighted force is

$$
\mathbf F=\sum_m p_m\mathbf F_m.
$$

For motion through the lattice, the ground populations obey

$$
\dot{\mathbf p}=A_{\rm pump}[\mathbf r(t)]\mathbf p.
$$

The low-velocity friction estimate is

$$
\beta_v=-\frac{F(+v)-F(-v)}{2v}.
$$

The reduced recoil-only diffusion tensor is

$$
\mathbf D_{pp}^{\rm recoil}=\frac{(\hbar k)^2\bar R_{\rm sc}}{2}\left[\sum_b w_b\hat{\mathbf k}_b\hat{\mathbf k}_b^T+\frac13\mathbf I_3\right].
$$

It deliberately omits internal-state and dipole-force fluctuations, so it is **not** combined with the coherent OBE friction to claim a quantitative PGC temperature.

---

# 0.11 External motion and recoil

The semiclassical trajectory equations are

$$
\dot{\mathbf r}=\mathbf v,
\qquad
m\dot{\mathbf v}=\mathbf F(\mathbf r,\mathbf v,t).
$$

A photon absorbed from beam $i$ gives

$$
\Delta\mathbf p_{\rm abs}=+\hbar\mathbf k_i.
$$

Isotropic spontaneous emission has zero mean recoil and component variance

$$
\left\langle\Delta p_j^2\right\rangle=\frac{(\hbar k)^2}{3}.
$$

The Bernoulli photon-event time step is constrained by

$$
R_{\rm tot}\Delta t\lesssim0.1.
$$

---

# 0.12 Experimental-sequence equations

A linear ramp uses

$$
y(f)=y_0+(y_1-y_0)f,
\qquad0\le f\le1.
$$

The smoothstep ramp used by the sequence engine replaces $f$ with

$$
s(f)=f^2(3-2f).
$$

After switch-off, one exponential coil response is

$$
G(t)=G_0e^{-(t-t_{\rm off})/\tau_{\rm coil}}.
$$

A residual field may include

$$
\mathbf B(t)=\mathbf B_{\rm DC}+\mathbf B_{\rm eddy}e^{-(t-t_{\rm off})/\tau_{\rm eddy}}+\mathbf B_{\rm AC}\sin(2\pi ft+\phi).
$$

---

# 0.13 Vapour, thermal-flux, and collision equations

The implemented natural-rubidium vapour-pressure fit is

$$
\log_{10}P[\mathrm{Pa}]=7.738-\frac{4215}{T}
$$

for the solid branch and

$$
\log_{10}P[\mathrm{Pa}]=7.193-\frac{4040}{T}
$$

for the liquid branch in the documented validity range.

Number density follows the ideal-gas relation

$$
n=\frac{P}{k_BT}.
$$

The one-sided equilibrium particle flux is

$$
\frac{\Phi}{A}=n\sqrt{\frac{k_BT}{2\pi m}}=\frac{n\langle v\rangle}{4}.
$$

With

$$
a=\frac{m}{2k_BT},
$$

the normalized speed distribution of particles crossing a surface is

$$
p_{\rm flux}(v)=2a^2v^3e^{-av^2}.
$$

The cosine-law incidence variable $\mu=\cos\theta$ has

$$
p(\mu)=2\mu,
\qquad0\le\mu\le1.
$$

A kinetic background-gas loss estimate is

$$
\gamma_{\rm bg}=n_{\rm bg}\sigma_{\rm loss}\langle v_{\rm rel}\rangle,
$$

with reduced mass

$$
\mu_r=\frac{m_{\rm trap}m_{\rm bg}}{m_{\rm trap}+m_{\rm bg}},
$$

and

$$
\langle v_{\rm rel}\rangle=\sqrt{\frac{8k_BT}{\pi\mu_r}}.
$$

---

# 0.14 Capture statistics and loading equations

A trajectory is classified as captured only if it remains inside the configured acceptance radius and below the configured speed threshold continuously for at least the dwell time. This is an explicit numerical definition, not a claim that the acceptance sphere is the vacuum-cell wall.

The trajectory-derived loading rate is

$$
R_{\rm load}=\Phi_{\rm incident}P_{\rm capture}.
$$

For a binomial capture estimate $\hat p=k/n$, the code uses the Wilson interval. With normal quantile $z$,

$$
p_{\rm centre}=\frac{\hat p+z^2/(2n)}{1+z^2/n},
$$

$$
\Delta p=\frac{z}{1+z^2/n}\sqrt{\frac{\hat p(1-\hat p)}{n}+\frac{z^2}{4n^2}}.
$$

The one- and two-body loading model is

$$
\dot N=R_{\rm load}-\gamma N-\frac{\beta_2}{V_{2,\rm eff}}N^2.
$$

For one-body loss only,

$$
N(t)=N_{\rm ss}+(N_0-N_{\rm ss})e^{-\gamma t},
\qquad N_{\rm ss}=\frac{R_{\rm load}}{\gamma}.
$$

With two-body loss the positive steady state is

$$
N_{\rm ss}=\frac{2R_{\rm load}}{\gamma+\sqrt{\gamma^2+4(\beta_2/V_{2,\rm eff})R_{\rm load}}}.
$$

For a Gaussian cloud with RMS widths $\sigma_x,\sigma_y,\sigma_z$,

$$
\int n^2dV=\frac{N^2}{8\pi^{3/2}\sigma_x\sigma_y\sigma_z},
$$

so

$$
V_{2,\rm eff}=8\pi^{3/2}\sigma_x\sigma_y\sigma_z.
$$

---

# 0.15 Collective-MOT equations

The Gaussian cloud density is

$$
n(\mathbf r)=\frac{N}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z}\exp\left[-\frac12\sum_i\frac{x_i^2}{\sigma_i^2}\right].
$$

The peak density is

$$
n_0=\frac{N}{(2\pi)^{3/2}\sigma_x\sigma_y\sigma_z}.
$$

The central optical depth along one axis is

$$
OD_i=\sigma_{\rm opt}\,\mathcal N_i,
$$

where $\mathcal N_i$ is the corresponding column density.

The mean-field multiple-scattering coefficient is

$$
Q=\frac{\sigma_L\sigma_R I_{\rm tot}}{4\pi c},
$$

and the radial repulsive force is modelled as

$$
F_{\rm rep}(r)=Q\frac{N_{\rm enclosed}(r)}{r^2}.
$$

The single-reabsorption proxy is

$$
P_{\rm reabs}=1-e^{-OD_R}.
$$

The associated recoil-diffusion proxy per Cartesian axis is

$$
D_{pp}^{\rm trap}=\frac{(\hbar k)^2R_{\rm sc}}{3}P_{\rm reabs}.
$$

Balancing a linear restoring force against the mean-field repulsion gives the large-MOT density scale

$$
n_{\rm lim}=\frac{3\kappa}{4\pi Q}.
$$

For a thermal harmonic cloud the single-particle width is

$$
\sigma_{\rm th}=\sqrt{\frac{k_BT}{\kappa}}.
$$

---

# 0.16 Foundation equations also present in the repository

These utilities are not the main MOT solver, but they are part of the repository and are useful sanity checks:

$$
\sigma_v=\sqrt{\frac{k_BT}{m}},
$$

$$
x(t)=x_0+v_0t+\frac12at^2,
$$

$$
z_R=\frac{\pi w_0^2}{\lambda},
\qquad w(z)=w_0\sqrt{1+(z/z_R)^2},
$$

$$
U(r,z)=-U_0\left(\frac{w_0}{w(z)}\right)^2\exp\left[-\frac{2r^2}{w^2(z)}\right],
$$

$$
\omega_r=\sqrt{\frac{4U_0}{mw_0^2}},
\qquad\omega_z=\sqrt{\frac{2U_0}{mz_R^2}},
$$

$$
\Delta z_{\rm sag}=\frac{g}{\omega^2}.
$$

---

# 0.17 Which equation should I use for which question?

| Question | Equation/model | Why |
|---|---|---|
| Is the MOT force restoring and damping? | effective scattering force | fastest transparent model |
| Where does repump population go? | 24-state population rate equations | optical pumping without coherence |
| Are quantum coherences important? | Lindblad OBE | retains off-diagonal density-matrix elements |
| What does a transverse stray field do internally? | vector 24-state OBE | full vector Zeeman mixing |
| Why can cooling go below the Doppler limit? | PGC light shifts + pumping | exposes Sisyphus mechanism |
| Is an atom captured? | Newton/RK45 trajectory | external motion is semiclassical |
| How many thermal atoms enter the capture region? | one-sided flux distribution | correct surface-crossing ensemble |
| What is the loading curve? | $\dot N=R-\gamma N-(\beta_2/V)N^2$ | separates loading and loss |
| When does independent-atom physics break down? | collective Gaussian mean field | tractable density-dependent extension |

The rest of the tutorial now derives these equations in physical order, explains why each approximation is accepted, shows the corresponding calculation, and then explains why the next model is needed.
