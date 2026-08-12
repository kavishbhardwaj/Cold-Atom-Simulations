# Part IV — Cross-cutting interpretation, approximations, and validation

The preceding chapters now discuss each result immediately after the equations and decisions that produced it. This chapter therefore serves as a cross-cutting guide: how the pieces connect, which approximations were accepted, how validation language is used, and what conclusions remain outside the present model.

---
# 16. The reference simulation flow, step by step

A student can think of one standard calculation as the following pipeline.

## Step 1: load atomic data

Read 87Rb D2 mass, wavelength, lifetime, hyperfine constants, angular momenta and saturation intensity.

Compute:

$$
k=\frac{2\pi}{\lambda},\quad
\Gamma=1/\tau,\quad
v_r=\hbar k/m,\quad
T_D=\hbar\Gamma/(2k_B).
$$

## Step 2: build the hyperfine/Zeeman basis

Generate $|F,m_F\rangle$ states, hyperfine offsets, $g_F$, dipole strengths and spontaneous branching.

## Step 3: construct six beams

For each beam calculate Gaussian intensity, saturation, propagation vector and polarization components.

## Step 4: construct the magnetic field

Use either the ideal quadrupole or a physical coil model plus residual/bias fields.

## Step 5: choose fidelity

- **effective model:** fast force/trajectory/capture calculation;
- **rate equation:** populations and repump without coherence;
- **24-state OBE:** coherences, vector Zeeman physics and beam-resolved moving-atom force;
- **reduced PGC:** intuitive Sisyphus mechanism;
- **collective model:** density-dependent extension after a single-atom model is chosen.

## Step 6: compute internal state and force

Depending on the selected model, calculate scattering rates, stationary populations or density matrix, then obtain $\mathbf F(\mathbf r,\mathbf v,t)$.

## Step 7: propagate motion

Integrate Newton's equations deterministically or with photon-event recoil.

## Step 8: classify capture

For vapour loading, launch atoms from the acceptance surface using the flux-weighted thermal distribution, integrate each trajectory, and apply the explicit capture criterion.

## Step 9: calculate loading

$$
R=\Phi_\mathrm{incident}P_\mathrm{capture}.
$$

Then integrate

$$
\dot N=R-\gamma N-\beta\int n^2dV.
$$

## Step 10: add experimental timing if required

Replace static laser/field settings with ramps, coil decay, eddy currents and a molasses/TOF timeline.

## Step 11: interpret at the correct fidelity

A force curve is a force curve. A capture rate is not automatically a calibrated atom number. A friction coefficient is not automatically a temperature unless diffusion is calculated at compatible fidelity.

That final interpretation rule is one of the most important design principles of the project.

---
# 17. Major approximations and why they were accepted

| Approximation / decision | Why it was used | What it prevents us from claiming |
|---|---|---|
| Semiclassical external motion | MOT de Broglie scale is much smaller than apparatus scale; efficient for trajectories | No fully quantum motional wavepacket dynamics |
| Effective cycling-transition force | Fast, transparent baseline for apparatus/capture scans | No hyperfine optical pumping or PGC |
| Stationary rate equations | Efficient multilevel population solution when internal relaxation is fast | No coherences, dark states or Sisyphus interference |
| Block-rotating 24-state OBE | Removes numerically useless GHz cooling-repump beat | Very far-off-resonant cross-ground terms are treated by controlled RWA |
| Secular cross-F magnetic approximation in OBE block frame | Valid for weak MOT fields relative to GHz hyperfine splitting | Must be revisited toward hyperfine Paschen-Back fields |
| Phase averaging for incoherent beams | Real unrelated beams should not possess one arbitrary fixed phase | Requires convergence testing over phase ensembles |
| Reduced population PGC model | Gives an interpretable Sisyphus mechanism at modest cost | No quantitative coherent PGC temperature |
| Isotropic spontaneous recoil | Correct zero mean and simple diffusion baseline | Does not include complete dipole radiation pattern/internal force noise |
| Spherical capture surface | Separates capture dynamics from unspecified chamber geometry | Not a full vapour-cell wall/dispenser model |
| User/literature collision parameters | Avoids invented atom-number predictions | Absolute N requires apparatus-specific loss data |
| Gaussian collective cloud | Efficient mean-field density model | No exact radiative transfer or arbitrary cloud deformation |

---
# 18. Validation philosophy

The repository uses three different words deliberately:

**Implemented**: code exists and passes internal tests.

**Externally verified**: a matched calculation agrees with an independent package such as QuTiP or PyLCP.

**Experimentally calibrated**: measured apparatus inputs and quantitative agreement with a specified experiment support the prediction.

The current strongest external validations are the two-level OBE, normalized two-beam force and 87Rb ground Zeeman spectrum.

The principal high-fidelity validation still pending is a matched PyLCP comparison of the full 87Rb D2 multilevel force/populations, followed by a quantitatively matched polarization-gradient benchmark and a fully specified 87Rb MOT/PGC experiment.

---
# 19. What the repository deliberately does not claim yet

The code should **not yet** be interpreted as giving a universal or calibrated answer for:

- the final sub-Doppler temperature of a real 87Rb apparatus;
- the exact temperature increase caused by a given stray magnetic field;
- the absolute atom number of an arbitrary vapour-cell MOT;
- chamber-wall recycling/desorption;
- exact multiple-scattering radiative transfer;
- D1 gray molasses and Raman dark states;
- a complete high-field Paschen-Back model;
- a laboratory digital twin without measured beam, coil, vacuum and timing inputs.

These are not hidden shortcomings. They are explicit boundaries used to keep the simulation scientifically interpretable.

---
