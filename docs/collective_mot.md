# Optional collective MOT mean field

The validated single-atom force is unchanged. `physics.collective` adds an
optional Gaussian continuum described by atom number, three RMS widths,
temperature/velocity widths, peak density and central optical depths. The exact
Gaussian integral `integral n^2 dV=N^2/(8 pi^(3/2) sigma_x sigma_y sigma_z)`
drives the light-assisted two-body term. Background-gas and hot-Rb losses remain
separate one-body inputs; beta is rejected unless accompanied by an explicit
experimental/literature provenance string.

The first multiple-scattering model follows the Coulomb-like approximation of
Walker, Sesko, and Wieman, *Phys. Rev. Lett.* **64**, 408 (1990), DOI
`10.1103/PhysRevLett.64.408`: `F=Q N_enclosed/r^2`, with
`Q=sigma_L sigma_R I/(4 pi c)`. Both cross sections and intensity are physical
inputs, not an arbitrary repulsion constant. Balancing this force with a linear
MOT restoring force gives the familiar large-MOT constant-density limit. A
Gaussian thermal width is used at small N and the larger density-limited uniform-
sphere RMS width at large N.

Beer-Lambert central optical depth supplies attenuation and a shadow-force
scale. `1-exp(-OD_R)` estimates single-photon reabsorption probability and an
additional isotropic recoil diffusion. This distinguishes direct laser
scattering from radiation trapping, but it is **not exact radiative transfer**:
frequency redistribution, polarization, multilevel cross sections, repeated
scattering, anisotropic escape and cloud deformation are omitted.

`CollectiveLoading` integrates
`dN/dt=R-(gamma_background+gamma_hotRb)N-beta integral(n^2)dV` while updating
the quasi-static radius and peak density. QUICK and RESEARCH modes differ only
in ODE tolerances. With beta=0 and no scattering mean field it reduces to the
existing independent-atom loading curve; this is regression tested.

The committed result uses explicitly labeled user scenarios for temperature,
restoring coefficient, beta and reabsorption cross section. It demonstrates
loading saturation, cloud expansion, density limitation, optical depth,
repulsion, shadowing and radiation-trapping diffusion. It does not infer a
temperature from the old PGC proxy. Quantitative comparison to Walker et al. is
limited to the constant-density/expansion trend because the reference apparatus
does not reproduce their measured intensities, cloud shape or effective
reabsorption cross section. Status: **LITERATURE-TREND VERIFIED**, not
experimentally validated.
