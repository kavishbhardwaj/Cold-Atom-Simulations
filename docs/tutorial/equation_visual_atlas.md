# Equation visual atlas

This page is a visual companion to the tutorial. The small SVG figures below are **pedagogical plots generated directly from the stated equations**; they are not experimental measurements. Repository simulation results are shown separately and are linked beside the relevant equation.

## Atomic structure

Hyperfine equation:

$$
\frac{E_{\rm hfs}}{h}=\frac{A_{\rm hfs}}{2}K+B_{\rm hfs}\frac{\frac34K(K+1)-I(I+1)J(J+1)}{2I(2I-1)J(2J-1)}.
$$

![Hyperfine levels](figures/hyperfine_energy_levels.svg)

Full vector Zeeman result from the repository:

![Exact vector Zeeman spectra](../../results/atomic_structure/exact_zeeman_spectra.png)

## Gaussian beam

$$
I(r)=I_0e^{-2r^2/w^2}.
$$

![Gaussian intensity](figures/gaussian_beam_profile.svg)

## Doppler cooling

$$
\mathbf F_i=\hbar\mathbf k_iR_i.
$$

Two counterpropagating red-detuned beams:

![Doppler damping](figures/doppler_force_vs_velocity.svg)

## MOT restoring force

A position-dependent Zeeman shift breaks the force balance:

![MOT restoring force](figures/mot_restoring_force.svg)

The full reference six-beam force map is:

![Effective MOT force map](../../results/effective_mot/force_map_x_vx.png)

## Lindblad spontaneous decay

$$
\dot\rho=-\frac{i}{\hbar}[H,\rho]+\mathcal D[C]\rho.
$$

With no drive, $\rho_{ee}=e^{-\Gamma t}$:

![Lindblad decay](figures/lindblad_spontaneous_decay.svg)

## Driven optical Bloch equations

The driven OBE produces damped Rabi oscillations:

![Rabi oscillations](figures/obe_rabi_oscillations.svg)

The steady-state excited population is

$$
\rho_{ee}=\frac{s/2}{1+s+(2\delta/\Gamma)^2}.
$$

![OBE steady state](figures/obe_steady_state_lorentzian.svg)

Independent QuTiP/PyLCP validation:

![External validation](../../results/validation/independent_software_comparison.png)

## Multilevel optical pumping

The 24-population rate equation changes both the force and internal-state distribution:

![Multilevel force](../../results/multilevel/effective_vs_multilevel_force.png)

![Manifold populations](../../results/multilevel/manifold_populations.png)

## Polarization-gradient cooling

State-dependent light shifts and pumping:

![PGC light shifts](../../results/polarization_gradient/light_shifts_pumping.png)

Velocity-dependent sub-Doppler force:

![PGC force](../../results/polarization_gradient/subdoppler_force_velocity.png)

Vector residual-field diagnostic:

![Residual field](../../results/polarization_gradient/vector_residual_obe.png)

## Thermal vapour entering a capture surface

Bulk Maxwell speeds are not the same as the flux-weighted surface-crossing distribution:

![Thermal flux](figures/thermal_flux_distribution.svg)

The actual trajectory-derived response is:

![Capture response](../../results/capture_loading/capture_response_map.png)

## Loading and loss

$$
\dot N=R_{\rm load}-\gamma N-(\beta_2/V_{2,\rm eff})N^2.
$$

![Loading equation illustration](figures/loading_loss_dynamics.svg)

Repository loading/loss scenario:

![Loading result](../../results/capture_loading/loading_loss_sensitivity.png)

## Experimental timing and collective physics

![Sequence](../../results/sequence/sequence_timeline.png)

![Collective MOT](../../results/collective_mot/collective_mot_diagnostics.png)

The visual atlas is intentionally paired with the [notation/equation inventory](00_notation_and_equation_inventory.md): the atlas shows what the equations do, while the inventory defines every symbol and records the exact equation used by each solver.
