# Tutorial: From Rubidium Atoms to a MOT Simulation

This documentation is written to teach the physics **as a worked calculation**, not to list APIs.

> **Start here:** [**A Worked Textbook Derivation of the MOT Model**](textbook_derivation.md)
>
> This is now the primary learning document. It constructs the calculation from atomic constants through the explicit 24-state basis, the full \(24\times24\) density matrix, the \(576\times576\) sparse Liouvillian, optical forces, sub-Doppler physics, trajectories, vapour loading, and collective effects. Every major symbol is defined before it is used and the numerical object solved by the code is shown explicitly.

A shorter narrative version is available as the [continuous equation-to-result walkthrough](continuous_walkthrough.md).

Three companions make the textbook easier to audit and visualize:

- [**Notation and complete governing-equation inventory**](00_notation_and_equation_inventory.md) — symbol glossary plus the core physical equations that directly determine repository results. It distinguishes easy-to-confuse symbols such as \(A_{\rm hfs}\) versus the rate-generator matrix, nuclear spin \(I\) versus optical intensity, and damping \(\beta_v\) versus two-body loss \(\beta_2\).
- [**Equation visual atlas**](equation_visual_atlas.md) — pedagogical plots for the hyperfine equation, Gaussian beam, Doppler damping, MOT restoring force, Lindblad decay, Rabi oscillations, OBE steady state, thermal flux and loading, alongside the actual generated repository results.
- [**Scientific results gallery**](../../results/README.md) — the configured simulation outputs and their fidelity/provenance notes.

> **Scope / provenance.** This repository is independent after-hours work developed from personal scientific interest and kept as a reproducible record and backup. Laboratory control, acquisition, and other lab codes are not kept here.

## What the textbook means by “show the full calculation”

For the 87Rb D2 coherent model, the textbook does not stop at the phrase “24-state OBE.” It explicitly:

1. lists the 8 ground and 16 excited \(|F,m_F\rangle\) basis states in the same order used by the code;
2. shows the exact \(24\times24\) density-matrix block structure;
3. explains populations, ground coherences, excited coherences, and optical coherences;
4. constructs the block Hamiltonian \(h=H/\hbar\);
5. derives the beam-resolved Rabi couplings and moving optical phases;
6. constructs every spontaneous-emission collapse operator;
7. vectorizes \(\rho\) into 576 components;
8. derives the \(576\times576\) sparse Liouvillian;
9. explains stationary and time-dependent solution strategies;
10. derives the Hamiltonian-gradient force used by the multilevel OBE.

Printing 576 individual symbols would obscure rather than clarify the physics, so the textbook gives the complete state ordering and indexing rule \(\rho_{ij}=\langle i|\rho|j\rangle\). That uniquely defines every matrix element while keeping the physical block structure visible.

## Deeper reference chapters

After the textbook, these chapters provide more detail on individual subsystems:

1. [Part I — Physical system and apparatus](01_physical_system_and_apparatus.md) — \(^{87}\)Rb D2 choice, hyperfine/Zeeman basis, fields, coils, six physical beams, polarization and coherence.
2. [Part II — MOT force models, OBEs, and sub-Doppler physics](02_mot_force_models_and_obes.md) — effective MOT force, multilevel rate equations, two-level and 24-state OBEs, PGC and residual magnetic fields.
3. [Part III — Motion, loading, and collective physics](03_motion_loading_and_collective_physics.md) — RK45 motion, recoil, time sequence, vapour flux, capture, loading/loss and collective effects.
4. [Part IV — Cross-cutting interpretation, approximations, and validation](04_results_validation_and_scope.md) — complete pipeline, approximation boundaries and claim language.
5. [Part V — Digital twin, learning, and reproduction](05_digital_twin_learning_and_reproduction.md) — experiment-specific inputs, clean-checkout reproduction, equation-to-code map and reproducibility checklist.

## The modelling principle

The repository intentionally does **not** use one monolithic “exact MOT model.” It uses a hierarchy:

**atomic constants and basis** → **six beams + magnetic field** → **effective scattering force** → **multilevel rate equations** → **coherent OBEs** → **polarization-gradient physics** → **classical motion/recoil** → **vapour capture/loading** → **time sequence** → **collective-cloud effects**.

Each layer exists because it answers a different physical question at a different computational cost. A result is interpreted only at the fidelity of the model that generated it.

## What “all equations” means here

The textbook and [equation inventory](00_notation_and_equation_inventory.md) include the **governing physical equations that directly determine the scientific outputs**, together with numerical/statistical definitions that materially affect those outputs, such as phase averaging, rotating-frame approximations, capture statistics, loading/loss and sequence transients.

They do not promote trivial array reshaping, plotting coordinates, or ordinary software bookkeeping into “physics equations.” If an implementation equation can materially change a scientific result, it belongs in the inventory and is explained in the textbook or a reference chapter.
