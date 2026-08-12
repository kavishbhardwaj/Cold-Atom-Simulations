# Model hierarchy and capability boundaries

The package exposes complementary approximations, not chronological releases.
The effective model is fast and interpretable; population rate equations add
hyperfine pumping; OBEs retain coherence; the adiabatic polarization-gradient
model resolves phase-dependent dressed potentials but not ground coherences.
None is automatically experimentally calibrated.

The two-level OBE remains an analytical benchmark. The sparse multilevel OBE
constructs a generated hyperfine/Zeeman Hamiltonian and spontaneous-collapse
channels, but full six-beam steady-state grids are research mode because the
24-state 87Rb D2 density matrix has 576 complex components. D1 gray molasses is
explicitly unsupported until a Λ-system Raman benchmark is independently
validated.

`quick` means CI-scale deterministic grids. `research` means users must refine
basis, optical phases, OBE tolerances, spatial grids, atom count and timestep.
