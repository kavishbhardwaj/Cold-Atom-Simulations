# Atomic-system scope: isotope and D line

## Current choice: 87Rb D2

The implemented MOT uses **87Rb 5S1/2→5P3/2 (D2), 780.241209686 nm**. Cooling
is referenced to F=2→F'=3 and repumping to F=1→F'=2. This is not because other
rubidium systems are unimportant. It is a controlled first target because:

- F=2→F'=3 provides a strong stretched-state cycling transition;
- the 6.835 GHz ground splitting cleanly separates cooling and repump lasers;
- the isotope and transition are widespread in teaching and research MOTs;
- authoritative constants and hyperfine strengths are tabulated by Steck;
- one complete isotope/line is more auditable than four partially mixed models.

Natural rubidium contains both stable isotopes, so a real vapour cell can show
85Rb and 87Rb spectral features. Isotope selectivity comes from laser frequency,
not from the vapour containing only one isotope. Phase 5 loading must therefore
allow isotope abundance or enriched-source composition explicitly.

## Why D1 is not interchangeable with D2

D1 is 5S1/2→5P1/2 near 795 nm; D2 is 5S1/2→5P3/2 near 780 nm. They differ in
excited-state angular momentum, hyperfine manifolds, Landé factors, branching,
line strength and wavelength. A D2 transition table cannot be relabeled D1.
D1 is valuable for gray molasses and Λ-system coherence, but those effects are
precisely where population-only rates are inadequate. Credible D1 support needs
its own hyperfine data plus the Level-C multilevel Hamiltonian and decay graph.

## Why 85Rb is not a parameter substitution

85Rb has nuclear spin I=5/2 rather than 87Rb I=3/2, a 3.036 GHz ground splitting,
and ground manifolds F=2,3 rather than F=1,2. Its D-line hyperfine states,
Clebsch–Gordan graph, cooling/repump assignments and off-resonant leakage must be
rebuilt. Mass and wavelength replacement alone would be scientifically wrong.

## Support matrix

| Isotope/line | Constants registry | Level A | Level B | Level C |
|---|---:|---:|---:|---:|
| 87Rb D2 | yes | yes | 24-population D2 | reduced stretched transition |
| 87Rb D1 | yes | not configured | not implemented | not implemented |
| 85Rb D2 | yes | not configured | not implemented | not implemented |
| 85Rb D1 | yes | not configured | not implemented | not implemented |

“Constants registry” means wavelength, lifetime, mass, nuclear spin and ground
splitting have provenance. It does not mean a solver silently supports the line.
`AtomicLine.model_support` exposes that distinction programmatically.

## Cohesive extension path

1. Add isotope/line-specific hyperfine offsets, Landé factors and reduced line
   strengths from the corresponding Steck/NIST source.
2. Generate the basis and dipole graph from nuclear/electronic angular momentum,
   rather than duplicating 87Rb D2 tables.
3. Validate branching and line sums for that exact line.
4. Define its cooling/repump or molasses configuration explicitly.
5. Compare rate equations and OBEs under matched assumptions before adding the
   system to a public “supported” table.
