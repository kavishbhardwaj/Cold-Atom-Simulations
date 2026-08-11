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
The registry records the approximate natural fractions (85Rb 72.17%, 87Rb
27.83%) so that future loading work cannot silently assume equal abundance.
Those fractions are provenance data only in Phase 3; they are not yet used to
predict loading.

## Why D1 is not interchangeable with D2

D1 is 5S1/2→5P1/2 near 795 nm; D2 is 5S1/2→5P3/2 near 780 nm. They differ in
excited-state angular momentum, hyperfine manifolds, Landé factors, branching,
line strength and wavelength. A D2 transition table cannot be relabeled D1.
D1 is valuable for gray molasses and Λ-system coherence, but those effects are
precisely where population-only rates are inadequate. Credible D1 support needs
its own hyperfine data plus the Level-C multilevel Hamiltonian and decay graph.
Using D1 and D2 simultaneously is physically meaningful for specialized cooling,
state preparation or probing, but it creates a multi-frequency, multi-manifold
Hamiltonian. A conventional 87Rb vapour-cell MOT normally uses D2 cooling and a
D2 repump; adding D1 is therefore a new experiment model, not an automatic
accuracy upgrade.

## Why 85Rb is not a parameter substitution

85Rb has nuclear spin I=5/2 rather than 87Rb I=3/2, a 3.036 GHz ground splitting,
and ground manifolds F=2,3 rather than F=1,2. Its D-line hyperfine states,
Clebsch–Gordan graph, cooling/repump assignments and off-resonant leakage must be
rebuilt. Mass and wavelength replacement alone would be scientifically wrong.

## Support matrix

| Isotope/line | Constants registry | Level A | Level B | Level C | Level D |
|---|---:|---:|---:|---:|---:|
| 87Rb D2 | yes | yes | 24-population D2 | reduced stretched transition | adiabatic F=2→F′=3 PGC |
| 87Rb D1 | yes | not configured | not implemented | generic reduced two-level | not implemented |
| 85Rb D2 | yes | not configured | not implemented | generic reduced two-level | not implemented |
| 85Rb D1 | yes | not configured | not implemented | generic reduced two-level | not implemented |

“Constants registry” means wavelength, lifetime, mass, nuclear spin and ground
splitting have provenance. It does not mean a solver silently supports the line.
`AtomicLine.model_support` exposes that distinction programmatically.
The reduced two-level OBE can use any registered lifetime and detuning because
it deliberately discards hyperfine structure; this is useful for analytical
line-scale comparisons but is not a D1/85Rb MOT implementation.
The registry also derives wave number, recoil velocity, recoil temperature and
two-level Doppler temperature from each line's wavelength, lifetime and isotope
mass. These comparisons help expose scale changes without implying dynamic
support.

## Practical comparison

| System | Typical role relevant here | What must change from the reference |
|---|---|---|
| 87Rb D2 | Standard F=2→F'=3 MOT cooling; F=1→F'=2 repump | Implemented Level A/B reference |
| 87Rb D1 | Λ gray molasses, Raman/EIT-style preparation and probing | 5P1/2 F'=1,2 basis, D1 strengths, coherent multi-frequency OBE |
| 85Rb D2 | Standard F=3→F'=4 MOT cooling; F=2→F'=3 repump | I=5/2 basis, 85Rb hyperfine offsets, strengths and branching |
| 85Rb D1 | 85Rb D1 molasses/coherent preparation | Both the isotope basis and D1 coherent couplings change |

D2 is especially convenient for a conventional alkali MOT because its upper
fine-structure manifold contains the `F→F+1` stretched cycling transition. On
D1 the maximum excited F is not one larger than the maximum ground F, so there
is no analogous closed stretched cycling line. D1 can still exert trapping and
cooling forces, but leakage and coherent Λ physics are central rather than a
small correction. This is why “support both lines” requires a genuine model
extension, not a wavelength switch.

## Cohesive extension path

1. Add isotope/line-specific hyperfine offsets, Landé factors and reduced line
   strengths from the corresponding Steck/NIST source.
2. Generate the basis and dipole graph from nuclear/electronic angular momentum,
   rather than duplicating 87Rb D2 tables.
3. Validate branching and line sums for that exact line.
4. Define its cooling/repump or molasses configuration explicitly.
5. Compare rate equations and OBEs under matched assumptions before adding the
   system to a public “supported” table.
