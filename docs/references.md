# Scientific references

Only sources actually used in Phase 1 are listed here.  URLs are supplied so
that constants and equations can be audited.

1. D. A. Steck, *Rubidium 87 D Line Data*, revision 2.3.2 (2021),
   <https://steck.us/alkalidata/rubidium87numbers.pdf>. Atomic mass, D2 vacuum
   wavelength, lifetime, linewidth, saturation-intensity convention, hyperfine
   intervals, Landé factors, and relative hyperfine strengths.
2. NIST, *Atomic Spectra Database*, Rb I,
   <https://physics.nist.gov/PhysRefData/ASD/>. Independent level/wavelength
   provenance; no claim of a new fit to NIST data is made here.
3. E. L. Raab et al., “Trapping of Neutral Sodium Atoms with Radiation
   Pressure,” *Physical Review Letters* **59**, 2631 (1987),
   <https://doi.org/10.1103/PhysRevLett.59.2631>. Six-beam MOT principle.
4. H. J. Metcalf and P. van der Straten, *Laser Cooling and Trapping*
   (Springer, 1999), <https://doi.org/10.1007/978-1-4612-1470-0>. Effective
   two-level scattering force, Doppler cooling, recoil, and MOT linearization.
5. J. Dalibard and C. Cohen-Tannoudji, “Laser cooling below the Doppler limit
   by polarization gradients,” *JOSA B* **6**, 2023 (1989),
   <https://doi.org/10.1364/JOSAB.6.002023>. Why Phase-1 radiation pressure
   cannot represent polarization-gradient or Sisyphus cooling.
6. S. E. Galica, L. Aldridge, and E. E. Eyler, “PyLCP: A Python package for
   computing laser cooling physics,” *New Journal of Physics* **23**, 065002
   (2021), <https://doi.org/10.1088/1367-2630/abf9d8>. Planned independent
   rate-equation/OBE reference backend; Phase 1 does not claim PyLCP validation.
7. D. J. Griffiths, *Introduction to Electrodynamics*, 4th ed. (Pearson,
   2013), Sec. 5.2. Biot–Savart law.  The implementation here performs its own
   converged midpoint segmented-wire quadrature and copies no third-party code.
8. B. P. Anderson and M. A. Kasevich, “Loading a vapor-cell magneto-optic trap
   using light-induced atom desorption,” *Physical Review A* **63**, 023404
   (2001), <https://doi.org/10.1103/PhysRevA.63.023404>. Context for the future
   loading model only; Phase 1 predicts neither loading rate nor atom number.
9. A. R. Edmonds, *Angular Momentum in Quantum Mechanics*, 2nd ed. (Princeton
   University Press, 1960), <https://doi.org/10.1515/9781400884186>. Racah
   factorial convention used to generate the Phase-2 Clebsch–Gordan factors.
10. D. A. Steck, *Rubidium 85 D Line Data*, revision 2.3.2 (2021),
    <https://steck.us/alkalidata/rubidium85numbers.pdf>. 85Rb D1/D2 wavelength,
    lifetime, mass and hyperfine provenance used by the atomic-line registry.
11. G. Lindblad, “On the generators of quantum dynamical semigroups,”
    *Communications in Mathematical Physics* **48**, 119–130 (1976),
    <https://doi.org/10.1007/BF01608499>. Master-equation dissipator used by the
    reduced Phase-3 optical-Bloch backend.
