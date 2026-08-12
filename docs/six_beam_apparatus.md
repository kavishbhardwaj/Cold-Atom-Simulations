# Six physical laser beams

`GaussianBeam` represents one physical travelling wave. The six-beam MOT is
therefore a collection of six independently replaceable objects, not one global
power, polarization or frequency. Direction and origin encode pointing and
offset; `waist`/`waist_y` encode an elliptical 1/e2 intensity profile; focus,
wavelength, power, detuning, AOM offset, optical phase, linewidth and coherence
group are independent. `propagation_mode="gaussian"` includes `w(z)`, separate
Rayleigh ranges, wavefront curvature and Gouy phase. The default `collimated`
mode preserves the previous fast approximation.

An explicit two-component Jones input and ordered `JonesElement` train support
linear polarizers, quarter/half-wave plates, arbitrary retarders, plate angle,
retardance error and viewport birefringence (an arbitrary retarder). The output
is normalized and mapped into the beam's transverse laboratory basis. Laser
linewidth is angular FWHM and is active as Lorentzian optical-coherence
broadening in effective and rate-equation scattering models.

`SixBeamApparatus` supports independent beams, three assigned pairs, and three
retroreflected inputs. `Retroreflection` applies power loss, phase, mirror tilt,
Jones change/double-pass elements, and a coherence group. Selected pairs may be
coherent standing waves while other pairs remain mutually incoherent.

## Polarization language

A MOT contains three counterpropagating pairs. Calling them globally “three
sigma+ and three sigma-” is incomplete: sigma+, pi and sigma- are defined with
respect to a **local quantization axis**, not merely propagation. For the ideal
apparatus decomposed relative to laboratory z, each transverse circular beam is
25% sigma-, 50% pi, 25% sigma+, while the z beams are pure opposite spherical
components. `local_polarizations(axis)` calculates the actual fractions for
every beam and any supplied local axis.

## Controlled sensitivity scan

`generate_six_beam_apparatus_results.py` stores a per-beam CSV, a 3-D apparatus
diagram, and one-at-a-time effective-MOT scans. In its specified recipe, +/-10%
power error in one x beam displaces the solved centre by roughly -0.74 to
+0.69 mm. A 5 mrad pointing error displaces it by about 36 micrometres. The QWP
angle sweep changes the sampled restoring force continuously but modestly. These
are effective-model apparatus sensitivities, not coherent PGC, capture, or
experimental uncertainty estimates; simultaneous correlated imperfections can
behave differently.
