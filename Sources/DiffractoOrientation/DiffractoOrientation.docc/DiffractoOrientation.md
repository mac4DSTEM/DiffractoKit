# ``DiffractoOrientation``

Reference libraries and vector matching for phase maps, and the value types of an orientation result.

## Overview

DiffractoOrientation builds, from the kinematic reflections of DiffractoScattering, the reciprocal-lattice vectors a
detector would show for each phase, zone axis and in-plane rotation (`PhaseReferenceLibrary`), matches measured diffraction
vectors against them to map phases over a scan (`PhaseVectorMatcher`), and holds the value types of an orientation
result: Euler angles in py4DSTEM's and SciPy's convention (`EulerAngles`), the detector frame of a zone axis
(`ZoneAxisFrame`) and one position's orientation result (`OrientationResult`).

## ACOM template matching

``ACOMTemplateLibrary`` builds, from a crystal's kinematic reflections, a polar template for every beam direction sampled in
the fundamental zone of its Laue class; ``ACOMMatcher`` correlates a measured pattern's peaks with all of them (the in-plane
rotation from one FFT along the azimuth, a second pass for the mirror half of the zone) and returns an ``OrientationResult``.
This is the simplified route of the suite's apps, not py4DSTEM's SO(3) correlogram (there is no sub-pixel refinement, shell
grouping or multi-match). Choose the excitation error with ``ACOMExcitationModel``: `.py4DSTEM(wavelengthAngstrom:)` for patterns
in py4DSTEM's orientation convention, `.appCurvature` for the apps', `.flat` when the in-plane angle may be ambiguous by a half
turn. docs/parity.md (P25) lists what is pinned to the app source and to py4DSTEM.

A whole scan is matched with ``ACOMMatcher/matchAll(peaks:width:height:originX:originY:invAngstromPerPixel:selection:maxWorkers:isCancelled:progress:)``:
the positions are dealt over worker threads, ``ACOMScanSelection`` picks all of them, one per block (a preview) or a square,
`isCancelled` is polled before every position and `progress` reports the finished fraction. The result is an ``ACOMScanResults`` grid.
`tools/acom-benchmark/run.sh` times it.

## Polarity

A kinematic pattern cannot tell +c from -c on a polar structure. Without anomalous scattering |F(hkl)|^2 equals
|F(-h-k-l)|^2 (Friedel's law), so the intensities of a structure and of its inversion are identical, in this kit and in
py4DSTEM. A kinematic match therefore fixes an orientation only up to this inversion: no consumer may claim a polarity (the
sense of a polar axis, or which of a pair of enantiomorphs is present) from a kinematic match, however good the match. The
orientation and phase results say which reflections matched, not which way a polar crystal points.

## Symmetry and IPF keys

``CrystalOrientationSymmetry`` constructs the proper subgroup of every ``DiffractoCore/LaueClass`` in its documented
standard setting. ``CrystalOrientationSymmetry/reduce(_:)`` returns one minimum-angle rotation representative;
``CrystalOrientationSymmetry/sampleFundamentalZone(count:)`` samples directions deterministically.

Choose an ``IPFColorKey`` explicitly. `.native` folds a pole using the original m-3m, 6/mmm or b-unique 2/m key;
other native keys throw a named error. `.py4DSTEM` expects the already reduced Python family pole in the automatic
sector (the original pole for 2/m), and reproduces the pinned full-mask beam key at Float precision. Python's
independent direction-family reduction is not ported; it returns three families rather than one rotation.
