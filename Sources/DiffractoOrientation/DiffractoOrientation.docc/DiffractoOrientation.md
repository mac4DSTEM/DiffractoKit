# ``DiffractoOrientation``

Reference libraries and vector matching for phase maps, and the value types of an orientation result.

## Overview

DiffractoOrientation builds, from the kinematic reflections of DiffractoScattering, the reciprocal-lattice vectors a
detector would show for each phase, zone axis and in-plane rotation (`PhaseReferenceLibrary`), matches measured diffraction
vectors against them to map phases over a scan (`PhaseVectorMatcher`), and holds the value types of an orientation
result: Euler angles in py4DSTEM's and SciPy's convention (`EulerAngles`), the detector frame of a zone axis
(`ZoneAxisFrame`) and one position's orientation result (`OrientationResult`).

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
