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
