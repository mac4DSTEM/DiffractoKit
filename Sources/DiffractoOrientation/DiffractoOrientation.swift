//
//  DiffractoOrientation.swift
//  Role: the module overview. DiffractoOrientation holds (registration docs/registration-slice-4.md):
//
//    - `PhaseReferenceLibrary` and its value types (`PhaseDefinition`, `PhaseReferenceSettings`, `PhaseOrientationReference`,
//      `ReferenceVector`, `LatticeVector`, `OrientationRelationship`): for each phase, zone axis and in-plane rotation, the
//      reciprocal-lattice vectors a detector would show, from the kit's kinematic reflections;
//    - `PhaseVectorMatcher` and its value types (`PhaseVectorSettings`, `PhaseVectorResolution`, `PhaseVerdict`,
//      `PhaseVectorResult`, `PhaseMap`, `PhaseDiskClaim`, `PeakCoordinates`): the vector-matched phase map of a scan, after
//      Thronsen et al., Ultramicroscopy 255 (2024) 113861;
//    - `EulerAngles` (py4DSTEM's orientation matrix to and from SciPy's Bunge angles), `ZoneAxisFrame` (the detector frame of
//      a zone axis) and `OrientationResult` (one position's orientation result).
//
//  Planned, not here (slice 5, against pyxem / orix and py4DSTEM's crystal orientation mapping): template matching,
//  symmetry reduction and fundamental zones per Laue class, IPF colour keys.
//  Depends on DiffractoCore and DiffractoScattering.
//
