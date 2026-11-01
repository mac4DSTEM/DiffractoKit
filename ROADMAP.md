# Roadmap

Where DiffractoKit is going next, in order. `docs/parity.md` says what is held today and its known limits, and the GitHub
issues track each piece of work. The numbers below are the ones `docs/parity.md` cites as "ROADMAP item n".

1. **Diagnose before fixing.** Phase matching (`minimumVectors` = 2 makes "matrix" by exclusion; the cross-phase winner
   ignores completeness), the ACOM mirror (a mirrored win carries `inPlaneAngle` + π with no flag; some off-grid
   self-recoveries land 36 or 60 degrees away) and CIF symmetry (a group named by Hermann–Mauguin symbol alone, a partial
   operation list). Each gets its cause stated first; what the diagnosis supports is fixed, the rest stays an open item.
2. **Space-group facts, public.** Space group symbol, number and Hall symbol through spglib's definitions (no Wyckoff
   positions first); a public `ReflectionFamily` initialiser; `RingIndexing`'s "accepted" documented as a likelihood ratio,
   not family membership, or changed.
3. **orix parity beyond the IPF key.** `Rotation` and `Orientation` types, misorientation, the orientation fundamental
   zone and its sampling, pole-figure points (stereographic projection), and the `.native` IPF key for every Laue class.
4. **Kinematic diffraction patterns.** A diffsims-like pattern API (the spots of a beam direction in one call) with
   Debye–Waller factors, and a streamed q-calibration reflection list (today it peaks near 800 MB for a 30 Å cell at
   k_max 2.5 Å⁻¹).
5. **py4DSTEM ACOM extras, only if needed.** `match_single_pattern` refinement and `num_matches_return` > 1; until then the
   parity table says the scores are not comparable.

Not planned: CIF export, grain clustering, `.ang` export, database CIF sources (COD, ICSD), plotting.
