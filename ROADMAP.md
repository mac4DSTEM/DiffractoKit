# Roadmap

Where DiffractoKit is going next, in order. `docs/parity.md` says what is held today and its known limits, and the GitHub
issues track each piece of work. The numbers below are the ones `docs/parity.md` cites as "ROADMAP item n".

1. **Diagnose before fixing.** A number that reads wrong gets its cause stated first; what the diagnosis supports is
   fixed, the rest stays an issue. Open today: the ACOM beam direction lands 3 to 6 degrees off the [001] and [011] zone
   axes on the AlMgSi demo where excitation widths nearly tie (#65). The phase-matching, ACOM-mirror and CIF-symmetry
   items of October are closed (#7, #8, #9). Listing on the Swift Package Index (#58) is open beside it.
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
