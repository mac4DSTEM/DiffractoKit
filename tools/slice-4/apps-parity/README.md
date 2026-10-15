# slice 4 apps-parity (DiffractoKit, tier 1)

`run.sh` compiles the app's own vector-matched phase mapping from mac4DSTEM 7d6fdbd1 (`git show`, read only), runs it on the
phases and planted scans of `main.swift`, and writes every number as bit patterns to
`Tests/DiffractoOrientationTests/Fixtures/slice-4-apps.json` (the -O build, the app's Release). The same sources at -Onone are dumped to
`.build/apps-parity-4/dump-Onone.json` and compared. `run.sh --hash-check` also runs -O once with a random Dictionary seed.
About 1 to 4 minutes depending on machine load (the app's `map` uses every core), 13.5 MB. The dump holds no clock and no path: it is deterministic.

## What is compiled, and what is cut out

Whole files, verbatim: `Crystal.swift`, `ScatteringFactors.swift`, `PhaseReferenceLibrary.swift`, `PhaseVectorMatching.swift`,
`AnalysisCancellationToken.swift`. No hand-written shim holds numeric logic: the three missing types are cut out of larger app
files by line range, verbatim, and the first line of each range is checked by `run.sh`.

| Cut out | From | Why not the whole file |
|---|---|---|
| `BraggPeak` (34-45), `BraggVectors` (689-705, closing brace added) | `DiskDetection.swift` | imports Accelerate and Metal; `calibrated(with:)` needs the whole analysis layer |
| `EulerAngles`, `ACOMOrientation` (1-2, 12-107) | `OrientationResult.swift` | the rest is `OrientationMap`, the symmetries and the IPF key (slice 5) |

## Dump layout

`D(x)` = `{"v": decimal|null, "b": 16-hex Double bits}`, `F(x)` the same with 8-hex Float bits. `reflections` (the app's
`Crystal.reflections(kMax: 1.6)` in its order), `libraries` (every entry; an entry's vectors are `[h,k,l,qx,qy,length,rel]` bit
patterns, `vd` the decimals at rotation 0, `hash` an FNV-1a 64 over every field of the entry: see `entryHash` in `main.swift`; the
`main` library dumps all 2520 entries in full, `al-all-axes` and `oblique` only the in-plane steps k = 0, 1, 17, 90, 179 in full and
all of them by hash), `scans` (every input vector as Float pixel bit patterns, origin, scale, planted truth, and the matcher's
experimental vectors), `runs` (every `PhaseVectorMatcher.map` result), `fieldEffects`, `matrixFits`, `zoneAxisFits`, `claims`,
`resolution`, `inPlaneSteps`, `nearest`, `euler`, `detectorBases`, `lowIndexZoneAxes`.

Libraries: `main` = Al (matrix) + Au + θ′ + β″ at the default `PhaseReferenceSettings`, zone axes [001] [011] [111] (θ′ also [100], β″ also
[010]); `main-or` the same with the orientation relationship (002)θ′ ∥ (200)Al on θ′ (its entries are bit-identical to `main`'s:
only `map`'s candidate filter moves); `al-all-axes` Al with `zoneAxes: []` (maximumEntries raised to 20000) beside Au [001]; `oblique` Al, Mg
and a triclinic cell. θ′ is the six-site cell of the app's own orientation-relationship test (a = b = 4.04, c = 5.80 Å).

Scans: (a) 16 x 16 Al [001] at 6 degrees with θ′ [100] at the orientation-relationship angle, θ′ [011], β″ [010] and [001], Au [011]
planted by region, uniform 0.004 Å⁻¹ jitter, two random vectors at some positions, an empty position and a direct-beam-only
position; (b) 3-12 random vectors per position; (c) one Friedel pair of a random candidate reference vector per position;
(d) twelve special positions (empty, beam only, exactly at the 0.15 Å⁻¹ beam radius, NaN and infinite peaks, one vector, exact
matrix / θ′ / β″ vectors, an exact Friedel pair, duplicates, far vectors) at an integer origin; (e) a malformed count; (f) zero width.
Generator: SplitMix64 written in `main.swift`; the dump carries every value it produced.

## Results worth knowing

- **-O vs -Onone.** `scans`, `matrixFits`, `runs` (every verdict, entry, score, count), `zoneAxisFits`, `claims`, `resolution`,
  `inPlaneSteps`, `nearest`, `detectorBases` are byte-identical. Four sections differ in the last bit:
  - `reflections`: Im F of 12 β″ and 14 triclinic reflections (|Δ| ≤ 2e-18) and |F|² of 2 triclinic ones (slice 2's amendment 1);
  - `libraries`: for every zone axis, only the in-plane step k = 13 (26 degrees) has vectors that differ (one ulp in q), 14 of 2520 entries
    in `main`, because at -O the app's `rotate` computes `cos(theta)` and `sin(theta)` together (a `__sincos_stret` call) and
    that pair differs from the separate libm calls at exactly this angle (measured in isolation: of the 180 steps, only k = 13). A
    kit `rotate` that wants the -O bits must also take the cos and sin of one angle in one function;
  - `euler`: `crystalToLabMatrix` of 2 of 270 cases (the angles themselves are identical), the same cause.
- **Dictionary seed.** `classify` ranks the phases from a `Dictionary` (`bestPerPhase.sorted`); an exact score tie between two
  phases would depend on the per-process hash seed. `run.sh` pins `SWIFT_DETERMINISTIC_HASHING=1`; a run without it (`--hash-check`)
  produced a byte-identical dump, so no tie occurs in these scans. The `.knownVariants` path sorts by entry index first and has no such dependence.
- **Input generation does not depend on the build**: the cos / sin in `main.swift` go through `@_optimize(none)` wrappers.
