#!/bin/zsh
# DiffractoKit slice 4, oracle (registration docs/registration-slice-4.md section 5 "The apps' compiled code" and gate G1):
# the apps' own vector-matched phase mapping, fetched by `git show` from mac4DSTEM 7d6fdbd1 (read only: nothing is built,
# written or checked out there), compiled with main.swift into a scratch tool in $K/.build/apps-parity-4 (gitignored) twice,
# with swiftc -O (the apps' Release) and -Onone (what a debug `swift test` builds), run, and the -O dump written to
# Tests/DiffractoOrientationTests/Fixtures/slice-4-apps.json (the -Onone dump stays in .build/apps-parity-4 and is compared).
#
# Source files compiled verbatim: Crystal.swift, ScatteringFactors.swift, PhaseReferenceLibrary.swift, PhaseVectorMatching.swift,
# and the standalone AnalysisCancellationToken.swift. Pieces cut out of larger app files BY LINE RANGE (verbatim, no edit,
# the first line of each range is checked): BraggPeak and BraggVectors (the struct heads, without `calibrated`) from
# DiskDetection.swift, which imports Accelerate and Metal and needs the whole analysis layer; EulerAngles and ACOMOrientation
# (lines 1-2 and 12-107) from OrientationResult.swift, whose other types are slice 5. No hand-written shim holds numeric logic.
#
#   tools/slice-4/apps-parity/run.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
HERE="$ROOT/tools/slice-4/apps-parity"
APPS="$ROOT/../mac4DSTEM"
COMMIT="7d6fdbd1"
WORK="$ROOT/.build/apps-parity-4"
FIX="$ROOT/Tests/DiffractoOrientationTests/Fixtures"
WHOLE=(mac4DSTEM/Core/Crystal/Crystal.swift mac4DSTEM/Core/Crystal/ScatteringFactors.swift
       mac4DSTEM/Core/Crystal/PhaseReferenceLibrary.swift mac4DSTEM/Core/Crystal/PhaseVectorMatching.swift
       mac4DSTEM/Core/Compute/AnalysisCancellationToken.swift)

mkdir -p "$WORK/src" "$FIX"
full="$(git -C "$APPS" --no-optional-locks rev-parse "$COMMIT^{commit}")"
echo "apps-parity-4: mac4DSTEM $full"
show() { git -C "$APPS" --no-optional-locks show "$COMMIT:$1"; }
for f in "${WHOLE[@]}"; do
  show "$f" > "$WORK/src/${f:t}"
  echo "apps-parity-4: $f  sha256 $(shasum -a 256 "$WORK/src/${f:t}" | cut -d' ' -f1)"
done

# line-range cuts, each guarded by the text its first line must carry
cut_range() {   # file from to expected-first-line-substring out
  local text; text="$(show "$1" | sed -n "$2,$3p")"
  [[ "$(print -r -- "$text" | head -1)" == *"$4"* ]] || { echo "apps-parity-4: $1 lines $2-$3 no longer start with '$4'"; exit 1; }
  print -r -- "$text" >> "$5"
}
rm -f "$WORK/src/CutOut.swift" "$WORK/src/OrientationCut.swift"
{ echo "import Foundation"; echo; } > "$WORK/src/CutOut.swift"
cut_range mac4DSTEM/Core/Analysis/DiskDetection.swift 34 45 "struct BraggPeak" "$WORK/src/CutOut.swift"
echo >> "$WORK/src/CutOut.swift"
cut_range mac4DSTEM/Core/Analysis/DiskDetection.swift 689 705 "struct BraggVectors" "$WORK/src/CutOut.swift"
echo "}" >> "$WORK/src/CutOut.swift"
cut_range mac4DSTEM/Core/Analysis/OrientationResult.swift 1 2 "import Foundation" "$WORK/src/OrientationCut.swift"
cut_range mac4DSTEM/Core/Analysis/OrientationResult.swift 12 107 "py4DSTEM/orix-compatible" "$WORK/src/OrientationCut.swift"
echo "apps-parity-4: cut-outs sha256 $(shasum -a 256 "$WORK/src/CutOut.swift" | cut -d' ' -f1) $(shasum -a 256 "$WORK/src/OrientationCut.swift" | cut -d' ' -f1)"
cp "$HERE/main.swift" "$WORK/src/main.swift"

SRCS=("$WORK/src/Crystal.swift" "$WORK/src/ScatteringFactors.swift" "$WORK/src/PhaseReferenceLibrary.swift"
      "$WORK/src/PhaseVectorMatching.swift" "$WORK/src/AnalysisCancellationToken.swift" "$WORK/src/CutOut.swift"
      "$WORK/src/OrientationCut.swift" "$WORK/src/main.swift")
# the apps' types use `package` access. One swiftc at a time; whole-module, as SwiftPM and the app build.
for mode in O Onone; do
  swiftc -$mode -wmo -swift-version 6 -package-name appsparity -o "$WORK/apps-parity-4-$mode" "${SRCS[@]}"
done

# SWIFT_DETERMINISTIC_HASHING pins the per-process Dictionary seed: the app's classify() ranks phases from a Dictionary, and
# an exact score tie between two phases would otherwise depend on the seed (see the README section). The -O dump is also
# made once WITHOUT it and compared, to record whether any output depends on the seed.
export SWIFT_DETERMINISTIC_HASHING=1
SECONDS=0
"$WORK/apps-parity-4-O" "$WORK/dump-O.json" > "$WORK/run-O.log" 2>&1
"$WORK/apps-parity-4-Onone" "$WORK/dump-Onone.json" > "$WORK/run-Onone.log" 2>&1
echo "apps-parity-4: both dumps written in ${SECONDS} s"
cat "$WORK/run-O.log"

cp "$WORK/dump-O.json" "$FIX/slice-4-apps.json"
ls -l "$FIX"/slice-4-apps.json
if cmp -s "$WORK/dump-O.json" "$WORK/dump-Onone.json"; then
  echo "apps-parity-4: -O and -Onone dumps are byte-identical (every bit pattern agrees)"
else
  echo "apps-parity-4: -O and -Onone DIFFER (the -O dump is the reference); first differing bytes:"
  cmp "$WORK/dump-O.json" "$WORK/dump-Onone.json" | head -1 || true
fi
if [[ "${1:-}" == "--hash-check" ]]; then
  unset SWIFT_DETERMINISTIC_HASHING
  "$WORK/apps-parity-4-O" "$WORK/dump-O-seed2.json" > /dev/null 2>&1
  if cmp -s "$WORK/dump-O.json" "$WORK/dump-O-seed2.json"; then
    echo "apps-parity-4: a run with a random Dictionary seed is byte-identical to the pinned run"
  else
    echo "apps-parity-4: a run with a random Dictionary seed DIFFERS from the pinned run"; cmp "$WORK/dump-O.json" "$WORK/dump-O-seed2.json" | head -1 || true
  fi
fi
