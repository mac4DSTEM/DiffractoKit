#!/bin/zsh
# DiffractoKit slice 3b, oracle (registration docs/registration-slice-3.md section 3): the apps' own ellipse fit, fetched by
# `git show` from mac4DSTEM 7d6fdbd1 (read only: nothing is built, written or checked out there), compiled with main.swift
# into a scratch tool in $K/.build/apps-parity-3b (gitignored) twice, with swiftc -O (the apps' Release) and -Onone (what a
# debug `swift test` builds: the two differ in the last bits of some profile fits, as slice 2's amendment 1 found), run on the
# planted rings, and the dumps written to Tests/DiffractoGeometryTests/Fixtures/slice-3b-apps.json (-O) and
# slice-3b-apps-onone.json. The clean sweep rings' raw Float32 files go to .build/apps-parity-3b/sweep for tools/slice-3b/sweep.py.
#
# Then the KIT-CHECK step: the kit's own sources (EllipseCalibration.swift, DiffractionPattern.swift, and the two files that hold
# the model enum) are compiled with the same main.swift (plus kit-shim.swift) the same two ways, and each kit dump must be
# byte-identical to the apps' dump of the same optimisation level. This is the bit-identity result independent of SwiftPM; the
# Swift test EllipseAppsParityTests repeats it inside `swift test` (Debug against the -Onone dump, release against the -O one).
#
#   tools/slice-3b/apps-parity/run.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
HERE="$ROOT/tools/slice-3b/apps-parity"
APPS="$ROOT/../mac4DSTEM"
COMMIT="7d6fdbd1"
WORK="$ROOT/.build/apps-parity-3b"
FIX="$ROOT/Tests/DiffractoGeometryTests/Fixtures"
SRC="$ROOT/Sources/DiffractoGeometry"
FILES=(mac4DSTEM/Core/Analysis/EllipseCalibration.swift mac4DSTEM/Core/Data/DiffractionPattern.swift)

mkdir -p "$WORK/src" "$FIX"
full="$(git -C "$APPS" --no-optional-locks rev-parse "$COMMIT^{commit}")"
echo "apps-parity-3b: mac4DSTEM $full"
for f in "${FILES[@]}"; do
  git -C "$APPS" --no-optional-locks show "$COMMIT:$f" > "$WORK/src/${f:t}"
  echo "apps-parity-3b: $f  sha256 $(shasum -a 256 "$WORK/src/${f:t}" | cut -d' ' -f1)"
done
cp "$HERE/main.swift" "$WORK/src/main.swift"

# the apps' types use `package` access; both files import Foundation only. The two optimisation levels build and run side by side.
swiftc -O -swift-version 6 -package-name x -o "$WORK/apps-parity-3b-O" \
  "$WORK/src/EllipseCalibration.swift" "$WORK/src/DiffractionPattern.swift" "$WORK/src/main.swift"
swiftc -Onone -swift-version 6 -package-name x -o "$WORK/apps-parity-3b-Onone" \
  "$WORK/src/EllipseCalibration.swift" "$WORK/src/DiffractionPattern.swift" "$WORK/src/main.swift"
# the kit, same main.swift, same two levels (whole-module, as SwiftPM builds)
KIT=("$SRC/EllipseCalibration.swift" "$SRC/DiffractionPattern.swift" "$SRC/RingCalibration.swift" "$SRC/ElectronWavelength.swift"
     "$HERE/kit-shim.swift" "$HERE/main.swift")
swiftc -O -wmo -swift-version 6 -package-name x -o "$WORK/kit-parity-3b-O" "${KIT[@]}"
swiftc -Onone -wmo -swift-version 6 -package-name x -o "$WORK/kit-parity-3b-Onone" "${KIT[@]}"

rm -rf "$WORK/sweep"
"$WORK/apps-parity-3b-O" "$WORK/dump-O.json" "$WORK/sweep" 2> "$WORK/apps-O.log" &
"$WORK/apps-parity-3b-Onone" "$WORK/dump-Onone.json" 2> "$WORK/apps-Onone.log" &
"$WORK/kit-parity-3b-O" "$WORK/kit-O.json" 2> "$WORK/kit-O.log" &
"$WORK/kit-parity-3b-Onone" "$WORK/kit-Onone.json" 2> "$WORK/kit-Onone.log" &
wait

cp "$WORK/dump-O.json" "$FIX/slice-3b-apps.json"
cp "$WORK/dump-Onone.json" "$FIX/slice-3b-apps-onone.json"
ls -l "$FIX"/slice-3b-apps.json "$FIX"/slice-3b-apps-onone.json
for level in O Onone; do
  if cmp -s "$WORK/dump-$level.json" "$WORK/kit-$level.json"; then
    echo "apps-parity-3b: kit -$level equals the apps' -$level in every output bit"
  else
    echo "apps-parity-3b: kit -$level DIFFERS from the apps' -$level" ; exit 1
  fi
done
if cmp -s "$WORK/dump-O.json" "$WORK/dump-Onone.json"; then
  echo "apps-parity-3b: the apps' -O and -Onone dumps are byte-identical"
else
  echo "apps-parity-3b: the apps' -O and -Onone dumps DIFFER (optimisation level moves the last bits; each is the reference for its own build)"
fi
