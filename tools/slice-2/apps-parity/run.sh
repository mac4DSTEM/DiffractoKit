#!/bin/zsh
# DiffractoKit slice 2, oracle 3 (registration docs/registration-slice-2.md section 2 and 4.4): the apps' own reflection code,
# fetched by `git show` from mac4DSTEM 7d6fdbd1 (read only: nothing is built, written or checked out there), compiled with
# main.swift into a scratch tool in $K/.build/apps-parity (gitignored), run, and its dump written to
# Tests/DiffractoScatteringTests/Fixtures/slice-2-apps.json (compiled -O, the apps' Release). The same sources compiled
# -Onone are dumped too and compared (the dump stays in .build/apps-parity; registration amendment 1 records that they
# differ only in Im F below 1e-15 of beta'').
#
#   tools/slice-2/apps-parity/run.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
HERE="$ROOT/tools/slice-2/apps-parity"
APPS="$ROOT/../mac4DSTEM"
COMMIT="7d6fdbd1"
WORK="$ROOT/.build/apps-parity"
FIX="$ROOT/Tests/DiffractoScatteringTests/Fixtures"
FILES=(mac4DSTEM/Core/Crystal/Crystal.swift mac4DSTEM/Core/Crystal/ScatteringFactors.swift)

mkdir -p "$WORK/src" "$FIX"
full="$(git -C "$APPS" --no-optional-locks rev-parse "$COMMIT^{commit}")"
echo "apps-parity: mac4DSTEM $full"
for f in "${FILES[@]}"; do
  git -C "$APPS" --no-optional-locks show "$COMMIT:$f" > "$WORK/src/${f:t}"
  echo "apps-parity: $f  sha256 $(shasum -a 256 "$WORK/src/${f:t}" | cut -d' ' -f1)"
done
cp "$HERE/main.swift" "$WORK/src/main.swift"

# the apps' types use `package` access; both files import Foundation and simd only
for mode in O Onone; do
  swiftc -$mode -swift-version 6 -package-name appsparity -o "$WORK/apps-parity-$mode" \
    "$WORK/src/Crystal.swift" "$WORK/src/ScatteringFactors.swift" "$WORK/src/main.swift"
  "$WORK/apps-parity-$mode" "$WORK/dump-$mode.json"
done

cp "$WORK/dump-O.json" "$FIX/slice-2-apps.json"
if cmp -s "$WORK/dump-O.json" "$WORK/dump-Onone.json"; then
  echo "apps-parity: -O and -Onone dumps are byte-identical (every bit pattern agrees)"
else
  echo "apps-parity: -O and -Onone DIFFER (see registration amendment 1; the -O dump is the reference)"
fi
ls -l "$FIX"/slice-2-apps.json
