#!/bin/zsh
# PROVENANCE: macTEM before v1.0.0 (history squashed) tools/lib/fetch-py4dstem-scattering.sh adapted
# Adapted: the same py4DSTEM lock commit and the same sha256, the destination is this kit's gitignored References/ (the
# script sits two levels deeper than macTEM's, in tools/pins/structure-factors/, and the root is three levels up).
# END PROVENANCE
# The py4DSTEM electron-scattering-factor table lock - fetched, not vendored.
#
# tools/pins/structure-factors/pins.py builds its reference f_e from THIS table, not from any Swift table: the pins are independent of the
# port only if the table comes from elsewhere. The file is py4DSTEM/process/utils/scattering_factors.txt at the lock commit
# (103 rows, row Z-1 = a0 b0 a1 b1 ... a4 b4); py4DSTEM's single_atom_scatter.py reads columns 0:10:2 as a and 1:10:2 as b.
# Fetched into the gitignored References/. A present file is used only when its sha256 is the pin; any other is REFUSED.
set -euo pipefail

PY4DSTEM_LOCK_COMMIT="f050d207851f8accd3aa550365baaf1f38230428"   # dev, 2026-03-26, version 0.14.19 (mac4DSTEM's lock)
TABLE_SHA256="be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4"
URL="https://raw.githubusercontent.com/py4dstem/py4DSTEM/$PY4DSTEM_LOCK_COMMIT/py4DSTEM/process/utils/scattering_factors.txt"
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DEST="$ROOT/References/py4dstem-scattering/scattering_factors.txt"

have() { shasum -a 256 "$DEST" 2>/dev/null | cut -d' ' -f1; }

if [[ -f "$DEST" && "$(have)" == "$TABLE_SHA256" ]]; then
  echo "fetch-table: lock present (sha256 $TABLE_SHA256)"
  exit 0
fi
mkdir -p "$(dirname "$DEST")"
tmp="$(mktemp "$DEST.XXXXXX")"
trap 'rm -f "$tmp"' EXIT
echo "fetch-table: fetching scattering_factors.txt at py4DSTEM $PY4DSTEM_LOCK_COMMIT" >&2
curl -fsSL --max-time 60 "$URL" -o "$tmp"
got="$(shasum -a 256 "$tmp" | cut -d' ' -f1)"
if [[ "$got" != "$TABLE_SHA256" ]]; then
  echo "fetch-table: REFUSED - the fetched table has sha256 $got, the lock is $TABLE_SHA256" >&2
  exit 1
fi
mv "$tmp" "$DEST"
trap - EXIT
echo "fetch-table: ready"
