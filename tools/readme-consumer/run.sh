#!/bin/zsh
# DiffractoKit readme-consumer: prove that the README's "First ten minutes" is what it says it is. Exit 0 only if
#   1. Sources/readme-consumer/main.swift equals, byte for byte, the README's program (the first ```swift block after the
#      heading "First ten minutes") and Examples/diffracto-example/main.swift;
#   2. the separate package builds in DEBUG against the public API with no warning, and no test-only target is built;
#   3. the program's output equals the README's stated output (the first plain ``` block after "Output of the program").
#
#   run.sh
#
#   READMEC_JOBS=N   passes `-j N` to swift build.
set -euo pipefail
here=${0:A:h}
root=${here:h:h}
die() { print -u2 "readme-consumer: $*"; exit 1 }
mkdir -p $here/.build
code=$here/.build/readme-code.swift
expected=$here/.build/readme-output.txt

# The README's program and output, cut out with awk (fences are exact lines).
awk '/^## First ten minutes/ {on=1} grab && /^```$/ {grab=0; done=1; next} on && /^```swift$/ && !done {grab=1; next} grab {print}' \
  $root/README.md >$code
awk '/^Output of the program/ {on=1} grab && /^```$/ {grab=0; done=1; next} on && /^```$/ && !done {grab=1; next} grab {print}' \
  $root/README.md >$expected
[[ -s $code ]] || die "no swift block found under 'First ten minutes' in README.md"
[[ -s $expected ]] || die "no output block found under 'Output of the program' in README.md"
cmp -s $code $here/Sources/readme-consumer/main.swift || die "main.swift differs from the README's program (diff $code $here/Sources/readme-consumer/main.swift)"
cmp -s $code $root/Examples/diffracto-example/main.swift || die "Examples/diffracto-example/main.swift differs from the README's program"

jobs=()
[[ -n ${READMEC_JOBS:-} ]] && jobs=(-j $READMEC_JOBS)
cd $here
log=.build/readme-consumer-build.log
# --disable-dependency-cache: SwiftPM's shared repository cache can print a 'skipping cache' warning on a CI runner,
# which the warning count would take for the kit's own (DiffractoKit 0.4.1's CI).
if ! swift build --disable-dependency-cache -c debug $jobs >$log 2>&1; then
  grep -E "error" $log | head -20 >&2
  die "BUILD FAILED (full log: $here/$log)"
fi
warnings=$(grep -c "warning:" $log || true)
(( warnings == 0 )) || { grep "warning:" $log | head >&2; die "$warnings warning(s) in the build" }
if grep -qE "DiffractoCoreTests|DiffractoScatteringTests|DiffractoGeometryTests|DiffractoOrientationTests|DiffractoKitParityTests|DiffractoKitTestSupport" $log; then
  die "the build log names a test-only target"
fi
bin=$(swift build --disable-dependency-cache -c debug --show-bin-path)
actual=.build/readme-actual.txt
set +e
$bin/readme-consumer >$actual 2>&1
rc=$?
set -e
(( rc == 0 )) || { cat $actual; die "the program FAILED (exit $rc)" }
if ! diff $expected $actual >.build/readme-output.diff; then
  cat .build/readme-output.diff >&2
  die "the output differs from the README's stated output (diff above: < README, > actual)"
fi
head=$(git -C $root --no-optional-locks rev-parse --short HEAD 2>/dev/null || print unknown)
print "readme-consumer: passed — the README's program, byte for byte, built in Debug with no warnings and no test-only target; its output equals the README's ($(wc -l <$actual | tr -d ' ') lines) (DiffractoKit $head)"
