#!/bin/zsh
# DiffractoKit consumer check: build tools/consumer-check (a separate package that depends on DiffractoKit by path) in
# DEBUG, as a consumer's Debug build does with a dependency that may carry no unsafe flags, and run it. Exit 0 only if the
# build is warning-free, no test-only target was built, and every check passed.
#
#   run.sh [--controls]
#
#   (no flag)   build, run the checks, print "consumer-check: all passed — N checks ..." and exit 0.
#   --controls  then the negative controls (controls.py, python3 standard library only; SpectroKit's pattern): for each
#               public symbol the consumer uses, an edited COPY of the repository (rsync into
#               <root>/.build/consumer-check/mutants/, never the working tree) has that one declaration reverted from
#               `public` to `package`; the consumer is built against the copy (one warm build directory, an incremental
#               build per control) and must FAIL with the compiler's protection-level diagnostic naming that symbol.
#               Prints "consumer-check: negative controls: N of M caught" and exits non-zero unless all are caught.
#
#   CONSUMER_CHECK_JOBS=N   passes `-j N` to every swift build (the plain build and each control).
#   CONSUMER_CHECK_PYTHON   the python3 for --controls.
set -euo pipefail
here=${0:A:h}
root=${here:h:h}
die() { print -u2 "consumer-check: $*"; exit 1 }

controls=0
for arg in "$@"; do
  case $arg in
    --controls) controls=1 ;;
    -h|--help) sed -n '2,/^set /p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown argument $arg (usage: run.sh [--controls])" ;;
  esac
done
jobs=()
[[ -n ${CONSUMER_CHECK_JOBS:-} ]] && jobs=(-j $CONSUMER_CHECK_JOBS)
cd $here
mkdir -p .build
log=.build/consumer-check-build.log
if ! swift build -c debug $jobs >$log 2>&1; then
  grep -E "error" $log | head -20 >&2
  die "BUILD FAILED (full log: $here/$log)"
fi
warnings=$(grep -c "warning:" $log || true)
(( warnings == 0 )) || { grep "warning:" $log | head >&2; die "$warnings warning(s) in the build" }
if grep -qE "DiffractoCoreTests|DiffractoScatteringTests|DiffractoGeometryTests|DiffractoOrientationTests" $log; then
  die "the build log names a test-only target"
fi
bin=$(swift build -c debug --show-bin-path)
out=.build/consumer-check-run.log
set +e
$bin/consumer-check >$out 2>&1
rc=$?
set -e
cat $out
(( rc == 0 )) || die "FAILED (exit $rc)"
summary=$(grep "^consumer-check summary:" $out)
checks=${${summary#consumer-check summary: }%% checks*}
head=$(git -C $root --no-optional-locks rev-parse --short HEAD 2>/dev/null || print unknown)
print "consumer-check: all passed — $checks checks through DiffractoKit's public API only, Debug build, no warnings, no test-only target built (DiffractoKit $head)"
(( controls )) || exit 0

# Negative controls. The repository is never edited: controls.py works in an rsync'd copy under .build.
PY=${CONSUMER_CHECK_PYTHON:-python3}
command -v $PY >/dev/null || die "no $PY for --controls (set CONSUMER_CHECK_PYTHON)"
exec $PY -I -B $here/controls.py $root
