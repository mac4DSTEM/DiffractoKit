#!/bin/zsh
# Builds the benchmark package and prints the kit's JSON lines. Run under ../../../tools/memcap.sh, with the machine otherwise
# quiet (the numbers are this Mac at that moment). DIFFRACTO_BENCH_LARGE=1 for the nightly sizes.
set -euo pipefail
cd ${0:A:h}
swift build -c release 2>&1 | tail -3 >&2
print -u2 "load: $(uptime | sed 's/.*load averages*: //')"
.build/release/diffracto-bench inputs
