#!/bin/zsh
# DiffractoKit ACOM benchmark: builds the separate package in release and prints the timings of `ACOMMatcher.matchAll`
# (library generation, one worker, 2 / 4 / 8 workers, the cost of the cancellation poll). Heavy CPU: run it with the machine
# otherwise quiet and note the load; the numbers are a measurement of this Mac at that moment, not a pin.
#
#   run.sh [width] [height]     scan size, default 64 x 64
set -euo pipefail
here=${0:A:h}
cd $here
swift build -c release 2>&1 | tail -3
print "load: $(uptime | sed 's/.*load averages*: //')"
.build/release/acom-benchmark "$@"
