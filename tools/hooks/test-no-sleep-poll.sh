#!/bin/bash
# Proof for no-sleep-poll.sh: feeds synthetic hook JSON, prints PASS/FAIL per case, exits non-zero on any failure.
cd "$(dirname "$0")" || exit 1
fail=0
check() { # expected-exit command
  local json; json=$(/usr/bin/python3 -I -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.argv[1]}}))' "$2")
  ./no-sleep-poll.sh <<<"$json" 2>/dev/null; local rc=$?
  if [[ $rc == "$1" ]]; then echo "PASS (exit $rc): $2"; else echo "FAIL (exit $rc, wanted $1): $2"; fail=1; fi
}
check 2 'sleep 60'
check 2 'sleep 30'
check 2 'sleep 45 && echo done'
check 2 'echo a; /bin/sleep 120'
check 2 'sleep 1m'
check 2 'sleep 0.5h'
check 2 'while true; do x; sleep 2; done'
check 2 'until grep -q x f; do sleep 1; done'
check 0 'echo "sleep 60"'   # quoted: passes (sleep must follow whitespace or an operator)
check 0 'sleep 1'
check 0 'sleep 5 && ls'
check 0 'sleep 29.9'
check 0 'git log'
check 0 'swift test'
check 0 'tools/test.sh --filter Sleeper'
# The machine-wide-kill guard (2026-10-15).
check 2 'pkill -f swift-test'
check 2 'pkill -f swift-build; echo x'
check 2 'pkill swiftpm-testing-helper'
check 2 'killall swift-frontend'
check 2 'echo a && sudo killall python3'
check 2 'kill $(pgrep -f xctest)'
check 2 'pgrep -f swift-test | xargs kill'
check 2 'kill -9 `pidof python3`'
check 2 'ps aux | grep swift | awk "{print \$2}" | xargs kill'
check 0 'pkill -f /private/tmp/claude-501/-Users-x/8f99abbe-a8c5-43a5-ad58-b674279022fc/scratchpad/dd'
check 0 'pkill -f 8f99abbe-a8c5-43a5-ad58-b674279022fc'
check 0 'pkill -P 12345'
check 0 'kill 12345'
check 0 'kill -TERM 12345 12346'
check 0 'pgrep -f swift-test'
check 0 'ps -axo pid,command | grep swift'
check 0 'git log --grep killall'
exit $fail
