#!/bin/bash
# Proof for no-sleep-poll.sh: feeds synthetic hook JSON, prints PASS/FAIL per case, exits non-zero on any failure.
# The commands are only TEXT piped to the hook; nothing here is executed.
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
check 0 'echo "sleep 60"'   # quoted mention, no shell invoked: passes
check 0 'sleep 1'
check 0 'sleep 5 && ls'
check 0 'sleep 29.9'
check 0 'git log'
check 0 'swift test'
check 0 'tools/test.sh --filter Sleeper'
check 0 'echo done'
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
check 0 'kill 12345'
check 0 'kill -TERM 12345 12346'
check 0 'kill -TERM 123 456'
check 0 'kill -0 12345'
check 0 'kill 12345 2>/dev/null || true'
check 0 'pgrep -f swift-test'
check 0 'ps -axo pid,command | grep swift'
check 0 'ls # pkill -f swift'
check 0 'echo hi # sleep 100'
# 2026-10-15 hardening (refuter 0.2.0): killall/pkill not at a separator.
check 2 'env pkill -f swift'
check 2 'command pkill -f swift'
check 2 'nohup pkill -f swift'
check 2 'exec pkill -f swift'
check 2 'time pkill -f swift'
check 2 'nice pkill -f swift'
check 2 'FOO=1 pkill -f swift'
check 2 "bash -c 'pkill -f swift'"
check 2 'sh -c "pkill -f swift"'
check 2 'bash -c "killall Xcode"'
check 2 'eval pkill -f swift'
check 2 'if true; then pkill -f swift; fi'
check 2 'for i in 1; do pkill -f swift; done'
check 2 '! pkill -f swift'
check 2 '"pkill" -f swift'
check 2 "'pkill' -f swift"
check 2 '\pkill -f swift'
check 2 '\killall Xcode'
check 2 'sudo -n pkill -f swift'
check 2 'sudo -u root pkill -f swift'
check 2 'xargs pkill -f < /tmp/names'
check 2 'echo swift | xargs pkill -f'
check 2 'echo swift | xargs -I{} pkill -f {}'
check 2 'python3 -c "import os; os.system('"'"'pkill -f swift'"'"')"'
check 2 'python3 -c "import os; os.system(\"killall Xcode\")"'
check 2 $'echo a\npkill -f swift'
check 2 'pkill \
  -f swift'
check 2 'git log --grep killall'   # accepted false positive: a word after --grep is refused like a command
# Pinning: -P is no pin; a pin in a comment or elsewhere on the line is no pin.
check 2 'pkill -P 12345'
check 2 'pkill -P 1 -f swift'
check 2 'pkill -P 1'
check 2 'pkill -9 -f . -P1'
check 2 'pkill -f swift #/claude-501/'
check 2 'pkill -f . # 8f99abbe-a8c5-43a5-ad58-b674279022fc'
check 2 'pgrep -f swift | xargs kill #/claude-501/'
check 2 'kill $(pgrep -f swift) ; echo /claude-501/'
check 2 'pkill -f swift; echo /private/tmp/claude-501/x'
check 2 "bash -c 'pkill -f swift'; echo 8f99abbe-a8c5-43a5-ad58-b674279022fc"
check 0 "bash -c 'pkill -f /private/tmp/claude-501/x/dd'"
check 0 'pgrep -f /private/tmp/claude-501/x/dd | xargs kill'
check 0 'kill $(pgrep -f /private/tmp/claude-501/x/dd)'
# kill without pgrep or ps: pid 0, negative pids, other sources of pids.
check 2 'kill -- -1'
check 2 'kill -9 -1'
check 2 'kill -KILL -1'
check 2 'kill 0'
check 2 'kill -9 0'
check 2 'kill -s KILL -1'
check 2 'kill -1'
check 2 'kill -9 -12345'
check 2 'kill -9 $(ps -ax -o pid=)'
check 2 'kill $(lsof -t -c swift)'
check 2 'lsof -t -c swift | xargs kill'
check 2 'xargs kill < /tmp/pids'
check 2 'kill $(cat /tmp/pids)'
check 2 "kill \$(launchctl list | awk '{print \$1}')"
check 2 "top -l1 | awk '{print \$1}' | xargs kill"
check 2 'PIDS=$(pgrep -f swift); kill $PIDS'
check 2 'python3 -c "import os,signal; os.kill(-1, signal.SIGKILL)"'
check 2 'python3 -c "import os; os.killpg(0, 9)"'
check 2 'perl -e "kill 9, -1"'
check 2 'launchctl kill SIGKILL system/com.apple.foo'
check 2 'launchctl bootout gui/501'
check 0 'launchctl list'
check 0 'kill $!'
# sleep guard: durations that are not plain literals under 30 s, and sleeps inside other interpreters.
check 2 'sleep 1e2'
check 2 'sleep infinity'
check 2 'sleep "45"'
check 2 "sleep '45'"
check 2 'sleep $((60))'
check 2 'sleep 45.'
check 2 'sleep 0x40'
check 2 'sleep $SECS'
check 2 "bash -c 'sleep 100'"
check 2 'sh -c "sleep 100"'
check 2 "perl -e 'sleep 100'"
check 2 'python3 -c "import time; time.sleep(100)"'
check 2 'read -t 100 x'
check 2 'read -rt 100 x'
check 2 'read -t100 x'
check 2 'env sleep 100'
check 2 'sleep 100 # comment'
check 2 "bash -c 'while true; do sleep 1; done'"
check 0 'sleep 5'
check 0 'sleep 0.5'
check 0 'sleep "2"'
check 0 'sleep 5 2>/dev/null'
check 0 "bash -c 'sleep 2'"
check 0 'read -t 5 x'
check 0 'read -r line'
check 0 'grep -n sleep tools/test.sh'
check 0 'cat tools/hooks/no-sleep-poll.sh'
exit $fail
