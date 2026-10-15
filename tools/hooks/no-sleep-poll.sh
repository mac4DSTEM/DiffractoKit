#!/bin/bash
# PreToolUse hook (matcher "Bash"): refuses sleep and poll loops. Reads the hook JSON on stdin, takes .tool_input.command.
# Refused with exit 2 (stderr goes to the agent): (1) any `sleep` of 30 s or more (`sleep 45`, `sleep 1m`, `sleep 0.5h`,
# `/bin/sleep 120`, also inside && / ; chains); (2) a `while` or `until` loop with `sleep` anywhere in the command, any length.
# Everything else passes (`sleep 1`, `sleep 5`, `git log`). Matched as text, not parsed: a sleep right after a
# quote passes, a quoted sleep after a space is refused (accepted false positive).
# Wait on background-task notifications (Bash run_in_background) or the Monitor tool instead.
# Also refused (2026-10-15: a lane ran machine-wide `pkill -f swift-test` and stopped other sessions' jobs; the guard follows
# macTEM 08ebfa0's): (3) `killall` and `pkill` (with or without -f) unless pinned to this machine's session scratch, that is
# their arguments name a /claude-501/ path, a session UUID, or `-P <pid>`; (4) a kill by pattern: `kill` fed by pgrep, pidof
# or ps (`kill $(pgrep ...)`, `pgrep ... | xargs kill`) under the same pinning rule. Kill by pid only: `kill 12345`.
exec /usr/bin/python3 -I -c '
import json, re, sys
try:
    cmd = json.load(sys.stdin).get("tool_input", {}).get("command", "")
except Exception:
    sys.exit(0)
if not isinstance(cmd, str):
    sys.exit(0)
msg = ("Refused by tools/hooks/no-sleep-poll.sh: %s. Do not sleep or poll. Run the job with Bash run_in_background "
       "and wait for its completion notification, or use the Monitor tool (an until-loop inside Monitor is its job).")
PINNED = r"(/claude-501/|[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}|-P\s*[0-9]+)"
kmsg = ("Refused by tools/hooks/no-sleep-poll.sh: %s can stop other sessions\x27 processes. Never pkill or killall machine-wide; "
        "kill the pids you started (kill <pid>), or pin the pattern to a /claude-501/ path, a session UUID or -P <pid>.")
flat = cmd.replace("\n", ";")
for m in re.finditer(r"(?:^|[;&|(`{])\s*(?:sudo\s+)?(?:/\S*/)?(killall|pkill)\b([^;&|)`]*)", flat):
    if not re.search(PINNED, m.group(2)):
        print(kmsg % ("`%s%s`" % (m.group(1), m.group(2).rstrip())), file=sys.stderr); sys.exit(2)
if re.search(r"(?:^|[\s;&|(`{])(?:/\S*/)?(?:kill|xargs(?:\s+-\S+)*\s+kill)\b", flat) and re.search(r"\b(?:pgrep|pidof)\b|\bps\s+-?[a-zA-Z]*\b.*\|", flat):
    if not re.search(PINNED, flat):
        print(kmsg % "a kill by pattern (pgrep, pidof or ps feeding kill)", file=sys.stderr); sys.exit(2)
sleeps = list(re.finditer(r"(?:^|[\s;&|()`{])(?:/\S*/)?sleep\s+([0-9]*\.?[0-9]+)([smhd]?)(?![\w.])", cmd))
if re.search(r"\b(?:while|until)\b", cmd) and re.search(r"(?:^|[\s;&|()`{])(?:/\S*/)?sleep\b", cmd):
    print(msg % "a while/until loop containing sleep", file=sys.stderr); sys.exit(2)
for m in sleeps:
    secs = float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]
    if secs >= 30:
        print(msg % ("sleep of %g s (limit: under 30 s)" % secs), file=sys.stderr); sys.exit(2)
sys.exit(0)
'
