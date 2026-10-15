#!/bin/bash
# PreToolUse hook (matcher "Bash"): refuses sleep and poll loops, and machine-wide kills. Reads the hook JSON on stdin,
# takes .tool_input.command. Exit 2 refuses (stderr goes to the agent). Matched as text, not parsed: a false positive is
# accepted (a quoted or mentioned word is refused when it looks like a command), a false negative is not.
# Wait on background-task notifications (Bash run_in_background) or the Monitor tool instead of sleeping.
#
# Normalisation first: shell comments are stripped (a `#` at the start or after whitespace, outside quotes, to the end of
# the line: a pin or a command hidden in a comment counts for nothing), `\<newline>` joined, `\x` and empty quote pairs
# inside a word removed.
# (1) sleep: refused unless every argument is a plain non-negative decimal literal under 30 s (suffix s m h d allowed),
#     quotes removed: `sleep infinity`, `sleep $((60))`, `sleep 0x40`, `sleep 1e2`, `sleep "45"`, `sleep 45.`, `sleep 30`
#     are refused; `sleep 5`, `sleep 0.5`, `sleep "2"` pass. Also as `env sleep`, `/bin/sleep`, `x.sleep(N)`, inside
#     `bash -c '...'`, `perl -e '...'` and similar (a quoted `sleep` counts when a shell, perl, python... is invoked; a plain
#     `echo "sleep 60"` passes), and `read -t N` with N >= 30. Any `while`/`until` with a sleep anywhere is refused.
# (2) killall and pkill, as a word anywhere (after env, command, nohup, exec, time, nice, sudo, xargs, eval, bash -c,
#     python os.system, in if/for/!, backslash-escaped, quoted), unless their own argument list (to the next ; & | newline
#     or the closing quote) names a /claude-501/ path or a session UUID. `-P <pid>` is no pin (`pkill -P 1`). Accepted
#     false positive: `git log --grep killall`.
# (3) kill: refused with pid 0 or a negative pid (`kill 0`, `kill -9 -1`, `kill -- -1`, `kill -KILL -1`), with no pid
#     (`kill -1`), fed by a command substitution or by xargs unless that substitution or pipeline names a /claude-501/
#     path or a session UUID, and `kill $var` next to a pgrep, pidof or ps pipe without a pin. `kill 12345` and
#     `kill -TERM 123 456` pass.
# (4) launchctl kill, bootout, remove, unload, stop, kickstart; `.kill(`, `killpg` in python, node, ruby, perl one-liners.
# Left allowed on purpose: a command built at run time (`P=pk; ${P}ill`, base64, `"pk""ill"` split inside a word), perl
# `select(undef,undef,undef,N)`, `timeout N`, `kill $pid` of a plain variable (no pgrep or ps beside it).
read -r -d '' PY <<'PYEOF'
import json, re, sys
try:
    cmd = json.load(sys.stdin).get("tool_input", {}).get("command", "")
except Exception:
    sys.exit(0)
if not isinstance(cmd, str):
    sys.exit(0)

SMSG = ("Refused by tools/hooks/no-sleep-poll.sh: %s. Do not sleep or poll. Run the job with Bash run_in_background "
        "and wait for its completion notification, or use the Monitor tool (an until-loop inside Monitor is its job).")
KMSG = ("Refused by tools/hooks/no-sleep-poll.sh: %s can stop other sessions' processes. Never pkill, killall or kill by "
        "pattern or by a pid list machine-wide; kill the literal pids you started (kill <pid>), or pin the pattern to a "
        "/claude-501/ path or a session UUID.")
PIN = re.compile(r"/claude-501/|[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}")


def refuse(text):
    print(text, file=sys.stderr)
    sys.exit(2)


def strip_comments(s):
    """Drop `#...` to the end of the line where `#` starts a word outside quotes. q: None, ', " or $' (ANSI-C)."""
    out, q, i, n = [], None, 0, len(s)
    while i < n:
        c = s[i]
        if q == "'":
            out.append(c)
            if c == "'":
                q = None
        elif q in ('"', "$'"):
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(s[i + 1])
                i += 1
            elif c == q[-1]:
                q = None
        else:
            if c == "\\" and i + 1 < n:
                out.append(c)
                out.append(s[i + 1])
                i += 1
            elif c == "$" and s[i + 1:i + 2] == "'":
                q = "$'"
                out.append("$'")
                i += 1
            elif c in "'\"":
                q = c
                out.append(c)
            elif c == "#" and (i == 0 or s[i - 1] in " \t\n"):
                while i < n and s[i] != "\n":
                    i += 1
                continue
            else:
                out.append(c)
        i += 1
    return "".join(out)


def quote_states(s):
    """st[i] = the quote state before character i (None, ', " or $')."""
    st, q, i, n = [], None, 0, len(s)
    while i < n:
        c = s[i]
        st.append(q)
        if q == "'":
            if c == "'":
                q = None
        elif q in ('"', "$'"):
            if c == "\\" and i + 1 < n:
                st.append(q)
                i += 1
            elif c == q[-1]:
                q = None
        else:
            if c == "\\" and i + 1 < n:
                st.append(q)
                i += 1
            elif c == "$" and s[i + 1:i + 2] == "'":
                st.append(q)
                q = "$'"
                i += 1
            elif c in "'\"":
                q = c
        i += 1
    return st


t = strip_comments(cmd)
t = t.replace("\\\n", "")
t = re.sub(r"\\([A-Za-z])", r"\1", t)
t = re.sub(r"(?<=\w)(?:\"\"|'')|(?:\"\"|'')(?=\w)", "", t)
ST = quote_states(t)
SEG = re.compile(r";|\n|&&|\|\||(?<![<>])&(?!>)")


def args_of(m, stops):
    """The argument text after word match m: to a character in stops, or to the quote that closes the one the word is in."""
    start, q = m.end(), ST[m.start()]
    if q and m.start() > 0 and t[m.start() - 1] == q[-1] and start < len(t) and t[start] == q[-1]:
        q, start = None, start + 1  # the word is quoted alone: "pkill" -f x
    end = start
    while end < len(t):
        c = t[end]
        if c in stops or (q and c == q[-1]):
            break
        end += 1
    return t[start:end]


def segment(pos):
    start = 0
    for m in SEG.finditer(t):
        if m.start() >= pos:
            return t[start:m.start()]
        start = m.end()
    return t[start:]


def word(name, after=r"(?![\w-])"):
    return re.finditer(r"(?<![\w-])" + name + after, t)


def substitutions(a):
    """The texts of the $( ) and backtick substitutions in a (an unterminated one runs to the end), and a with them blanked."""
    subs, rest, i = [], [], 0
    while i < len(a):
        if a.startswith("$(", i):
            depth, j = 1, i + 2
            while j < len(a) and depth:
                depth += (a[j] == "(") - (a[j] == ")")
                j += 1
            subs.append(a[i + 2:j - 1 if depth == 0 else j])
            rest.append("X")
            i = j
        elif a[i] == "`":
            j = a.find("`", i + 1)
            j = len(a) if j < 0 else j
            subs.append(a[i + 1:j])
            rest.append("X")
            i = j + 1
        else:
            rest.append(a[i])
            i += 1
    return subs, "".join(rest)


# (2) killall, pkill
for name in ("killall", "pkill"):
    for m in word(name):
        a = args_of(m, ";&|\n)`")
        if not PIN.search(a):
            refuse(KMSG % ("`%s%s`" % (name, a.rstrip())))

# (4) launchctl, in-language kills
if re.search(r"(?<![\w-])launchctl\b[^;&|\n]*?(?<![\w-])(?:kill|bootout|remove|unload|stop|kickstart)(?![\w-])", t):
    refuse(KMSG % "`launchctl kill/bootout/remove/unload/stop/kickstart`")
if re.search(r"killpg|\.kill\s*\(", t):
    refuse(KMSG % "an os.kill, killpg or process.kill call")

# (3) kill
for m in word("kill"):
    a = args_of(m, ";&\n")
    subs, rest = substitutions(a)
    for s in subs:
        if not PIN.search(s):
            refuse(KMSG % "`kill` fed by a command substitution")
    seg = segment(m.start())
    if re.search(r"(?<![\w-])xargs(?![\w-])", seg) and not PIN.search(seg):
        refuse(KMSG % "`xargs ... kill`")
    toks = [x.strip("()[]{},;'\"") for x in rest.split("|")[0].split()]
    toks = [x for x in toks if x and not re.match(r"^[0-9]*[<>]", x)]
    pids, sig, after_dd, listing, i = [], None, False, False, 0
    while i < len(toks):
        x = toks[i]
        if after_dd:
            pids.append(x)
        elif x == "--":
            after_dd = True
        elif x.startswith("-"):
            if x in ("-l", "-L"):
                listing = True
                break
            if sig is not None or pids:
                refuse(KMSG % ("`kill` with a negative pid (`%s`)" % x))
            sig = x
            if x in ("-s", "-n"):
                i += 1
        else:
            pids.append(x)
        i += 1
    if listing:
        continue
    for p in pids:
        if p.startswith("-") or re.fullmatch(r"\+?0+", p):  # pid 0 and negative pids (process groups, -1 = everything)
            refuse(KMSG % ("`kill` with pid %s" % p))
    if not pids and sig is not None and re.fullmatch(r"-[0-9]+", sig):
        refuse(KMSG % ("`kill %s` without a pid" % sig))
    if any("$" in p for p in pids) and not PIN.search(t):
        if re.search(r"(?<![\w-])(?:pgrep|pidof)(?![\w-])|(?<![\w-])ps\s+-?[a-zA-Z]*\b[^\n]*\|", t):
            refuse(KMSG % "a kill by pattern (pgrep, pidof or ps feeding kill)")

# (1) sleep
INVOKER = re.compile(r"(?<![\w-])(?:sh|bash|zsh|dash|ksh|eval|exec|xargs|perl|python[0-9.]*|ruby|node|osascript|expect|"
                     r"tclsh|awk)(?![\w-])|\s-c\b")
LIT = re.compile(r"([0-9]+(?:\.[0-9]+)?|\.[0-9]+)([smhd]?)")


def judge(tok):
    """None when tok is fine as a sleep or read -t argument, else the reason."""
    x = tok.replace('"', "").replace("'", "").strip("()[]{},;")
    if not x or re.match(r"^[0-9]*[<>]", x):
        return None
    m = LIT.fullmatch(x)
    if m:
        secs = float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]
        return None if secs < 30 else "sleep of %g s (limit: under 30 s)" % secs
    if x[0] in "0123456789.$`%" or x.lower().startswith("inf"):
        return "sleep of a non-literal or unbounded duration (%s)" % x
    return None


sleeps = []
for m in word("sleep", r"(?![\w.\-/])"):
    prev = t[m.start() - 1] if m.start() > 0 else ""
    if prev in ("'", '"'):
        alone = m.end() < len(t) and t[m.end()] == prev
        if not alone and not INVOKER.search(t):
            continue  # echo "sleep 60": a mention, not a command
    sleeps.append(m)
if sleeps and re.search(r"(?<![\w-])(?:while|until)(?![\w-])", t):
    refuse(SMSG % "a while/until loop containing sleep")
for m in sleeps:
    for tok in args_of(m, ";&|\n)`").split():
        why = judge(tok)
        if why:
            refuse(SMSG % why)
for m in word("read", r"(?![\w.\-/])"):
    toks = args_of(m, ";&|\n)`").split()
    for i, tok in enumerate(toks):
        mm = re.fullmatch(r"-[A-Za-z]*?t([0-9.]*)", tok)
        if mm:
            why = judge(mm.group(1) or (toks[i + 1] if i + 1 < len(toks) else ""))
            if why:
                refuse(SMSG % ("read -t: " + why))
sys.exit(0)
PYEOF
exec /usr/bin/python3 -I -c "$PY"
