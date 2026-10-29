#!/usr/bin/env python3
"""Pins for RingIndexing.

Port of macTEM 9482831 tools/indexing-parity/reference.py (read by `git show`) onto the KIT's own ring lists: the families
after the |F| cut that slice 2 pinned (Tests/DiffractoScatteringTests/Fixtures/numpy-structure-factor-pins.json, written by
tools/pins/structure-factors/pins.py), at k_max 1.2. The decision is written here from its definition in DiffractoGeometry/RingIndexing.swift's header (NumPy, plus a plain-Python
scalar re-implementation that must agree on every plant), NOT from any Swift:

  z = |q - g| / sigma, delta = g / q - 1, candidates ordered by (z, g); within = z <= k (not z^2 <= k);
  ambiguous when >= 2 within (payload: all of them, in (z, g) order) or exactly one within whose margin
  (z_next^2 - z^2) / 2 to the 2nd candidate is < ln 19 (payload [best, next]); r = sqrt(k^2 + 2 ln 19);
  beyondKMax when q + R > k_max, R = max(r sigma, |q - g| of best and next if one is within, of the nearest if none);
  noneWithinTolerance (nearest, or none for an empty list); else accepted (best, next, margin).
  Verdict order: invalidInput, ambiguous, beyondKMax, noneWithinTolerance, accepted.

Plants. T1: every family of every list, u in {1e-4, 1e-3, 3e-3, 1e-2, 3e-2}, s in {-2.9, -2, -1, -0.5, 0, 0.5, 1, 2, 2.9},
q = g (1 + s u), sigma = u q, k = 3 (closed forms P1: z = |s| / (1 + s u), delta = -s u / (1 + s u)). T4: the refusal plants
macTEM registered, on the kit's lists (invalid numbers, beyondKMax, forbidden positions, ambiguous pairs, the margin plant,
the one-family and the empty list). Lists: Au, Al, Si, Mg, ideal hcp, Weak-cut (macTEM's six, so the T1 counts equal its
registered prediction) and beta'' (macTEM's accepted counts). Non-finite numbers are written as the strings
"NaN", "Infinity", "-Infinity".

Writes Tests/DiffractoGeometryTests/Fixtures/numpy-ring-indexing-pins.json BEFORE the Swift test is run. Deterministic (no clock,
no random numbers). Run with ../References/python/pymatgen/bin/python -I -B tools/pins/ring-indexing/pins.py
Exit 0: every self-check holds. 1: a self-check failed (the file is still written, with the failures listed).
"""
import json
import math
import os
import sys

import numpy as np

assert np.__version__ == "2.5.4", f"numpy {np.__version__}: the lock is 2.5.4"

HERE = os.path.dirname(os.path.abspath(__file__))
KIT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
PINS2 = os.path.join(KIT, "Tests", "DiffractoScatteringTests", "Fixtures", "numpy-structure-factor-pins.json")
OUT = os.path.join(KIT, "Tests", "DiffractoGeometryTests", "Fixtures", "numpy-ring-indexing-pins.json")

K = 3.0
LN19 = math.log(19.0)
KMAX = 1.2
U1 = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2)
S1 = (-2.9, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 2.9)
# macTEM's registered T1 counts per u over its six crystals (accepted as own, accepted as another family, ambiguous,
# beyondKMax, none): Au, Al, Si, Mg, Mg-ideal, Weak-cut; and beta'' ('accepted in all' per u).
T1_PREDICT = {1e-4: (711, 0, 0, 36, 0), 1e-3: (691, 0, 20, 36, 0), 3e-3: (647, 0, 64, 36, 0),
              1e-2: (408, 0, 309, 30, 0), 3e-2: (131, 19, 572, 11, 14)}
T1_BETA_ACCEPTED = {1e-4: 2584, 1e-3: 508, 3e-3: 184, 1e-2: 58, 3e-2: 25}
# macTEM ADR 024, at s = 0 (q = g of each family), sigma = u q: the number of families that are accepted
BETA_S0_ACCEPTED = {3e-3: 23, 1e-2: 8, 3e-2: 4}
BETA_S0_NAMES_1E2 = ["{2 0 0}", "{0 0 1}", "{2 0 -1}", "{2 0 1}", "{1 1 0}", "{1 1 -1}", "{3 1 0}", "{3 1 -1}"]
HAZARD_Z, HAZARD_MARGIN, HAZARD_REACH = 1e-6, 1e-6, 1e-9      # macTEM's redraw distances

FAILURES = []


def check(cond, msg):
    if not cond:
        FAILURES.append(msg)
        print("FAIL:", msg)
    return cond


def num(x):
    """JSON-safe number: non-finite floats as the strings the Swift decoder maps back."""
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Infinity" if x > 0 else "-Infinity"
    return x


# -------------------------------------------------------------------------------------------------------- the decision
def ok(x):
    return math.isfinite(x) and x > 0


def decide(g, q, sigma, k, kmax):
    """The text, NumPy over the family array g (increasing). Returns a dict (fields absent when not part of the verdict)."""
    q, sigma, k, kmax = float(q), float(sigma), float(k), float(kmax)
    if not (ok(q) and ok(sigma) and ok(k) and ok(kmax)):
        return {"v": "invalidInput"}
    F = len(g)
    r = math.sqrt(k * k + 2 * LN19)
    if F == 0:
        reach = q + r * sigma
        return {"v": "beyondKMax", "reach": reach} if reach > kmax else {"v": "noneWithinTolerance"}
    z = np.abs(q - g) / sigma
    order = np.lexsort((g, z))                      # by z, ties by g
    nw = int((z <= k).sum())
    first = int(order[0])
    second = int(order[1]) if F > 1 else -1
    if nw >= 2:
        return {"v": "ambiguous", "p": [int(i) for i in order[:nw]]}
    margin = None
    if second >= 0:
        margin = (z[second] * z[second] - z[first] * z[first]) / 2
        if nw == 1 and margin < LN19:
            return {"v": "ambiguous", "p": [first, second]}
    if nw == 0:
        reach = q + max(r * sigma, abs(q - g[first]))
        if reach > kmax:
            return {"v": "beyondKMax", "reach": reach}
        return {"v": "noneWithinTolerance", "nr": first, "nz": float(z[first]), "nd": float(g[first] / q - 1)}
    reach = q + max(r * sigma, abs(q - g[first]), abs(q - g[second]) if second >= 0 else 0.0)
    if reach > kmax:
        return {"v": "beyondKMax", "reach": reach}
    out = {"v": "accepted", "b": first, "z": float(z[first]), "d": float(g[first] / q - 1)}
    if second >= 0:
        out.update({"n": second, "m": float(margin)})
    return out


def decide_scalar(g, q, sigma, k, kmax):
    """Plain Python, one candidate list at a time, as the text reads (cross-check of `decide`)."""
    if not (ok(q) and ok(sigma) and ok(k) and ok(kmax)):
        return {"v": "invalidInput"}
    cands = sorted(((abs(q - gj) / sigma, gj, j) for j, gj in enumerate(g)))
    within = [c for c in cands if c[0] <= k]
    nxt = cands[1] if len(cands) > 1 else None
    if len(within) >= 2:
        return {"v": "ambiguous", "p": [c[2] for c in within]}
    if len(within) == 1 and nxt is not None:
        if (nxt[0] * nxt[0] - within[0][0] * within[0][0]) / 2 < LN19:
            return {"v": "ambiguous", "p": [within[0][2], nxt[2]]}
    r = math.sqrt(k * k + 2 * LN19)
    named = []
    if cands:
        named.append(cands[0])
        if len(within) == 1 and nxt is not None:
            named.append(nxt)
    R = max([r * sigma] + [abs(q - c[1]) for c in named])
    if q + R > kmax:
        return {"v": "beyondKMax", "reach": q + R}
    if not within:
        if cands:
            return {"v": "noneWithinTolerance", "nr": cands[0][2], "nz": cands[0][0], "nd": cands[0][1] / q - 1}
        return {"v": "noneWithinTolerance"}
    b = within[0]
    out = {"v": "accepted", "b": b[2], "z": b[0], "d": b[1] / q - 1}
    if nxt is not None:
        out.update({"n": nxt[2], "m": (nxt[0] * nxt[0] - b[0] * b[0]) / 2})
    return out


def agree(a, b):
    if a["v"] != b["v"]:
        return False
    for key in ("p", "b", "n", "nr"):
        if a.get(key) != b.get(key):
            return False
    for key in ("z", "d", "m", "nz", "nd", "reach"):
        if (key in a) != (key in b):
            return False
        if key in a and abs(a[key] - b[key]) > 1e-12 * max(1.0, abs(b[key])):
            return False
    return True


def hazards(g, q, sigma, k, kmax):
    """Distances of a plant from the three edges the Swift side could round across: z to k, the margin to ln 19 (one
    within), q + R to k_max."""
    z = np.abs(q - g) / sigma
    hk = float(np.min(np.abs(z - k))) if len(g) else math.inf
    order = np.lexsort((g, z))
    hm = math.inf
    if len(g) > 1 and int((z <= k).sum()) == 1:
        hm = abs(float((z[order[1]] ** 2 - z[order[0]] ** 2) / 2 - LN19))
    d = decide(g, q, sigma, k, kmax)
    hr = abs(d["reach"] - kmax) if "reach" in d else math.inf
    if d["v"] in ("accepted", "noneWithinTolerance") and len(g):
        r = math.sqrt(k * k + 2 * LN19) * sigma
        top = max(r, abs(q - g[order[0]]), abs(q - g[order[1]]) if (d["v"] == "accepted" and len(g) > 1) else 0.0)
        hr = abs(q + top - kmax)
    return hk, hm, hr


# ----------------------------------------------------------------------------------------------------------- the lists
def main():
    with open(PINS2) as fh:
        pins2 = json.load(fh)
    by_name = {s["name"]: s for s in pins2["structures"]}
    lists = []                                       # {name, source, kmax, names, g}

    def add_list(name, source, kmax, fams):
        lists.append({"name": name, "source": source, "kmax": kmax, "names": [f["name"] for f in fams],
                      "g": [f["g"] for f in fams]})

    base = [("Au", "Au"), ("Al", "Al"), ("Si", "Si"), ("Mg", "Mg"), ("Mg-ideal", "ideal_hcp"), ("Weak-cut", "Weak-cut"),
            ("beta''", "beta_double_prime")]
    for name, src in base:
        s = by_name[src]
        check(s["reflections"]["k_max"] == KMAX, f"{name}: slice-2 pins k_max {s['reflections']['k_max']}")
        add_list(name, src, KMAX, s["families"])
        gs = [f["g"] for f in s["families"]]
        check(all(a < b for a, b in zip(gs, gs[1:])), f"{name}: family g not strictly increasing")
    check([len(l["g"]) for l in lists] == [8, 8, 10, 24, 23, 10, 417], f"list sizes {[len(l['g']) for l in lists]}")
    au = next(l for l in lists if l["name"] == "Au")
    for kmax, name in ((0.45, "Au@0.45"), (0.3, "Au@0.3")):         # the one-family and the empty list (macTEM's T4)
        keep = [i for i, g in enumerate(au["g"]) if g <= kmax]
        lists.append({"name": name, "source": "Au", "kmax": kmax, "names": [au["names"][i] for i in keep],
                      "g": [au["g"][i] for i in keep]})
    L = {l["name"]: l for l in lists}
    check(len(L["Au@0.45"]["g"]) == 1 and len(L["Au@0.3"]["g"]) == 0, "Au sub-lists: one family at 0.45, none at 0.3")
    G = {name: np.array(l["g"], dtype=np.float64) for name, l in L.items()}

    def fam_of(name, hkl):
        s = by_name[L[name]["source"]]
        h, nh = list(hkl), [-x for x in hkl]
        for j, f in enumerate(s["families"]):
            if h in f["members"] or nh in f["members"]:
                return j
        raise KeyError((name, hkl))

    # ------------------------------------------------------------------------------------------------------------- T1
    rows = []
    all_plants = []                                  # (list, g array index, q, sigma, k) for cross-checks and hazards

    def record(row, g, kmax, q, sigma, k):
        d = decide(g, q, sigma, k, kmax)
        d2 = decide_scalar(list(map(float, g)), float(q), float(sigma), float(k), float(kmax))
        check(agree(d, d2), f"scalar vs NumPy decision differ: {row}")
        row.update(d)
        return d

    counts = {}                                      # (list, u) -> Counter-like dict
    hz = [math.inf, math.inf, math.inf]
    s0_accepted = {}                                 # (list, u) -> names accepted at s = 0
    for ci, (name, _) in enumerate(base):
        g = G[name]
        for j in range(len(g)):
            for u in U1:
                for s in S1:
                    q = float(g[j]) * (1 + s * u)
                    sigma = u * q
                    row = {"c": ci, "t": j, "u": u, "s": s, "q": q, "sg": sigma}
                    d = record(row, g, KMAX, q, sigma, K)
                    # P1 (closed forms): the truth family's z and delta under the plant
                    if d["v"] == "accepted" and d["b"] == j:
                        check(abs(d["z"] - abs(s) / (1 + s * u)) <= 1e-9 * max(1, d["z"]),
                              f"P1 z: {name} family {j} u {u} s {s}")
                        check(abs(d["d"] + s * u / (1 + s * u)) <= 1e-14, f"P1 delta: {name} family {j} u {u} s {s}")
                    hk, hm, hr = hazards(g, q, sigma, K, KMAX)
                    hz = [min(hz[0], hk), min(hz[1], hm), min(hz[2], hr)]
                    rows.append(row)
                    c = counts.setdefault((name, u), {"own": 0, "other": 0, "amb": 0, "bey": 0, "none": 0, "inv": 0})
                    v = d["v"]
                    if v == "accepted":
                        c["own" if d["b"] == j else "other"] += 1
                        if s == 0.0:
                            s0_accepted.setdefault((name, u), []).append(L[name]["names"][d["b"]])
                    else:
                        c[{"ambiguous": "amb", "beyondKMax": "bey", "noneWithinTolerance": "none",
                           "invalidInput": "inv"}[v]] += 1
    n_t1 = len(rows)
    check(n_t1 == (83 + 417) * 45, f"T1 has {n_t1} rows")
    print(f"T1: {n_t1} plants; closest |z - k| {hz[0]:.3e}, closest margin to ln 19 {hz[1]:.3e}, "
          f"closest q + R to k_max {hz[2]:.3e}")
    check(hz[0] > HAZARD_Z and hz[1] > HAZARD_MARGIN and hz[2] > HAZARD_REACH,
          "a T1 plant lies within macTEM's redraw distance of k, ln 19 or k_max")

    # macTEM's registered counts: the six base lists summed per u
    print("T1 over macTEM's six lists (own, other, ambiguous, beyond, none) | registered")
    six = ("Au", "Al", "Si", "Mg", "Mg-ideal", "Weak-cut")
    for u in U1:
        got = tuple(sum(counts[(n, u)][key] for n in six) for key in ("own", "other", "amb", "bey", "none"))
        same = got == T1_PREDICT[u]
        print(f"  u {u:<6g} {got} | {T1_PREDICT[u]}  {'EQUAL' if same else 'DIFFERENT'}")
        check(same, f"T1 counts over the six lists at u={u}: {got}, macTEM registered {T1_PREDICT[u]}")
    print("beta'' accepted in all per u | registered")
    for u in U1:
        c = counts[("beta''", u)]
        got = c["own"] + c["other"]
        print(f"  u {u:<6g} {got} (own {c['own']}, other {c['other']}, ambiguous {c['amb']}, beyond {c['bey']}, "
              f"none {c['none']}) | {T1_BETA_ACCEPTED[u]}")
        check(got == T1_BETA_ACCEPTED[u], f"beta'' accepted at u={u}: {got}, registered {T1_BETA_ACCEPTED[u]}")
    print("beta'' families accepted at s = 0 | registered (ADR 024)")
    for u in (3e-3, 1e-2, 3e-2):
        got = s0_accepted.get(("beta''", u), [])
        print(f"  u {u:<6g} {len(got)} | {BETA_S0_ACCEPTED[u]}")
        check(len(got) == BETA_S0_ACCEPTED[u], f"beta'' at s=0, u={u}: {len(got)} accepted, registered {BETA_S0_ACCEPTED[u]}")
    check(sorted(s0_accepted.get(("beta''", 1e-2), [])) == sorted(BETA_S0_NAMES_1E2),
          f"beta'' s=0 u=1e-2 names {s0_accepted.get(('beta' + chr(39) + chr(39), 1e-2))}")

    # ------------------------------------------------------------------------------------------------------------- T4
    t4 = []

    def plant(label, name, q, sigma, k, expect, truth=None):
        g = G[name]
        row = {"label": label, "l": name, "q": num(q), "sg": num(sigma), "k": num(k)}
        d = decide(g, q, sigma, k, L[name]["kmax"])
        d2 = decide_scalar(list(map(float, g)), float(q), float(sigma), float(k), float(L[name]["kmax"]))
        check(agree(d, d2), f"scalar vs NumPy decision differ: {label}")
        check(d["v"] == expect, f"T4 '{label}': {d['v']}, expected {expect}")
        if truth is not None and d["v"] == "accepted":
            check(d["b"] == truth, f"T4 '{label}': accepted family {d.get('b')}, truth {truth}")
        row.update(d)
        row["expect"] = expect
        t4.append(row)
        return d

    g111 = float(G["Au"][0])
    for badname, bad in (("0", 0.0), ("-1", -1.0), ("NaN", math.nan), ("+inf", math.inf), ("-inf", -math.inf)):
        plant(f"invalid q = {badname}", "Au", bad, 1e-3 * g111, K, "invalidInput")
        plant(f"invalid sigma = {badname}", "Au", g111, bad, K, "invalidInput")
        plant(f"invalid k = {badname}", "Au", g111, 1e-3 * g111, bad, "invalidInput")
    # beyondKMax
    j = fam_of("Mg", (2, 0, -5))
    check(L["Mg"]["names"][j] == "{2 0 -2 5}" and j == len(G["Mg"]) - 1, "Mg {20-25} is the last family")
    plant("beyond: Mg {20-25} at its own g, u = 1e-3", "Mg", G["Mg"][j], 1e-3 * G["Mg"][j], K, "beyondKMax")
    j = fam_of("Si", (6, 2, 0))
    plant("beyond: Si {620} at its own g, u = 1e-3", "Si", G["Si"][j], 1e-3 * G["Si"][j], K, "beyondKMax")
    plant("beyond: q = 1.25 on Au, u = 1e-3", "Au", 1.25, 1.25e-3, K, "beyondKMax")
    # the one-family list and the empty one
    g0 = G["Au@0.45"][0]
    plant("one-family list (Au, k_max 0.45): {111} at its own g, u = 1e-3 -> accepted, next nil", "Au@0.45", g0, 1e-3 * g0, K,
          "accepted", truth=0)
    plant("one-family list: q = 0.3, u = 1e-3 -> none, nearest {111}", "Au@0.45", 0.3, 3e-4, K, "noneWithinTolerance")
    plant("one-family list: q = 0.449, u = 1e-3 -> beyond", "Au@0.45", 0.449, 0.449e-3, K, "beyondKMax")
    plant("empty list (Au, k_max 0.3): q = 0.2, u = 1e-3 -> none, no nearest", "Au@0.3", 0.2, 2e-4, K, "noneWithinTolerance")
    plant("empty list: q = 0.299, u = 1e-3 -> beyond", "Au@0.3", 0.299, 0.299e-3, K, "beyondKMax")
    # noneWithinTolerance: forbidden positions
    a_si, c_mg, a_au = 5.4309, 5.2108, 4.0782
    plant("none: Si {200} (forbidden) at 2/a", "Si", 2 / a_si, 1e-3 * 2 / a_si, K, "noneWithinTolerance")
    plant("none: Si {222} (forbidden) at sqrt(12)/a", "Si", math.sqrt(12) / a_si, 1e-3 * math.sqrt(12) / a_si, K,
          "noneWithinTolerance")
    plant("none: Mg {0001} (forbidden) at 1/c", "Mg", 1 / c_mg, 1e-3 / c_mg, K, "noneWithinTolerance")
    cell = by_name["Mg"]["cell"]
    A = np.array([[cell["a"], 0, 0],
                  [cell["b"] * math.cos(math.radians(cell["gamma"])), cell["b"] * math.sin(math.radians(cell["gamma"])), 0],
                  [0, 0, cell["c"]]])
    g1121 = float(np.linalg.norm(np.array([1.0, 1.0, 1.0]) @ np.linalg.inv(A).T))
    plant("none: Mg {11-21} (forbidden)", "Mg", g1121, 1e-3 * g1121, K, "noneWithinTolerance")
    plant("none: Au {100} (forbidden) at 1/a", "Au", 1 / a_au, 1e-3 / a_au, K, "noneWithinTolerance")
    plant("none: Au {111} position measured against Al's list", "Al", g111, 1e-3 * g111, K, "noneWithinTolerance")
    # ambiguous
    ja, jb = fam_of("Mg", (1, 0, 5)), fam_of("Mg", (2, 1, 2))
    qa = G["Mg"][ja] + 0.4 * (G["Mg"][jb] - G["Mg"][ja])
    plant("ambiguous: Mg {10-15}/{21-32}, 40 % of the way, u = 1e-3", "Mg", qa, 1e-3 * qa, K, "ambiguous")
    ja, jb = fam_of("Mg", (3, 0, 2)), fam_of("Mg", (0, 0, 6))
    qa = G["Mg"][ja] + 0.4 * (G["Mg"][jb] - G["Mg"][ja])
    plant("ambiguous: Mg {30-32}/{0006}, 40 % of the way, u = 3e-3", "Mg", qa, 3e-3 * qa, K, "ambiguous")
    j331, j420 = fam_of("Au", (3, 3, 1)), fam_of("Au", (4, 2, 0))
    g331, g420 = G["Au"][j331], G["Au"][j420]
    plant("ambiguous: Au {331} at its own g, u = 1e-2 ({420} 2.6 % away)", "Au", g331, 1e-2 * g331, K, "ambiguous", truth=j331)
    sg = (g420 - g331) / 5.7
    d = plant("ambiguous: the margin plant, Au {331} -> {420}, sigma = gap/5.7, z = 2.5 (next 3.2, margin 2.0 < ln 19)", "Au",
              g331 + 2.5 * sg, sg, K, "ambiguous", truth=j331)
    check(d.get("p") == [j331, j420], "the margin plant names {331} and {420}")
    bq = "beta''"
    gp, gr = G[bq][fam_of(bq, (3, 1, 5))], G[bq][fam_of(bq, (7, 3, 0))]
    qm = (gp + gr) / 2
    d = plant("ambiguous: beta'' {3 1 5}/{7 3 0} (1.03e-5 apart), midway, u = 1e-3", bq, qm, 1e-3 * qm, K, "ambiguous")
    check({fam_of(bq, (3, 1, 5)), fam_of(bq, (7, 3, 0))} <= set(d.get("p", [])), "beta'' pair: both named")
    check(abs(abs(gp - gr) / min(gp, gr) - 1.03e-5) <= 5e-8, "beta'' pair is 1.03e-5 apart")
    jf = fam_of(bq, (2, 0, 0))
    d = plant("accepted: beta'' {2 0 0} at its own g, u = 1e-3, next {0 0 1}", bq, G[bq][jf], 1e-3 * G[bq][jf], K, "accepted",
              truth=jf)
    check(L[bq]["names"][d["n"]] == "{0 0 1}", "beta'' {2 0 0}: the next is {0 0 1}")
    # the boundary the NumPy side must not be sitting on
    for r in t4:
        if r["v"] != "invalidInput":
            hk, hm, hr = hazards(G[r["l"]], r["q"], r["sg"], r["k"], L[r["l"]]["kmax"])
            check(hk > 1e-9 and hm > 1e-9 and hr > 1e-12, f"T4 '{r['label']}' sits on a boundary {hk:.2e} {hm:.2e} {hr:.2e}")

    # ----------------------------------------------------------------------------------------------------------- write
    out = {
        "versions": {"python": sys.version.split()[0], "numpy": np.__version__},
        "definitions": {
            "text": "z = |q - g| / sigma; delta = g / q - 1; candidates by (z, g); within = z <= k; ambiguous: >= 2 within "
                    "(payload all, in (z, g) order) or one within with margin (z_next^2 - z^2) / 2 < ln 19 (payload "
                    "[best, next]); r = sqrt(k^2 + 2 ln 19); beyondKMax: q + max(r sigma, |q - g| of the named families) > "
                    "k_max; noneWithinTolerance; accepted. T1: q = g (1 + s u), sigma = u q, k = 3.",
            "k": K, "ln19": LN19, "k_max": KMAX, "u": list(U1), "s": list(S1),
            "tolerances": {"delta": 1e-14, "z": "1e-10 max(1, z)", "margin": "1e-9 max(1, margin)",
                           "reach": "1e-9 max(1, reach)"},
            "source_of_lists": "Tests/DiffractoScatteringTests/Fixtures/numpy-structure-factor-pins.json (families after the |F| cut)",
            "t1_columns": {"c": "list index", "t": "truth family", "u": "sigma/q", "s": "(q - g)/(g u)", "q": "q", "sg": "sigma",
                           "v": "verdict", "b": "best family", "z": "best z", "d": "best delta", "n": "next family",
                           "m": "margin", "p": "payload (ambiguous): family indices in (z, g) order", "nr": "nearest (none)",
                           "nz": "nearest z", "nd": "nearest delta", "reach": "q + R (beyondKMax)"},
        },
        "lists": lists,
        "t1_lists": [name for name, _ in base],
        "t1": rows,
        "t4": t4,
        "t1_counts": {f"{n}@{u}": c for (n, u), c in counts.items()},
        "beta_s0_accepted": {str(u): s0_accepted.get(("beta''", u), []) for u in U1},
        "hazards": {"min_abs_z_minus_k": hz[0], "min_abs_margin_minus_ln19": hz[1], "min_abs_reach_minus_kmax": hz[2]},
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, separators=(",", ":"), allow_nan=False)
    size = os.path.getsize(OUT)
    print(f"wrote {OUT}: {size / 1e6:.2f} MB, {len(rows)} T1 plants, {len(t4)} T4 plants, {len(lists)} lists")
    if FAILURES:
        print(f"{len(FAILURES)} self-check failure(s)")
        sys.exit(1)
    print("every self-check holds")


if __name__ == "__main__":
    main()
