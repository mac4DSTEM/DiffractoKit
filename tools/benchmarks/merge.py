"""Merges results/<date>.jsonl (kit and Python lines) into the per-function table of docs/benchmarks.md.

  python3 merge.py results/2026-10-21.jsonl > /tmp/table.md
Prints the table body; the header, the backend decisions and the notes are written by hand in docs/benchmarks.md.
"""
import json, sys
rows = [json.loads(l) for l in open(sys.argv[1]) if l.startswith("{")]
order = ["niggli_reduce", "reflections", "kinematic_reflections", "symmetry_detection", "cif_import", "orientation_reduce",
         "ellipse_fit_1d", "acom_library", "acom_match_scan"]
def per(r, name):
    """acom_match_scan is compared per pattern: the kit runs the whole scan, Python the first n patterns."""
    if name != "acom_match_scan" or r.get("median_ms") is None:
        return r
    r = dict(r)
    for k in ("median_ms", "min_ms", "max_ms"):
        r[k] = r[k] / r["n"]
    return r
fmt = lambda r: "nightly" if r is None or r.get("median_ms") is None else f"{r['median_ms']:.3g} ms ({r['min_ms']:.3g}-{r['max_ms']:.3g})"
print("| Function | Size | n | Kit | Python | Kit / Python (first Python row) |")
print("|---|---|---|---|---|---|")
for name in order:
    for size in ("small", "typical", "large"):
        kit = [per(r, name) for r in rows if r["side"] == "kit" and r["name"] == name and r["size"] == size]
        py = [per(r, name) for r in rows if r["side"] == "python" and r["name"] == name and r["size"] == size and "median_ms" in r]
        if not kit and not py:
            continue
        kcell = "<br>".join(f"{fmt(r)} [{r['backend']}]" if len(kit) > 1 else fmt(r) for r in kit) or "-"
        pcell = "<br>".join(f"{fmt(r)} [{r['backend']}]" if len(py) > 1 else fmt(r) for r in py) or "-"
        ratio = "-"
        kk = [r for r in kit if r.get("median_ms")]; pp = [r for r in py if r.get("median_ms")]
        if kk and pp:
            ratio = " / ".join(f"{k['median_ms'] / pp[0]['median_ms']:.2g}" for k in kk)
        n = (kit or py)[0]["n"]
        print(f"| `{name}` | {size} | {n} | {kcell} | {pcell} | {ratio} |")
