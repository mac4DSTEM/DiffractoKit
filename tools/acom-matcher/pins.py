"""DiffractoKit ACOM matcher: the pin generator (docs/parity.md row P25).

    tools/memcap.sh --gb 4 -- References/CuboKit/py4dstem-venv/bin/python -I -B tools/acom-matcher/pins.py

writes Tests/DiffractoOrientationTests/Fixtures/acom-matcher-pins.json. The py4DSTEM venv is the one of CuboKit (py4DSTEM
f050d207 = 0.14.19, NumPy, SciPy; versions are written into the fixture). Nothing here is generated from the new Swift.

Two kinds of pin, as docs/parity.md says:

  APP ROUTE ("transcription of the app source at mac4DSTEM c809679c; the three files' content is identical at dc9bdaa6, see the SHA-256s in the fixture"): app_route.py is a NumPy transcription of the apps'
  OrientationPlan / OrientationMatcher (polar image, ring FFTs, correlation, mirror pass, runner-up, orientation matrix). Per
  crystal (m-3m Al, 6/mmm Mg, 2/m monoclinic) it pins the plan's beam directions and one template (image and spots), and for
  planted patterns (the apps' own excitation model, jittered peaks, a Q calibration) the experimental polar image, the ring
  FFT of one, every template's score and bin, the winner, the runner-up, the in-plane angle and the orientation matrix.

  PY4DSTEM (Crystal f050d207, run): (1) the spots of a zone axis, generate_diffraction_pattern at lambda -> 0 (setup_diffraction
  (1e30)) and at 200 kV, which the kit's `project` must reproduce in the excitation models `.flat` and `.py4DSTEM`; (2) planted
  orientations: py4DSTEM's own peaks at 200 kV, the orientation py4DSTEM's orientation_plan + match_single_pattern returns for
  them, to compare with the route in the kit. The transcription is also run on the same peaks with py4DSTEM's excitation
  model (the kit's option; the numbers it gives are measurements for the notes, the test bounds come from them).
"""
import base64
import hashlib
import json
import math
import os
import subprocess
import sys

import numpy as np
import scipy
import py4DSTEM
from py4DSTEM.process.diffraction import Crystal
from scipy.spatial.transform import Rotation as Rot

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import app_route as A  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ORG = os.path.abspath(os.path.join(ROOT, "..", ".."))
OUT = os.path.join(ROOT, "Tests", "DiffractoOrientationTests", "Fixtures", "acom-matcher-pins.json")
APP = os.path.join(ORG, "mac4DSTEM")
APP_FILES = ["mac4DSTEM/Core/Crystal/OrientationPlan.swift", "mac4DSTEM/Core/Crystal/OrientationMatcher.swift",
             "mac4DSTEM/Core/Analysis/OrientationResult.swift"]
PY_SRC = os.path.join(ORG, "References", "CuboKit", "py4DSTEM-f050d207")
PY_LOCK = "f050d207851f8accd3aa550365baaf1f38230428"
SIGMA = 0.03 / math.sqrt(2)          # py4DSTEM's sigma for the kit's w = 0.03
SG_MAX = 0.1


def b64(a, dtype):
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode()


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def provenance():
    inst = os.path.dirname(py4DSTEM.__file__)
    for f in ["process/diffraction/crystal.py", "process/diffraction/crystal_ACOM.py"]:
        assert sha(os.path.join(inst, f)) == sha(os.path.join(PY_SRC, "py4DSTEM", f)), f
    head = open(os.path.join(PY_SRC, ".git", "HEAD")).read().strip()
    app_head = subprocess.run(["git", "-C", APP, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    return {
        "py4dstem": py4DSTEM.__version__, "py4dstem_lock": PY_LOCK, "py4dstem_tree_head": head, "numpy": np.__version__,
        "scipy": scipy.__version__, "python": sys.version.split()[0],
        "app_transcription": "transcription of the app source at mac4DSTEM " + app_head[:8]
                             + " (the three files below were not modified in its working tree at the time)",
        "app_files_sha256": {f: sha(os.path.join(APP, f)) for f in APP_FILES},
    }


# ---------------------------------------------------------------------------------------------------------------------
# crystals

def crystals():
    s3 = math.sqrt(3)
    return {
        "al": dict(name="Al fcc (m-3m)", laue="m-3m", cell=[4.05, 4.05, 4.05, 90, 90, 90], kmax=1.6,
                   pos=[[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]], z=[13] * 4, count=600,
                   zone=[[0, 1, 1], [1, 1, 1]], step=2.0, sampler=A.sample_cubic, fold=A._cubic_fold),
        "mg": dict(name="Mg hcp (6/mmm)", laue="6/mmm", cell=[3.21, 3.21, 5.21, 90, 90, 120], kmax=1.6,
                   pos=[[1 / 3, 2 / 3, .25], [2 / 3, 1 / 3, .75]], z=[12, 12], count=300,
                   zone=[[s3 / 2, .5, 0], [1, 0, 0]], step=3.0, sampler=A.sample_hex, fold=A._hex_fold),
        "mono": dict(name="synthetic monoclinic b-unique (2/m)", laue="2/m", cell=[4.2, 3.5, 5.0, 90, 100, 90], kmax=1.2,
                     pos=[[0, 0, 0], [.3, 0, .2], [.7, 0, .8], [.15, .5, .4], [.85, .5, .6]], z=[14, 8, 8, 8, 8], count=300,
                     zone="half", step=3.0, sampler=A.sample_mono, fold=A._mono_fold),
        # no transcription for this one (the apps fall back to unreduced beam sampling for it): py4DSTEM comparison only
        "rutile": dict(name="rutile-like tetragonal (4/mmm)", laue="4/mmm", cell=[4.59, 4.59, 2.96, 90, 90, 90], kmax=1.4,
                       pos=[[0, 0, 0], [.5, .5, .5], [.3, .3, 0], [.7, .7, 0], [.2, .8, .5], [.8, .2, .5]],
                       z=[22, 22, 8, 8, 8, 8], count=None, zone=[[1, 0, 0], [1, 1, 0]], step=3.0, sampler=None, fold=None),
    }


def build(c):
    cr = Crystal(np.array(c["pos"]), c["z"], c["cell"])
    cr.calculate_structure_factors(k_max=c["kmax"])
    g = cr.g_vec_all.T.copy()
    inten = cr.struct_factors_int.astype(np.float64).copy()
    hkl = cr.hkl.T.copy().astype(np.int64)
    return cr, g, inten, hkl


def plan_py4dstem(cr, c):
    cr.orientation_plan(zone_axis_range=np.array(c["zone"]) if not isinstance(c["zone"], str) else c["zone"],
                        angle_step_zone_axis=c["step"], angle_step_in_plane=c["step"], accel_voltage=200e3,
                        sigma_excitation_error=SIGMA, progress_bar=False)


def ang_deg(M1, M2, ops):
    best = 360.0
    for S in ops:
        c = (np.trace(M2.T @ S @ M1) - 1) / 2
        best = min(best, math.degrees(math.acos(max(-1.0, min(1.0, c)))))
    return best


def proper_ops(laue):
    import itertools
    if laue == "m-3m":
        return A.cubic_operators()
    ops = [np.eye(3)]
    if laue == "2/m":
        ops.append(np.diag([-1.0, 1, -1]))
    if laue == "6/mmm":
        ops = []
        for k in range(6):
            a = k * math.pi / 3
            ops.append(np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]]))
        for k in range(6):
            a = k * math.pi / 6
            x, y = math.cos(a), math.sin(a)
            ops.append(np.array([[2 * x * x - 1, 2 * x * y, 0], [2 * x * y, 2 * y * y - 1, 0], [0, 0, -1]]))
    if laue == "4/mmm":
        ops = []
        for k in range(4):
            a = k * math.pi / 2
            ops.append(np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]]))
        for k in range(4):
            a = k * math.pi / 4
            x, y = math.cos(a), math.sin(a)
            ops.append(np.array([[2 * x * x - 1, 2 * x * y, 0], [2 * x * y, 2 * y * y - 1, 0], [0, 0, -1]]))
    return ops


def euler_zxz(M):
    return Rot.from_matrix(M.T).as_euler("zxz")


# ---------------------------------------------------------------------------------------------------------------------

def app_model_peaks(g, inten, M, lam):
    """The pattern of orientation M under the apps' excitation model (s = g.z + lambda |g|^2 / 2): detector x = g.M[:,0],
    y = g.M[:,1], weight exp(-(s/w)^2) times |F|^2, no power (the matcher applies its own exponent)."""
    out = []
    for gi, ii in zip(g, inten):
        sg = float(gi @ M[:, 2]) + lam / 2 * float(gi @ gi)
        if abs(sg) > SG_MAX:
            continue
        out.append([float(gi @ M[:, 0]), float(gi @ M[:, 1]), ii * math.exp(-(sg / 0.03) ** 2)])
    return np.array(out)


def pixelise(q, rng, scale, origin, jitter, spurious):
    """Peaks as detector values (float32): position q/scale + origin with a Gaussian jitter (1/A) and `spurious` extra peaks.
    scale = 1/A per pixel (peaks in pixels) or 1 (peaks in 1/A)."""
    xy = q[:, :2] + rng.normal(0, jitter, q[:, :2].shape)
    inten = q[:, 2] * 1000.0 * rng.uniform(0.85, 1.15, len(q))
    rows = [[xy[i, 0] / scale + origin[0], xy[i, 1] / scale + origin[1], inten[i]] for i in range(len(q))]
    for _ in range(spurious):
        a = rng.uniform(0, 2 * math.pi)
        r = rng.uniform(0.15, 1.1)
        rows.append([r * math.cos(a) / scale + origin[0], r * math.sin(a) / scale + origin[1], float(np.median(inten) * 0.5)])
    return np.array(rows, dtype=np.float32)


def margins(plan, res, peaks, origin, inv):
    """How far the pinned choice is from a tie: the winner's score minus the best other template's, and the winner's in-plane
    correlation peak minus the best peak at least 4 bins away (so a Swift/NumPy rounding difference cannot flip a pin)."""
    scores = np.array(res["scores"], dtype=np.float64)
    t = res["template"]
    margin_t = scores[t] - np.max(np.delete(scores, t))
    polar = A.experimental_polar(plan, peaks, (np.float32(origin[0]), np.float32(origin[1])), inv)
    nr, na = plan["nr"], plan["na"]
    ere, eim = A.ring_ffts(polar, nr, na)
    E = ere.astype(np.float64) + 1j * eim.astype(np.float64)
    T = plan["fre"][t].astype(np.float64) + 1j * plan["fim"][t].astype(np.float64)
    corr = np.fft.ifft(((np.conj(E) if not res["mirrored"] else E) * T).sum(axis=0)).real
    b = int(np.argmax(corr))
    others = [corr[(b + d) % na] for d in range(4, na - 3)]
    return float(margin_t), float(corr[b] - max(others))


def main():
    prov = provenance()
    out = {"provenance": prov, "sigma": SIGMA, "sgMax": SG_MAX, "crystals": {}}
    cs = crystals()
    for key, c in cs.items():
        print("crystal", key, flush=True)
        cr, g, inten, hkl = build(c)
        cr.setup_diffraction(200e3)
        lam = float(cr.wavelength)
        nonzero = np.linalg.norm(g, axis=1) > 1e-9
        rec = {"name": c["name"], "laue": c["laue"], "cell": c["cell"], "kmax": c["kmax"], "wavelength": lam,
               "g": b64(g[nonzero], "<f8"), "intensity": b64(inten[nonzero], "<f8"), "count": int(nonzero.sum())}
        gz, iz = g[nonzero], inten[nonzero]
        ops = proper_ops(c["laue"])
        rng = np.random.default_rng(20261016 + len(key))

        # ---- the plan of the app route (transcription): axes, one template ------------------------------------------------
        route = None
        if c["sampler"] is not None:
            axes = c["sampler"](c["count"])
            bases = [A.detector_basis(n) for n in axes]
            plan = A.generate_plan(gz, iz, c["kmax"], axes, wavelength=lam)
            rec["appRoute"] = {"zoneAxisCount": c["count"], "axes": b64(np.array(axes), "<f8")}
            tpick = [0, len(axes) // 3] if key == "al" else [len(axes) // 3]
            rec["appRoute"]["templates"] = [
                {"index": t, "image": b64(plan["templates"][t], "<f4"),
                 "spots": b64(np.array(list(zip(*plan["spots"][t]))), "<f8"),
                 "ringFFTRe": b64(plan["fre"][t], "<f4") if key == "al" and t == 0 else None,
                 "ringFFTIm": b64(plan["fim"][t], "<f4") if key == "al" and t == 0 else None}
                for t in tpick]
            # planted patterns
            ipp = 0.005
            cases = []
            ncase = {"al": 6, "mg": 4, "mono": 4}[key]
            k = 0
            while k < ncase:
                Mt = Rot.random(random_state=int(rng.integers(1 << 30))).as_matrix()
                q = app_model_peaks(gz, iz, Mt, lam)
                calib = (k % 2 == 1)
                origin = (256.3 + k * 0.17, 255.7 - k * 0.11) if calib else (0.31 + k * 0.01, -0.22)
                inv = ipp if calib else 1.0
                peaks = pixelise(q, rng, inv, origin, 0.0 if k == 0 else 0.3 * ipp, 0 if k < 2 else 1)
                res = A.match(plan, axes, bases, peaks, (np.float32(origin[0]), np.float32(origin[1])), inv)
                red, _ = A.reduce_min_angle(res["matrix"], ops)
                margin_template, margin_inplane = margins(plan, res, peaks, origin, inv)
                if not (margin_template > 1e-4 and margin_inplane > 1e-4):
                    print('  (a near tie, drawing another pattern)', margin_template, margin_inplane)
                    continue
                case = {"peaks": b64(peaks, "<f4"), "origin": [float(np.float32(origin[0])), float(np.float32(origin[1]))],
                        "invAngstromPerPixel": inv, "trueMatrix": Mt.T.reshape(-1).tolist(),
                        "template": res["template"], "bin": int(res["bin"]), "mirrored": bool(res["mirrored"]),
                        "score": res["score"], "second": res["second"], "angle": res["angle"],
                        "matrix": res["matrix"].T.reshape(-1).tolist(),     # column-major: column 0 first
                        "scores": b64(res["scores"], "<f4"),
                        "errorToTruthDeg": ang_deg(res["matrix"], Mt, ops),
                        "marginTemplate": margin_template, "marginInPlane": margin_inplane}
                if key == "al":
                    case["euler"] = euler_zxz(red).tolist()
                if k == 0:
                    case["polar"] = b64(res["polar"], "<f4")
                    if key == "al":
                        case["ringFFTRe"] = b64(res["ere"], "<f4")
                        case["ringFFTIm"] = b64(res["eim"], "<f4")
                cases.append(case)
                k += 1
                print("  case", k, res["template"], res["mirrored"], round(res["score"], 3), "err", round(case["errorToTruthDeg"], 2),
                      flush=True)
            rec["appRoute"]["cases"] = cases
            rec["appRoute"]["calibrationNote"] = "odd cases pass invAngstromPerPixel = 0.005 with peaks in pixels; even cases pass 1 with peaks in 1/A"
            # the same transcription with py4DSTEM's excitation model (the kit's option), for the end-to-end measurements
            plan_p = A.generate_plan(gz, iz, c["kmax"], axes, wavelength=lam, py4dstem=True)
        else:
            axes = bases = plan_p = None

        # ---- py4DSTEM: spots of a zone axis, flat (lambda -> 0) and 200 kV --------------------------------------------------
        spots = []
        pick_axes = (A.sample_cubic(40) if key == "al" else
                     (c["sampler"](40) if c["sampler"] else
                      [np.array(v) / np.linalg.norm(v) for v in ([0, 0, 1], [1, 0, 0], [1, 1, 0], [.3, .2, .93], [.5, .5, .7])]))
        pick_axes = [np.array(p, dtype=np.float64) for p in pick_axes][::max(1, len(pick_axes) // 4)][:4]
        for mode, volts in (("flat", 1e30), ("curved", 200e3)):
            cr.setup_diffraction(volts)
            for n in pick_axes:
                basis = A.detector_basis(n)
                pl, _ = cr.generate_diffraction_pattern(orientation_matrix=basis, sigma_excitation_error=SIGMA,
                                                        tol_excitation_error_mult=SG_MAX / SIGMA, tol_intensity=0.0,
                                                        k_max=c["kmax"], return_orientation_matrix=True)
                d = pl.data
                m = ~((d["h"] == 0) & (d["k"] == 0) & (d["l"] == 0))
                spots.append({"mode": mode, "wavelength": float(cr.wavelength), "axis": n.tolist(),
                              "sgMax": float(SIGMA * (SG_MAX / SIGMA)),
                              "qx": b64(d["qx"][m], "<f8"), "qy": b64(d["qy"][m], "<f8"),
                              "intensity": b64(d["intensity"][m], "<f8"), "hkl": b64(np.stack([d["h"][m], d["k"][m], d["l"][m]], 1), "i1")})
        rec["py4dstemSpots"] = spots

        # ---- py4DSTEM: planted orientation, its own peaks, its own match ----------------------------------------------------
        cr.setup_diffraction(200e3)
        plan_py4dstem(cr, c)
        e2e = []
        nplant = {"al": 6, "mg": 4, "mono": 4, "rutile": 4}[key]
        for k in range(nplant):
            Mt = Rot.random(random_state=int(rng.integers(1 << 30))).as_matrix()
            pl = cr.generate_diffraction_pattern(orientation_matrix=Mt, sigma_excitation_error=SIGMA,
                                                 tol_excitation_error_mult=SG_MAX / SIGMA, tol_intensity=0.0, k_max=c["kmax"])
            d = pl.data
            m = ~((d["h"] == 0) & (d["k"] == 0) & (d["l"] == 0))
            o = cr.match_single_pattern(pl, verbose=False)
            Mp = o.matrix[0]
            rec_e = {"trueMatrix": Mt.T.reshape(-1).tolist(), "qx": b64(d["qx"][m], "<f8"), "qy": b64(d["qy"][m], "<f8"),
                     "intensity": b64(d["intensity"][m], "<f8"), "py4dstemMatrix": Mp.T.reshape(-1).tolist(),
                     "py4dstemCorr": float(o.corr[0]), "py4dstemErrorToTruthDeg": ang_deg(Mp, Mt, ops)}
            if plan_p is not None:
                peaks = np.stack([d["qx"][m], d["qy"][m], d["intensity"][m]], 1).astype(np.float32)
                r2 = A.match(plan_p, axes, bases, peaks, (np.float32(0), np.float32(0)), 1.0)
                rec_e["transcriptionPy4dstemModel"] = {
                    "matrix": r2["matrix"].T.reshape(-1).tolist(), "score": r2["score"], "template": r2["template"],
                    "errorToTruthDeg": ang_deg(r2["matrix"], Mt, ops), "errorToPy4dstemDeg": ang_deg(r2["matrix"], Mp, ops)}
            print("  e2e", k, "py4dstem err", round(rec_e["py4dstemErrorToTruthDeg"], 2),
                  rec_e.get("transcriptionPy4dstemModel", {}).get("errorToTruthDeg"), flush=True)
            e2e.append(rec_e)
        rec["planted"] = e2e
        out["crystals"][key] = rec

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print("wrote", OUT, os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
