"""DiffractoKit: the py4DSTEM half of the pin generator.

Not run by hand: tools/pins/phase-vectors/pins.py starts it with the py4DSTEM venv's interpreter, sends the phases as JSON on stdin and
reads the JSON on stdout. Run through pins.py the whole chain sits under ONE memory cap.

    References/python/py4dstem/bin/python -I -B tools/pins/phase-vectors/pins-py4dstem.py < request.json > answer.json

py4DSTEM f050d207 (0.14.19) is only called, never edited: Crystal(...), calculate_structure_factors, setup_diffraction,
generate_diffraction_pattern (keep_qz) and parse_orientation. Per phase it answers

  * its lattice (lat_real, lat_inv) and the structure-factor table it computed (hkl, g, the complex64 F, |F|^2 as float32);
  * per zone axis the orientation matrix of parse_orientation (columns x', y', z') and the diffraction pattern at
    (a) lambda -> 0: setup_diffraction(1e30 V), sigma = w / sqrt(2), tol_excitation_error_mult = slab / sigma,
        tol_intensity 0, k_max = kMax: with lambda ~ 1e-27 the excitation error is s = g.n and the weight
        exp(-s^2 / (2 sigma^2)) is exp(-(s/w)^2), the kit's flat excitation;
    (b) 200 kV (setup_diffraction(200e3)): the curved Ewald sphere, py4DSTEM as it is used (a parity row, not a gate).

Both runs list every peak py4DSTEM returns (qx, qy, qz = the g vector in the projected frame, the intensity, h, k, l),
including (0 0 0); pins.py drops it.
"""
import hashlib
import json
import math
import os
import sys

import numpy as np
import scipy
import py4DSTEM
from py4DSTEM.process.diffraction import Crystal

SRC = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "References", "CuboKit",
                                   "py4DSTEM-f050d207"))
LOCK = "f050d207851f8accd3aa550365baaf1f38230428"
TABLE_SHA256 = "be845d1d541cdb75ae16f8928cd9e56a9ea38af8a3cc43da9137e73cc72bb6f4"


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def provenance():
    """The installed py4DSTEM is the source tree at the lock: the files this script's numbers come from are byte equal."""
    inst = os.path.dirname(py4DSTEM.__file__)
    files = ["process/diffraction/crystal.py", "process/utils/single_atom_scatter.py", "process/utils/utils.py",
             "process/utils/scattering_factors.txt"]
    head = open(os.path.join(SRC, ".git", "HEAD")).read().strip()
    out = {"lock": LOCK, "tree_head": head, "files": {}}
    for f in files:
        a, b = sha(os.path.join(inst, f)), sha(os.path.join(SRC, "py4DSTEM", f))
        out["files"][f] = {"installed": a, "tree": b, "equal": a == b}
        assert a == b, f"installed py4DSTEM differs from the tree at the lock: {f}"
    t = out["files"]["process/utils/scattering_factors.txt"]["installed"]
    assert t == TABLE_SHA256, f"scattering table {t} is not the lock"
    out["table_sha256"] = t
    return out


def peaks(pl):
    d = pl.data
    return {"qx": [float(x) for x in d["qx"]], "qy": [float(x) for x in d["qy"]], "qz": [float(x) for x in d["qz"]],
            "intensity": [float(x) for x in d["intensity"]],
            "h": [int(x) for x in d["h"]], "k": [int(x) for x in d["k"]], "l": [int(x) for x in d["l"]]}


def main():
    req = json.load(sys.stdin)
    k_max, slab, width = req["k_max"], req["slab"], req["width"]
    sigma = width / math.sqrt(2.0)
    mult = slab / sigma
    out = {"versions": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                        "py4DSTEM": py4DSTEM.__version__},
           "provenance": provenance(),
           "call": {"k_max": k_max, "sigma_excitation_error": sigma, "tol_excitation_error_mult": mult,
                    "tol_intensity": 0.0, "tol_structure_factor": 1e-4, "lambda_to_zero_voltage_V": 1e30,
                    "curved_voltage_V": 200e3},
           "phases": {}}
    for ph in req["phases"]:
        pos = np.array([s["frac"] for s in ph["sites"]], dtype=np.float64)
        num = np.array([s["z"] for s in ph["sites"]], dtype=np.intp)
        occ = np.array([s["occ"] for s in ph["sites"]], dtype=np.float64)
        cr = Crystal(pos, num, np.array(ph["cell"], dtype=np.float64), occupancy=occ)
        cr.calculate_structure_factors(k_max=k_max, tol_structure_factor=1e-4)
        sf = cr.struct_factors
        rec = {"lat_real": cr.lat_real.tolist(), "lat_inv": cr.lat_inv.tolist(),
               "hkl": cr.hkl.T.astype(np.int64).tolist(), "g": cr.g_vec_all.T.tolist(),
               "F_re": sf.real.astype(np.float64).tolist(), "F_im": sf.imag.astype(np.float64).tolist(),
               "F2": cr.struct_factors_int.astype(np.float64).tolist(), "F_dtype": str(sf.dtype),
               "F2_dtype": str(cr.struct_factors_int.dtype), "axes": {}}
        for uvw in ph["axes"]:
            key = ",".join(str(x) for x in uvw)
            z = np.array(uvw)
            O = cr.parse_orientation(zone_axis_lattice=z)
            cr.setup_diffraction(1e30)
            lam0 = cr.wavelength
            p0, _ = cr.generate_diffraction_pattern(zone_axis_lattice=z, sigma_excitation_error=sigma,
                                                    tol_excitation_error_mult=mult, tol_intensity=0.0, k_max=k_max,
                                                    keep_qz=True, return_orientation_matrix=True)
            cr.setup_diffraction(200e3)
            lam200 = cr.wavelength
            p2, _ = cr.generate_diffraction_pattern(zone_axis_lattice=z, sigma_excitation_error=sigma,
                                                    tol_excitation_error_mult=mult, tol_intensity=0.0, k_max=k_max,
                                                    keep_qz=True, return_orientation_matrix=True)
            rec["axes"][key] = {"orientation": O.tolist(), "wavelength_lambda0": lam0, "wavelength_200kV": lam200,
                                "lambda0": peaks(p0), "kv200": peaks(p2)}
        out["phases"][ph["name"]] = rec
    json.dump(out, sys.stdout, allow_nan=False)


if __name__ == "__main__":
    main()
