"""Python side of the DiffractoKit benchmark: the same inputs (inputs/), the same sizes, best-of-five after one warm-up.

  ../../../refs/orix/bin/python    -I -B bench_python.py pymatgen    # pymatgen, spglib, orix
  ../../../refs/py4dstem/bin/python -I -B bench_python.py py4dstem   # py4DSTEM

One JSON line per measurement, the same fields as the kit's (side = "python"). Large sizes only with DIFFRACTO_BENCH_LARGE=1.
"""
import json, os, pathlib, sys, time
import numpy as np

here = pathlib.Path(__file__).resolve().parent
inputs = here / "inputs"
large = os.environ.get("DIFFRACTO_BENCH_LARGE") == "1"
which = sys.argv[1]


def measure(name, size, n, body, note="", backend="cpu"):
    if size == "large" and not large:
        print(json.dumps({"side": "python", "name": name, "size": size, "n": n, "backend": backend, "median_ms": None, "note": "nightly"}))
        return
    body()
    t = []
    for _ in range(5):
        t0 = time.perf_counter(); body(); t.append((time.perf_counter() - t0) * 1e3)
    t.sort()
    print(json.dumps({"side": "python", "name": name, "size": size, "n": n, "backend": backend, "median_ms": t[2],
                      "min_ms": t[0], "max_ms": t[4], "note": note}), flush=True)


if which == "pymatgen":
    import warnings; warnings.filterwarnings("ignore")
    import spglib
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.cif import CifParser
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from orix.quaternion import Orientation, Rotation, symmetry

    cells = json.load(open(inputs / "cells.json"))
    for size, n in (("small", 100), ("typical", 1000), ("large", 10000)):
        lattices = [Lattice.from_parameters(*c, vesta=True) for c in cells[:n]]
        measure("niggli_reduce", size, n, lambda: [l.get_niggli_reduced_lattice().a for l in lattices],
                note="pymatgen Lattice.get_niggli_reduced_lattice (spglib 2.8.0 niggli_reduce: see the next row)")
        mats = [np.array(l.matrix) for l in lattices]
        measure("niggli_reduce", size, n, lambda: [spglib.niggli_reduce(m) for m in mats], note="spglib.niggli_reduce", backend="spglib")

    for size, kmax in (("small", 1.0), ("typical", 2.0), ("large", 4.0)):
        al = Structure.from_file(str(inputs / "al_super1.cif"))
        rec = al.lattice.reciprocal_lattice_crystallographic
        measure("reflections", size, int(kmax * 10),
                lambda: rec.get_points_in_sphere([[0, 0, 0]], [0, 0, 0], kmax, zip_results=False),
                note=f"Lattice.reciprocal_lattice_crystallographic.get_points_in_sphere, k_max {kmax}")

    for size, nsc in (("small", 1), ("typical", 2), ("large", 3)):
        s = Structure.from_file(str(inputs / f"al_super{nsc}.cif"))
        measure("symmetry_detection", size, len(s),
                lambda: SpacegroupAnalyzer(s, symprec=1e-4, angle_tolerance=5).get_symmetry_operations(),
                note=f"SpacegroupAnalyzer(symprec 1e-4).get_symmetry_operations, supercell {nsc}^3 ({len(s)} atoms)")

    for size, f in (("small", "c2m.cif"), ("typical", "al_super2.cif"), ("large", "al_super3.cif")):
        text = (inputs / f).read_text()
        measure("cif_import", size, len(text), lambda: CifParser.from_str(text).parse_structures(primitive=False), note=f)

    rot = np.fromfile(inputs / "rotations.f64").reshape(-1, 3, 3)
    for size, n in (("small", 1000), ("typical", 10000), ("large", 100000)):
        R = rot[:n]
        def reduce():
            o = Orientation.from_matrix(R, symmetry=symmetry.Oh)
            return o.map_into_symmetry_reduced_zone().angle.sum()
        measure("orientation_reduce", size, n, reduce, note="orix Orientation.from_matrix + map_into_symmetry_reduced_zone (Oh), vectorised")

if which == "py4dstem":
    import warnings; warnings.filterwarnings("ignore")
    from py4DSTEM.process.diffraction import Crystal
    from py4DSTEM.process.calibration import fit_ellipse_1D
    from scipy.spatial.transform import Rotation as Rot

    for size, kmax in (("small", 1.0), ("typical", 2.0), ("large", 3.0)):
        from pymatgen.core import Lattice, Structure
        from pymatgen.io.cif import CifParser
        s = CifParser(str(inputs / "c2m.cif")).parse_structures(primitive=False)[0]
        cr = Crystal(s.frac_coords, np.array([sp.Z for sp in s.species]), np.array(s.lattice.parameters))
        measure("kinematic_reflections", size, int(kmax * 10),
                lambda: cr.calculate_structure_factors(k_max=kmax), note=f"py4DSTEM Crystal.calculate_structure_factors, C2/m, k_max {kmax}")

    for size, n in (("small", 128), ("typical", 256), ("large", 512)):
        img = np.fromfile(inputs / f"ring_{n}.f32", dtype=np.float32).reshape(n, n).astype(np.float64)
        c = (n - 1) / 2 + 0.3
        def fit():
            return fit_ellipse_1D(img, (c, c), (0.15 * n, 0.45 * n))
        try:
            fit()
        except Exception as e:  # signature differences are reported, not hidden
            print(json.dumps({"side": "python", "name": "ellipse_fit_1d", "size": size, "n": n, "median_ms": None, "note": f"error {e!r}"}))
            continue
        measure("ellipse_fit_1d", size, n, fit, note=f"py4DSTEM fit_ellipse_1D, {n} x {n}")

    al = Crystal([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]], [13] * 4, 4.0495)
    al.calculate_structure_factors(k_max=1.6)
    for size, step in (("small", 4.0), ("typical", 1.6), ("large", 0.85)):
        def plan():
            al.setup_diffraction(200e3)
            al.orientation_plan(zone_axis_range="auto", angle_step_zone_axis=step, angle_step_in_plane=360 / 128,
                                accel_voltage=200e3, sigma_excitation_error=0.03, progress_bar=False)
        plan()
        templates = int(al.orientation_num_zones)
        measure("acom_library", size, templates, plan,
                note=f"py4DSTEM orientation_plan, angle_step_zone_axis {step} deg = {templates} zone axes, 128 in-plane steps")
        if size == "typical":
            rng = np.random.default_rng(5)
            pls = [al.generate_diffraction_pattern(orientation_matrix=Rot.random(random_state=int(rng.integers(1 << 30))).as_matrix(),
                                                   sigma_excitation_error=0.03, tol_excitation_error_mult=3, k_max=1.6) for _ in range(8)]
            measure("acom_match_scan", "small", 8, lambda: [al.match_single_pattern(pl, verbose=False) for pl in pls],
                    note=f"py4DSTEM match_single_pattern x 8 patterns, {templates} zone axes, one core", backend="python 1 core")
            pls32 = pls * 4
            measure("acom_match_scan", "typical", 32, lambda: [al.match_single_pattern(pl, verbose=False) for pl in pls32],
                    note=f"py4DSTEM match_single_pattern x 32 patterns, {templates} zone axes, one core (a 32 x 32 scan is 32x this)",
                    backend="python 1 core")
