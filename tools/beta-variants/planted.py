#!/usr/bin/env python3
"""Planted truth for DiffractoKit task 35: the 12 orientation-relationship variants of beta'' (Mg5Si6, C2/m) in Al.

numpy only. Writes `variants.npy` (12, 3, 3), `measured.npy` (N, 3, 3) and `truth.json` into --out (default ./out beside this
script; do not commit the output). It runs no orix and claims no orix number: the Mac pin script (kit audit step 0) adds
`Misorientation.equivalent()` / `reduce()` at atol 1e-4 from the same rotation matrices.

Frames and conventions (rows and columns; matrices act on column vectors):
  * Al frame: Cartesian axes along the cubic axes [100], [010], [001]. Matrix proper point group = the 24 proper rotations
    of m-3m (432). Precipitate proper point group of 2/m = {E, C2(b)}: inversion is dropped because rotations only.
  * Precipitate Cartesian frame (monoclinic, unique axis b): e1 || a, e2 || b, e3 = e1 x e2 (|| c*). Convert in the wrapper if
    the kit uses another setting; the lattice-direction statement below does not depend on it.
  * R0 maps precipitate Cartesian vectors into the Al frame. The orientation relationship (Andersen et al. 1998, Acta Mater.
    46, 3283; Zandbergen et al. 1997, Science 277, 1221) as derived here from the lattice parameters:
        b_beta'' || [001]Al                      (needle axis, D24)
        a_beta'' || [3 2 0]Al                    (|[320]| a_Al = 14.60 A against a = 15.16 A, misfit 3.8 %)
        c_beta'' || [1 -3 0]Al, in the (001)Al plane
    The angle between [3 2 0] and [1 -3 0] is arccos(-3/sqrt(130)) = 105.25 degrees, against beta = 105.3 degrees, which is why
    this pair was chosen; right-handedness (a x b . c > 0) is asserted below. The script's R0 is therefore the OR the planted
    case uses; the study's own OR (scope question 1) replaces it by editing R0 only. VERIFY against the paper before
    quoting it: this derivation, not a table lookup, produced the [3 2 0] and [1 -3 0] pair.
  * Variants: R_k = O R0 over the 24 cubic rotations O, with R_i ~ R_j when R_i = R_j S for S in the precipitate group.
    Representatives are the first of each class in the fixed order of `cubic_rotations()`.

Truth carried in truth.json:
  * 24 / |{E, R0 C2 R0^T}| = 12 variants, because R0 C2(b) R0^T is the 180 degree rotation about [001]Al, a cubic operation.
  * needle axis (b in the Al frame) per variant: +-[100], +-[010], +-[001], four variants for each of the three axes.
  * pairwise variant-to-variant disorientation in the precipitate frame, min over S1, S2 in {E, C2(b)} of the angle of
    S1 R_i^T R_j S2, angle = arccos((trace - 1) / 2): only 0, 90, 120 and 180 degrees occur; counts per angle for the 132
    ordered pairs i != j are in the file.
  * the OR's own reduced disorientation angle (cubic on the left, 2/m on the right): arctan(2/3) = 33.690068 degrees.
  * N planted measured orientations g = R_k S delta (S random in the group, delta a seeded 0.5 degree rotation) with their
    true nearest variant k: because distinct variants are at least 90 degrees apart, the nearest by disorientation is
    unambiguous (margin > 89 degrees), so `nearestVariant` must return k exactly.
Tolerance rationale: matrices are built from sqrt(13) and 0/1 entries, so every number is exact to double round-off (1e-15);
compare variants set-wise (as classes under the right group) at atol 1e-9, angles at atol 1e-9 degrees. orix pins sit at
atol 1e-4 (test_misorientation.py), which the Mac pin step applies, not this file.

Run:  python3 -I planted.py [--out DIR]
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np

A_AL = 4.05
LAT = dict(a=15.16, b=4.05, c=6.74, beta=105.3)   # CrystalPresets.betaDoublePrime
C2_B = np.diag([-1.0, 1.0, -1.0])                  # C2 about e2 (b) in the precipitate Cartesian frame
SYM_P = [np.eye(3), C2_B]


def cubic_rotations():
    ops = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([1, -1], repeat=3):
            m = np.zeros((3, 3))
            for i, p in enumerate(perm):
                m[i, p] = signs[i]
            if np.linalg.det(m) > 0:
                ops.append(m)
    return ops


def angle_deg(m):
    return float(np.degrees(np.arccos(np.clip((np.trace(m) - 1) / 2, -1, 1))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "out"))
    ap.add_argument("--n", type=int, default=48, help="number of planted measured orientations")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    s13 = np.sqrt(13.0)
    ea = np.array([3.0, 2.0, 0.0]) / s13           # a || [3 2 0]Al
    eb = np.array([0.0, 0.0, 1.0])                  # b || [001]Al
    e3 = np.cross(ea, eb)                           # c* direction
    r0 = np.column_stack([ea, eb, e3])
    assert abs(np.linalg.det(r0) - 1) < 1e-14 and np.allclose(r0.T @ r0, np.eye(3), atol=1e-14)
    c_dir = np.array([1.0, -3.0, 0.0])
    ang_ac = np.degrees(np.arccos(ea @ c_dir / np.linalg.norm(c_dir)))
    assert abs(ang_ac - LAT["beta"]) < 0.1, ang_ac
    assert np.dot(np.cross(ea, eb), c_dir) > 0      # right-handed (a, b, c)

    ops = cubic_rotations()
    assert len(ops) == 24
    reps = []
    for o in ops:
        r = o @ r0
        if not any(np.allclose(r, rj @ s, atol=1e-12) for rj in reps for s in SYM_P):
            reps.append(r)
    assert len(reps) == 12, len(reps)
    variants = np.stack(reps)
    np.save(os.path.join(args.out, "variants.npy"), variants)

    needle = []
    for r in reps:
        b = np.round(r[:, 1]).astype(int)
        needle.append(b.tolist())
    axes = {}
    for b in needle:
        axes.setdefault(str(int(np.argmax(np.abs(b)))), []).append(b)
    assert all(len(v) == 4 for v in axes.values()) and len(axes) == 3

    ang = np.zeros((12, 12))
    for i in range(12):
        for j in range(12):
            ang[i, j] = min(angle_deg(s1 @ reps[i].T @ reps[j] @ s2) for s1 in SYM_P for s2 in SYM_P)
    values, counts = np.unique(np.round(ang[~np.eye(12, dtype=bool)], 9), return_counts=True)
    assert set(np.round(values, 6)) <= {90.0, 120.0, 180.0}, values
    assert np.allclose(np.diag(ang), 0)

    or_angle = min(angle_deg(o @ r0 @ s) for o in ops for s in SYM_P)
    assert abs(or_angle - np.degrees(np.arctan2(2, 3))) < 1e-9

    rng = np.random.default_rng(35)
    measured, truth_k = [], []
    for n in range(args.n):
        k = int(rng.integers(12))
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        th = np.radians(0.5)
        kx = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        delta = np.eye(3) + np.sin(th) * kx + (1 - np.cos(th)) * kx @ kx   # Rodrigues
        g = reps[k] @ SYM_P[int(rng.integers(2))] @ delta
        # nearest variant by disorientation angle, as the kit must compute it
        d = [min(angle_deg(rj.T @ g @ s) for s in SYM_P) for rj in reps]
        order = np.argsort(d)
        assert order[0] == k and d[order[1]] > 89
        measured.append(g)
        truth_k.append(k)
    np.save(os.path.join(args.out, "measured.npy"), np.stack(measured))

    truth = {
        "generator": "DiffractoKit/tools/beta-variants/planted.py",
        "upstream_pin": "none; no orix or diffsims run here. The Mac pin script adds Misorientation.equivalent()/reduce() at atol 1e-4.",
        "phase": {"preset": "CrystalPresets.betaDoublePrime", "space_group": "C2/m", "lattice": LAT},
        "orientation_relationship": {
            "b_beta||": "[001]Al", "a_beta||": "[3 2 0]Al", "c_beta||": "[1 -3 0]Al",
            "angle_between_320_and_1-30_deg": float(ang_ac), "beta_deg": LAT["beta"],
            "a_misfit": float(np.linalg.norm([3, 2, 0]) * A_AL / LAT["a"] - 1),
            "source": "derived from the lattice parameters; check against Andersen 1998 before quoting",
            "precipitate_cartesian_frame": "e1||a, e2||b, e3=e1xe2; R0 columns = e1, e2, e3 in the Al frame",
            "R0": r0.tolist()},
        "variant_count": 12,
        "variant_count_reason": "24 proper cubic rotations / |{E, 180 deg about [001]Al}| = 12; R0 C2(b) R0^T is 180 deg about [001]Al",
        "variants_file": "variants.npy (12, 3, 3), R_k maps precipitate Cartesian to the Al frame",
        "needle_axis_by_variant_al_frame": needle,
        "needle_variants_per_cubic_axis": {k: len(v) for k, v in sorted(axes.items())},
        "pairwise_disorientation_deg": ang.tolist(),
        "pairwise_angle_counts_ordered_pairs": {str(float(v)): int(c) for v, c in zip(values, counts)},
        "or_reduced_disorientation_deg": float(or_angle),
        "or_reduced_disorientation_closed_form": "arctan(2/3) = 33.690068 degrees",
        "measured": {"file": "measured.npy", "count": args.n, "true_variant": truth_k,
                     "noise": "seeded 0.5 degree rotation, S random in {E, C2(b)}, matrix orientation = identity",
                     "min_margin_to_second_variant_deg": 89},
        "tolerance_rationale": __doc__.split("Tolerance rationale:")[1].split("Run:")[0].strip(),
    }
    with open(os.path.join(args.out, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=2)
    print(f"12 variants; needles per axis {truth['needle_variants_per_cubic_axis']}; "
          f"ordered-pair angle counts {truth['pairwise_angle_counts_ordered_pairs']}; OR reduced angle {or_angle:.6f} deg; "
          f"angle(a,c)={ang_ac:.3f} deg; a misfit {truth['orientation_relationship']['a_misfit']*100:.1f} %")
    return 0


if __name__ == "__main__":
    sys.exit(main())
