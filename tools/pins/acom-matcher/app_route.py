"""A NumPy TRANSCRIPTION of the mac4DSTEM ACOM polar-correlation route, as the app source is (never generated from the Swift in
this kit).

Source: mac4DSTEM commit c809679c (the files below were last changed in deb2ec4f), read with `git show`/from the working tree:

    mac4DSTEM/Core/Crystal/OrientationPlan.swift   (project, buildPolar, blurAzimuthal, normalizeUnit, ringFFTs, generate)
    mac4DSTEM/Core/Crystal/OrientationMatcher.swift (prepareExperimentalFFT, match, selectOrientation, makeOrientationResult)
    mac4DSTEM/Core/Analysis/OrientationResult.swift (ACOMOrientation.matrix, the cubic / hexagonal / monoclinic fundamental-zone
                                                    samplers, the cubic reduction)

Arithmetic follows the Swift: polar images and FFTs in Float32 (NumPy float32; the FFT is evaluated in float64 and rounded to
float32, vDSP's differs by about 1e-7 relative), deposition weights in Double, rounding half away from zero (Swift
`.rounded()`), strict `>` ties. Sums run in the Swift's order where it matters; where NumPy sums pairwise (mean, norm) the
difference is a few float32 ulps, which is why the Swift tests state tolerances rather than bit equality for the images.

Only the numbers: no Metal, no scan selection, no cancellation, no progress.
"""
import math

import numpy as np

F32 = np.float32


def rnd(x):
    """Swift `Double.rounded()`: half away from zero."""
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


# ---------------------------------------------------------------------------------------------------------------------
# detector basis and orientation matrix (ZoneAxis.detectorBasis / ACOMOrientation.matrix)

def detector_basis(n):
    n = np.asarray(n, dtype=np.float64)
    ref = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, ref)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return np.stack([e1, e2, n], axis=1)  # columns


def orientation_matrix(basis, angle32):
    a = float(F32(angle32))
    c, s = math.cos(a), math.sin(a)
    inplane = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1.0]])  # columns (c,-s,0), (s,c,0), (0,0,1)
    return basis @ inplane


# ---------------------------------------------------------------------------------------------------------------------
# fundamental-zone samplers (OrientationResult.swift; chord = 1 - dot, first maximum wins)

def _fps(count, vertices, candidates, chord):
    vertices = [np.array(v, dtype=np.float64) for v in vertices]
    if count <= len(vertices):
        return vertices[:count]
    selected = list(vertices)
    minc = np.array([min(chord(c, v) for v in vertices) for c in candidates])
    while len(selected) < count:
        best = int(np.argmax(minc))  # Swift `indices.max(by: <)` keeps the first of equal maxima, as argmax
        nxt = candidates[best]
        selected.append(nxt)
        for i in range(len(candidates)):
            minc[i] = min(minc[i], chord(candidates[i], nxt))
        minc[best] = -np.inf
    return selected


def _fib(n):
    golden = math.pi * (3 - math.sqrt(5))
    out = []
    for i in range(n):
        z = 1 - 2 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1 - z * z))
        out.append((r * math.cos(golden * i), r * math.sin(golden * i), z))
    return out


def _cubic_fold(d):
    v = sorted(abs(x) for x in d)
    w = np.array([v[1], v[0], v[2]])
    return w / np.linalg.norm(w)


def _hex_fold(d):
    d = np.asarray(d, dtype=np.float64)
    v = d / np.linalg.norm(d)
    if v[2] < 0:
        v = -v
    basal = math.hypot(v[0], v[1])
    az = math.fmod(math.atan2(v[1], v[0]), math.pi / 3)
    if az < 0:
        az += math.pi / 3
    if az > math.pi / 6:
        az = math.pi / 3 - az
    w = np.array([basal * math.cos(az), basal * math.sin(az), v[2]])
    return w / np.linalg.norm(w)


def sample_cubic(count):
    if count <= 0:
        return []
    verts = [(0, 0, 1), tuple(np.array([1, 0, 1]) / math.sqrt(2)), tuple(np.array([1, 1, 1]) / math.sqrt(3))]
    if count <= 3:
        return [np.array(v) for v in verts[:count]]
    n = max(128, count * 12)
    cands = [_cubic_fold(c) for c in _fib(n)]
    return _fps(count, verts, cands, lambda a, b: 1 - float(np.dot(a, b)))


def sample_hex(count):
    if count <= 0:
        return []
    verts = [(0, 0, 1), (1, 0, 0), (math.cos(math.pi / 6), math.sin(math.pi / 6), 0)]
    if count <= 3:
        return [np.array(v, dtype=np.float64) for v in verts[:count]]
    n = max(192, count * 16)
    cands = [_hex_fold(c) for c in _fib(n)]
    return _fps(count, verts, cands, lambda a, b: 1 - float(np.dot(a, b)))


def _mono_fold(d):
    d = np.asarray(d, dtype=np.float64)
    v = d / np.linalg.norm(d)
    if v[1] < 0:
        v = -v
    basal = math.hypot(v[0], v[2])
    az = math.fmod(math.atan2(v[2], v[0]), math.pi)
    if az < 0:
        az += math.pi
    if az > math.pi - 1e-15:
        az = 0.0
    w = np.array([basal * math.cos(az), v[1], basal * math.sin(az)])
    return w / np.linalg.norm(w)


def sample_mono(count):
    if count <= 0:
        return []
    verts = [(0, 1, 0), (1, 0, 0), (0, 0, 1)]
    if count <= 3:
        return [np.array(v, dtype=np.float64) for v in verts[:count]]
    n = max(192, count * 16)
    golden = math.pi * (3 - math.sqrt(5))
    cands = []
    for i in range(n):
        y = 1 - 2 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1 - y * y))
        cands.append(_mono_fold((r * math.cos(golden * i), y, r * math.sin(golden * i))))

    def chord(a, b):
        return 1 - max(float(np.dot(a, b)), float(np.dot(a, np.array([-b[0], b[1], -b[2]]))))
    return _fps(count, verts, cands, chord)


# ---------------------------------------------------------------------------------------------------------------------
# polar image (OrientationPlan)

def project(g, inten, n, sg_width, sg_max, wavelength, power, radial_power=0.0, py4dstem=False):
    """OrientationPlan.project: g is (N,3) float64, inten (N,). Returns lists r, azim, weight (float64) in input order.

    py4DSTEM=True is NOT the app: it replaces the app's sg = g.n + lambda |g|^2 / 2 by py4DSTEM's Crystal.excitation_errors,
    sg = (2 g.n - lambda |g|^2) / (2 - 2 lambda g.n) (the kit's own option, for the planted-pattern comparison)."""
    basis = detector_basis(n)
    e1, e2 = basis[:, 0], basis[:, 1]
    curvature = (wavelength or 0.0) / 2
    r_out, a_out, w_out = [], [], []
    for gi, ii in zip(g, inten):
        if py4dstem:
            gz = float(np.dot(gi, n))
            sg = (2 * gz - wavelength * float(np.dot(gi, gi))) / (2 - 2 * wavelength * gz)
        else:
            sg = float(np.dot(gi, n)) + curvature * float(np.dot(gi, gi))
        if abs(sg) > sg_max:
            continue
        excited = ii * math.exp(-(sg / sg_width) * (sg / sg_width))
        w = excited if power == 1 else math.pow(excited, power)
        x = float(np.dot(gi, e1))
        y = float(np.dot(gi, e2))
        radius = math.sqrt(x * x + y * y)
        if radial_power != 0:
            w *= math.pow(radius, radial_power)
        r_out.append(radius)
        a_out.append(math.atan2(y, x))
        w_out.append(w)
    return r_out, a_out, w_out


def blur_azimuthal(img, nr, na, sigma):
    radius = max(1, rnd(3 * sigma))
    taps = np.zeros(2 * radius + 1, dtype=F32)
    s = F32(0)
    for i in range(-radius, radius + 1):
        v = np.exp(-F32(i * i) / F32(2 * sigma * sigma)).astype(F32)
        taps[i + radius] = v
        s = F32(s + v)
    taps = (taps / s).astype(F32)
    img2 = img.reshape(nr, na)
    out = np.zeros_like(img2)
    for t in range(-radius, radius + 1):
        out = (out + np.roll(img2, -t, axis=1) * taps[t + radius]).astype(F32)  # acc over t ascending; img[(a+t)%na]
    return out.reshape(-1)


def build_polar(r, azim, weight, nr, na, radial_scale, azim_blur, radial_kernel):
    img = np.zeros(nr * na, dtype=F32)
    two_pi = 2 * math.pi
    sigma_bins = radial_kernel / radial_scale
    reach = max(1, rnd(3 * sigma_bins)) if sigma_bins > 0 else 0
    for rr, aa, ww in zip(r, azim, weight):
        exact = rr / radial_scale
        rb = rnd(exact)
        if rb + reach < 0 or rb - reach >= nr:
            continue
        a = math.fmod(aa, two_pi)
        if a < 0:
            a += two_pi
        ab = rnd(a / two_pi * na) % na
        if reach == 0:
            if 0 <= rb < nr:
                img[rb * na + ab] = F32(img[rb * na + ab] + F32(ww))
            continue
        for b in range(rb - reach, rb + reach + 1):
            if 0 <= b < nr:
                d = (b - exact) / sigma_bins
                img[b * na + ab] = F32(img[b * na + ab] + F32(ww * math.exp(-0.5 * d * d)))
    if azim_blur > 0:
        img = blur_azimuthal(img, nr, na, azim_blur)
    img2 = img.reshape(nr, na)
    for ri in range(nr):
        m = F32(0)
        row = img2[ri]
        for v in row:  # sequential float32 sum, as the Swift loop
            m = F32(m + v)
        m = F32(m / F32(na))
        img2[ri] = (row - m).astype(F32)
    return img2.reshape(-1)


def normalize_unit(img):
    norm = F32(0)
    for v in img:
        norm = F32(norm + F32(v * v))
    norm = F32(np.sqrt(norm))
    if norm > 0:
        img = (img / norm).astype(F32)
    return img


def ring_ffts(polar, nr, na):
    f = np.fft.fft(polar.reshape(nr, na).astype(np.float64), axis=1)
    return f.real.astype(F32), f.imag.astype(F32)


# ---------------------------------------------------------------------------------------------------------------------
# plan

def generate_plan(g, inten, k_max, axes, *, nr=32, na=128, sg_width=0.03, sg_max=0.1, azim_blur=1.5, wavelength=None,
                  power=0.25, radial_power=0.0, radial_kernel=0.08, py4dstem=False):
    scale = k_max / nr
    plan = {"scale": scale, "nr": nr, "na": na, "templates": [], "spots": [], "fre": [], "fim": [], "power": power,
            "kernel": radial_kernel}
    for n in axes:
        r, a, w = project(g, inten, n, sg_width, sg_max, wavelength, power, radial_power, py4dstem)
        pol = normalize_unit(build_polar(r, a, w, nr, na, scale, azim_blur, radial_kernel))
        fre, fim = ring_ffts(pol, nr, na)
        plan["templates"].append(pol)
        plan["spots"].append((r, a, w))
        plan["fre"].append(fre)
        plan["fim"].append(fim)
    return plan


# ---------------------------------------------------------------------------------------------------------------------
# matcher

def experimental_polar(plan, peaks, origin, inv_per_pixel):
    """peaks: (N,3) float32 x, y, intensity; origin float32 pair."""
    r, a, w = [], [], []
    for x, y, i in peaks:
        dx = float(F32(x) - F32(origin[0]))
        dy = float(F32(y) - F32(origin[1]))
        rr = math.sqrt(dx * dx + dy * dy) * inv_per_pixel
        inten = float(max(F32(i), F32(0)))
        r.append(rr)
        a.append(math.atan2(dy, dx))
        w.append(inten if plan["power"] == 1 else math.pow(inten, plan["power"]))
    pol = build_polar(r, a, w, plan["nr"], plan["na"], plan["scale"], 1.5, plan["kernel"])
    return normalize_unit(pol)


def template_correlations(plan, polar):
    """Per template: (best score pass0, bin0, best score pass1, bin1). The experimental ring FFT E; pass 0 uses conj(E)*T,
    pass 1 uses E*T; the inverse transform is vDSP's unnormalised one scaled by 1/na, i.e. NumPy's ifft."""
    nr, na = plan["nr"], plan["na"]
    ere, eim = ring_ffts(polar, nr, na)
    E = ere.astype(np.float64) + 1j * eim.astype(np.float64)
    out = []
    for fre, fim in zip(plan["fre"], plan["fim"]):
        T = fre.astype(np.float64) + 1j * fim.astype(np.float64)
        res = []
        for conj in (True, False):
            prod = (np.conj(E) if conj else E) * T
            corr = np.fft.ifft(prod.sum(axis=0)).real
            b = int(np.argmax(corr))  # first maximum (strict > in the Swift)
            res.append((F32(corr[b]), b))
        out.append(res)
    return out, (ere, eim)


def select_orientation(axes, scores, bins, distinct_rad):
    best, bt = -np.inf, -1
    for t, s in enumerate(scores):
        if s > best:
            best, bt = s, t
    if bt < 0:
        return -1, 0, 0.0, 0.0
    cos_limit = math.cos(distinct_rad)
    winner = axes[bt]
    second = -np.inf
    for t, s in enumerate(scores):
        if t == bt:
            continue
        if float(np.dot(axes[t], winner)) > cos_limit:
            continue
        if s > second:
            second = s
    if second == -np.inf:
        second = 0.0
    return bt, bins[bt], best, second


def match(plan, axes, bases, peaks, origin, inv_per_pixel, distinct_deg=10.0):
    if len(peaks) < 2:
        return None
    polar = experimental_polar(plan, peaks, origin, inv_per_pixel)
    corr, (ere, eim) = template_correlations(plan, polar)
    scores, bins, mirrored = [], [], []
    for (s0, b0), (s1, b1) in corr:
        if s1 > s0:
            scores.append(s1); bins.append(b1); mirrored.append(True)
        else:
            scores.append(s0); bins.append(b0); mirrored.append(False)
    t, b, best, second = select_orientation(axes, scores, bins, distinct_deg * math.pi / 180)
    na = plan["na"]
    mir = mirrored[t]
    bin_ = (na - b) % na
    angle = F32(bin_ / na * 2 * math.pi)
    if mir:
        angle = F32(math.fmod(float(F32(angle + F32(math.pi))), float(F32(2 * F32(math.pi)))))
    M = orientation_matrix(bases[t], angle)
    if mir:
        M[:, 1] = -M[:, 1]
        M[:, 2] = -M[:, 2]
    return {"template": t, "bin": b, "score": float(max(best, 0)), "second": float(max(second, 0)),
            "mirrored": mir, "angle": float(angle), "matrix": M, "scores": np.array(scores, dtype=F32),
            "bins": np.array(bins), "mirror_flags": np.array(mirrored), "polar": polar, "ere": ere, "eim": eim}


# ---------------------------------------------------------------------------------------------------------------------
# cubic reduction (CubicOrientationSymmetry.reduce): minimum rotation angle over the 24 proper operators

def cubic_operators():
    ops = []
    import itertools
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product([-1.0, 1.0], repeat=3):
            m = np.zeros((3, 3))
            for row in range(3):
                m[row, perm[row]] = signs[row]
            if np.linalg.det(m) > 0.5:
                ops.append(m)
    return ops


def reduce_min_angle(M, ops):
    best, ba = None, None
    for s in ops:
        c = s @ M
        w = (np.trace(c) - 1) / 2
        ang = math.acos(max(-1.0, min(1.0, w)))
        if ba is None or ang < ba - 1e-12:
            best, ba = c, ang
    return best, ba
