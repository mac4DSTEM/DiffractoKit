"""Pins for RingOverlay (DiffractoGeometry), written from macTEM's RingOverlay.swift by transcription into pure Python
(math module only; no NumPy needed). The app source is the oracle (there is no Python reference): macTEM
`TEMCore/Diffraction/RingOverlay.swift`, read read-only. Nothing here is computed by the new Swift code.
Run: python3 -I -B tools/pins/ring-display/ring-overlay-pins.py
Writes Tests/DiffractoGeometryTests/Fixtures/ring-overlay-pins.json (floats as repr strings are not needed: json repr
round-trips doubles).
"""
import json, math
from pathlib import Path

def polyline(x, y, a, b, theta, count):
    st, ct = math.sin(theta), math.cos(theta)
    out = []
    for i in range(count):
        phi = 2 * math.pi * i / count
        u, v = a * math.cos(phi), b * math.sin(phi)
        out.append([x + u * st + v * ct, y + u * ct - v * st])
    return out

def distance(px, py, x, y, a, b, theta):
    dr, dc = py - y, px - x
    ct, st = math.cos(theta), math.sin(theta)
    u, v = dr * ct + dc * st, -dr * st + dc * ct
    samples = 3600
    t, best = 0.0, math.inf
    for i in range(samples):
        phi = 2 * math.pi * i / samples
        ex, ey = a * math.cos(phi) - u, b * math.sin(phi) - v
        d2 = ex * ex + ey * ey
        if d2 < best:
            best, t = d2, phi
    for _ in range(60):
        xx, yy = a * math.cos(t), b * math.sin(t)
        dx, dy = -a * math.sin(t), b * math.cos(t)
        g1 = (xx - u) * dx + (yy - v) * dy
        g2 = dx * dx + dy * dy + (xx - u) * (-xx) + (yy - v) * (-yy)
        if g2 <= 0:
            g2 = dx * dx + dy * dy
        step = g1 / g2
        t -= step
        if abs(step) < 1e-15:
            break
    return math.hypot(a * math.cos(t) - u, b * math.sin(t) - v)

rings = [
    dict(x=100.0, y=80.0, a=50.0, b=50.0, theta=0.0),
    dict(x=64.5, y=70.25, a=40.0, b=30.0, theta=0.0),
    dict(x=64.5, y=70.25, a=40.0, b=30.0, theta=0.7),
    dict(x=200.0, y=150.0, a=120.5, b=100.25, theta=2.5),
]
points = [(64.5, 70.25), (64.5 + 55.0, 70.25), (70.0, 90.0), (100.0, 80.0), (210.0, 140.0), (0.0, 0.0), (190.0, 190.0)]
cases = []
for r in rings:
    cases.append(dict(ring=r, polyline8=polyline(r['x'], r['y'], r['a'], r['b'], r['theta'], 8),
                      distances=[dict(px=px, py=py, d=distance(px, py, r['x'], r['y'], r['a'], r['b'], r['theta'])) for px, py in points]))
out = dict(source="macTEM TEMCore/Diffraction/RingOverlay.swift, transcribed to pure Python", cases=cases)
path = Path(__file__).resolve().parents[3] / 'Tests/DiffractoGeometryTests/Fixtures/ring-overlay-pins.json'
path.write_text(json.dumps(out, indent=1) + '\n')
print(len(cases), 'rings pinned')
