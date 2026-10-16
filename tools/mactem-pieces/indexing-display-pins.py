"""Pins for IndexingDisplay (DiffractoGeometry): the numbers a room shows for an accepted ring, from macTEM's
IndexingDisplay.swift (the app source is the oracle; no Python reference), transcribed into pure Python (math only). The
decision (which family is accepted, the refusal) is the kit's `RingIndexing`, pinned elsewhere (P9); here the family list is
the dyadic ring list g = 1, 1.5, 2, 3 (exact doubles) and the accepted family is stated by hand.
    sigma_q = q*tol/100;  z = |q-g|/sigma_q;  delta% = 100*(g/q-1);  |dd|/d = |g/q-1|;  margin = (z2^2 - z1^2)/2 with the
    next family the second-smallest z;  odds = exp(margin) (inf on overflow);  text '%.3g : 1' or '> 10^6 : 1' above 1e6.
Run: python3 -I -B tools/mactem-pieces/indexing-display-pins.py
Writes Tests/DiffractoGeometryTests/Fixtures/indexing-display-pins.json
"""
import json, math
from pathlib import Path

gs = [1.0, 1.5, 2.0, 3.0]
cases = []
for q, tol in [(1.1, 10.0), (1.01, 1.0), (1.45, 2.0), (2.0, 0.5), (2.9, 5.0), (1.1, 12.0)]:
    sigma = q * tol / 100
    zs = sorted((abs(q - g) / sigma, g) for g in gs)
    (z1, g1), (z2, g2) = zs[0], zs[1]
    margin = (z2 * z2 - z1 * z1) / 2
    try:
        odds = math.exp(margin)
    except OverflowError:
        odds = math.inf
    text = '> 10⁶ : 1' if odds > 1e6 else '%.3g : 1' % odds
    cases.append(dict(q=q, tolerancePercent=tol, sigmaQ=sigma, g=g1, z=z1, deltaPercent=100 * (g1 / q - 1),
                      relativeDeltaD=abs(g1 / q - 1), dMeasured=1 / q, dListed=1 / g1, nextG=g2, margin=margin,
                      odds=None if math.isinf(odds) else odds, oddsInfinite=math.isinf(odds), oddsText=text))
out = dict(source="macTEM TEMCore/Diffraction/IndexingDisplay.swift, transcribed to pure Python", families=gs, cases=cases)
path = Path(__file__).resolve().parents[2] / 'Tests/DiffractoGeometryTests/Fixtures/indexing-display-pins.json'
path.write_text(json.dumps(out, indent=1) + '\n')
for c in cases:
    print(c['q'], c['tolerancePercent'], c['g'], round(c['z'], 4), round(c['margin'], 4), c['oddsText'])
