#!/usr/bin/env python3
"""Confidence-class calibration (ECE / Brier) for a labelled adjudication.

PanIsoGuard emits ordinal confidence CLASSES, not native probabilities. This maps
each class to an empirical P(genuine) learned on a held-out *calibration* half of a
labelled set, then reports the Brier score and Expected Calibration Error (ECE) on
the *test* half (so calibration is not evaluated on the same data that defined it).

Usage: python calibrate.py <truth.tsv: isoform<TAB>genuine|false> <prefix.adjudicated.tsv>
"""
import sys, hashlib, math
from collections import defaultdict

truth = {}
for ln in open(sys.argv[1]):
    if ln.startswith("isoform"):
        continue
    p = ln.rstrip("\n").split("\t")
    truth[p[0]] = p[1]

cls, h = {}, None
for ln in open(sys.argv[2]):
    f = ln.rstrip("\n").split("\t")
    if h is None:
        h = {n: k for k, n in enumerate(f)}
        continue
    cls[f[h["isoform_id"]]] = f[h["confidence_class"]]

items = [(i, cls.get(i, "NA"), 1 if truth[i] == "genuine" else 0) for i in truth if i in cls]
half = lambda i: int(hashlib.md5(i.encode()).hexdigest(), 16) % 2
calib = [x for x in items if half(x[0]) == 0]
test = [x for x in items if half(x[0]) == 1]

num, den = defaultdict(int), defaultdict(int)
for _, c, y in calib:
    den[c] += 1; num[c] += y
prob = {c: (num[c] / den[c] if den[c] else 0.5) for c in den}

print("class -> P(genuine) [learned on calibration half]:")
for c in sorted(prob):
    print(f"  {c:28} {prob[c]:.3f}  (n={den[c]})")

brier, n = 0.0, 0
binp, binc, biny = defaultdict(float), defaultdict(int), defaultdict(int)
for _, c, y in test:
    p = prob.get(c, 0.5)
    brier += (p - y) ** 2; n += 1
    b = min(9, int(p * 10))
    binp[b] += p; binc[b] += 1; biny[b] += y
brier /= max(1, n)
ece = sum(abs(binp[b] / binc[b] - biny[b] / binc[b]) * binc[b] for b in binc) / max(1, n)
print(f"\ntest n={n}   Brier={brier:.4f}   ECE={ece:.4f}")
