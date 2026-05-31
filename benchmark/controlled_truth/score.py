#!/usr/bin/env python3
"""Score a PanIsoGuard adjudication against controlled truth.

Usage: python score.py <truth.tsv> <prefix.adjudicated.tsv>
"""
import sys
from collections import Counter

truth = {}
for ln in open(sys.argv[1]):
    if ln.startswith("isoform"):
        continue
    i, t = ln.rstrip("\n").split("\t")
    truth[i] = t

cls, header = {}, None
for ln in open(sys.argv[2]):
    f = ln.rstrip("\n").split("\t")
    if header is None:
        header = {n: k for k, n in enumerate(f)}
        continue
    cls[f[header["isoform_id"]]] = f[header["confidence_class"]]

def pred(c):
    if c in ("HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL"):
        return "genuine"
    if c in ("LOW_CONF_PARTIAL", "ARTIFACT"):
        return "false"
    return "abstain"  # AMBIGUOUS / HIGH_CONF_KNOWN

conf = Counter()
for i, t in truth.items():
    conf[(t, pred(cls.get(i, "MISSING")))] += 1

print(f"{'':9}pred_genuine  pred_false  abstain")
for t in ("genuine", "false"):
    print(f"{t:9}{conf[(t,'genuine')]:11}{conf[(t,'false')]:12}{conf[(t,'abstain')]:9}")

TP, FP = conf[("genuine", "genuine")], conf[("false", "genuine")]
FN, TN = conf[("genuine", "false")], conf[("false", "false")]
dec = TP + FP + FN + TN
prec = TP / (TP + FP) if TP + FP else 0.0
spec = TN / (TN + FP) if TN + FP else 0.0
rec = TP / (TP + FN) if TP + FN else 0.0
f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
acc = (TP + TN) / dec if dec else 0.0
print(f"\nDecisive calls (n={dec}): precision(genuine)={prec:.3f} "
      f"specificity(false)={spec:.3f} F1={f1:.3f} accuracy={acc:.3f}")
ab_g, ab_f = conf[("genuine", "abstain")], conf[("false", "abstain")]
g_tot, f_tot = ab_g + TP + FN, ab_f + TN + FP
print(f"Abstain (AMBIGUOUS): genuine {ab_g}/{g_tot} ({100*ab_g/g_tot:.1f}%), "
      f"false {ab_f}/{f_tot} ({100*ab_f/f_tot:.1f}%)")
