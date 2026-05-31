#!/usr/bin/env python3
"""Per-category confidence-class table for a synthetic-axes adjudication.
Usage: python score_axes.py <truth.tsv> <prefix.adjudicated.tsv>
"""
import sys
from collections import defaultdict, Counter

cat = {}
for ln in open(sys.argv[1]):
    if ln.startswith("isoform"):
        continue
    i, c, t = ln.rstrip("\n").split("\t")
    cat[i] = c

cls = {}
h = None
for ln in open(sys.argv[2]):
    f = ln.rstrip("\n").split("\t")
    if h is None:
        h = {n: k for k, n in enumerate(f)}
        continue
    cls[f[h["isoform_id"]]] = f[h["confidence_class"]]

table = defaultdict(Counter)
for i, c in cat.items():
    table[c][cls.get(i, "MISSING")] += 1

order = ["HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL", "LOW_CONF_PARTIAL",
         "PAN_REF_RESCUED_FALSE_NOVEL", "ARTIFACT", "AMBIGUOUS"]
labeldesc = {"A": "A genuine", "B": "B noncanon-artifact",
             "C": "C mapping-artifact", "D": "D reference-bias"}
print(f"{'category':22} " + "  ".join(k.replace('_CONF','').replace('_FALSE_NOVEL','') for k in order))
for c in ["A", "B", "C", "D"]:
    row = "  ".join(f"{table[c].get(k,0):>{max(3,len(k.replace('_CONF','').replace('_FALSE_NOVEL','')))}}" for k in order)
    print(f"{labeldesc[c]:22} {row}")
