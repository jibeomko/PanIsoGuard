#!/usr/bin/env python3
"""BED of the regions the variant axis can read in a run: the span of every isoform that has a
junction outside the annotation, +- 1 kb, merged. From a run's junction table and the isoform GTF.

Usage: novel_isoform_regions.py JUNCTIONS.tsv ISOFORMS.gtf [JUNCTIONS.tsv ISOFORMS.gtf ...] > regions.bed
"""
import collections, csv, sys

PAD = 1000
by = collections.defaultdict(list)
for jx, gtf in zip(sys.argv[1::2], sys.argv[2::2]):
    span = {}
    for line in open(gtf):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "transcript":
            span[f[8].split('transcript_id "')[1].split('"')[0]] = (f[0], int(f[3]) - 1, int(f[4]))
    for r in csv.DictReader(open(jx), delimiter="\t"):
        if r["known"] != "false": continue
        for t in r["isoforms"].split(","):
            if t in span: c, a, b = span[t]; by[c].append((max(0, a - PAD), b + PAD))
for c in sorted(by):
    v = sorted(by[c]); cur = list(v[0])
    for a, b in v[1:]:
        if a <= cur[1]: cur[1] = max(cur[1], b)
        else: print(c, cur[0], cur[1], sep="\t"); cur = [a, b]
    print(c, cur[0], cur[1], sep="\t")
