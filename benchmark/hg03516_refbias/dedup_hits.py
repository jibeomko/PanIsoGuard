#!/usr/bin/env python3
"""Collect CREATED (reference-bias) junctions from the per-chrom/hap scans, dedup across
pat/mat, and write refbias_unique.tsv (chrom, intron_start, intron_end, strand, haplotypes).

Usage: dedup_hits.py <work_dir>
"""
import sys, glob, os
from collections import defaultdict

W = sys.argv[1]
hits = defaultdict(set)   # (chrom, a, b, strand) -> {pat, mat}
for gtf in glob.glob(f"{W}/scan/chr*/vt.gtf"):
    base = gtf.split('/')[-2]            # e.g. chr1.pat
    chrom, hap = base.split('.')[0], base.split('.')[-1]
    truth = {l.split('\t')[0]: l.strip().split('\t')[1] for l in open(f"{W}/scan/{base}/vt.truth")}
    ex = defaultdict(list); strand = {}
    for ln in open(gtf):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        iso = f[8].split('transcript_id "')[1].split('"')[0]
        ex[iso].append((int(f[3]), int(f[4]))); strand[iso] = f[6]
    for iso, exs in ex.items():
        if truth.get(iso) != "CREATED":
            continue
        e = sorted(exs)
        a, b = e[0][1], e[1][0] - 1       # intron [a, b]
        hits[(chrom, a, b, strand[iso])].add(hap)

with open(f"{W}/refbias_unique.tsv", "w") as o:
    o.write("chrom\tintron_start\tintron_end\tstrand\thaplotypes\n")
    for (c, a, b, s), haps in sorted(hits.items()):
        o.write(f"{c}\t{a}\t{b}\t{s}\t{','.join(sorted(haps))}\n")
print(f"wrote refbias_unique.tsv: {len(hits)} unique reference-bias junctions")
