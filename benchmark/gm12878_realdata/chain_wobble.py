#!/usr/bin/env python3
"""Whole-chain read support on real ONT direct-RNA: exact coordinates vs a per-splice-site
tolerance. For every caller-novel chain, count the primary reads at its first donor site
whose intron list contains the chain contiguously, with every splice site within +-k bp
(k = 0 is what PanIsoGuard counts), and the offset distribution of the sites among reads
that match within 10 bp.

Usage: chain_wobble.py <novel_chains.gtf> <ont.bam> <pig_ncallers.json> --out stats.json
"""
import collections, json, statistics, sys
import pysam

gtf, bam_path, ncj = sys.argv[1], sys.argv[2], sys.argv[3]
out_path = sys.argv[sys.argv.index("--out") + 1]
nc = json.load(open(ncj))
ex = collections.defaultdict(list); chrom = {}
for l in open(gtf):
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    t = f[8].split('transcript_id "')[1].split('"')[0]; ex[t].append((int(f[3]) - 1, int(f[4]))); chrom[t] = f[0]
chains = {}
for t, e in ex.items():
    e = sorted(e)
    if len(e) > 1: chains[t] = tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))


def read_introns(r):
    out = []; pos = r.reference_start
    for op, n in r.cigartuples:
        if op == 3: out.append((pos, pos + n)); pos += n
        elif op in (0, 2, 7, 8): pos += n
    return out


def match_at(chain, ri, w):
    k = len(chain)
    for i in range(len(ri) - k + 1):
        if all(abs(ri[i + j][0] - chain[j][0]) <= w and abs(ri[i + j][1] - chain[j][1]) <= w for j in range(k)):
            return i
    return None


bam = pysam.AlignmentFile(bam_path)
W = (0, 1, 2, 3, 5, 10)
counts = {w: {} for w in W}; offsets = collections.Counter(); spanning = {}
for t, c in chains.items():
    reads = [read_introns(r) for r in bam.fetch(chrom[t], max(0, c[0][0] - 1), c[0][0] + 1)
             if not (r.is_secondary or r.is_supplementary or r.is_unmapped)]
    reads = [ri for ri in reads if ri]
    spanning[t] = len(reads)
    for w in W: counts[w][t] = sum(match_at(c, ri, w) is not None for ri in reads)
    for ri in reads:
        i = match_at(c, ri, 10)
        if i is None: continue
        for j in range(len(c)):
            offsets[abs(ri[i + j][0] - c[j][0])] += 1; offsets[abs(ri[i + j][1] - c[j][1])] += 1

tot = sum(offsets.values())
stats = {
    "n_chains": len(chains),
    "reads_at_first_donor_median": statistics.median(spanning.values()),
    "by_tolerance_bp": {str(w): {"median_full_reads": statistics.median(counts[w].values()),
                                 "chains_ge1": sum(v >= 1 for v in counts[w].values()),
                                 "chains_ge3": sum(v >= 3 for v in counts[w].values())} for w in W},
    "splice_site_offset_fraction_within_10bp": {str(k): round(v / tot, 4) for k, v in sorted(offsets.items())},
    "exact_full_ge3_by_consensus": {
        "single_caller": [sum(1 for t, v in counts[0].items() if v >= 3 and nc.get(t, 1) < 2), sum(1 for t in chains if nc.get(t, 1) < 2)],
        "ge2_caller": [sum(1 for t, v in counts[0].items() if v >= 3 and nc.get(t, 1) >= 2), sum(1 for t in chains if nc.get(t, 1) >= 2)]},
}
json.dump(stats, open(out_path, "w"), indent=2)
print(json.dumps(stats, indent=1))
