#!/usr/bin/env python3
"""Build a truth-derived short-read SJ.tab + an intron-chain truth map for the
SQANTI-SIM benchmark.

Truth (no caller involved):
  simulated transcripts        = source transcripts in read_to_isoform.tsv
  genuine-novel (= "deleted")  = simulated transcripts NOT in the reduced GTF
                                 (SQANTI-SIM removed them, so a caller must rediscover them)
  known                        = simulated transcripts still in the reduced GTF

The SJ.tab lists the introns of ALL simulated transcripts (they are expressed, so
short reads would corroborate them) -> a genuine-novel isoform's junctions are
SR-supported, while a FLAIR/alignment artifact's novel junction is not.

Usage: prep_truth_sj.py <full.gtf> <reduced.gtf> <read_to_isoform.tsv> <out_prefix>
Writes: <out_prefix>.SJ.tab  and  <out_prefix>.truth.tsv (chain_key<TAB>label)
"""
import sys
from collections import defaultdict

full_gtf, reduced_gtf, r2i, out = sys.argv[1:5]

def tid_of(attr):
    return attr.split('transcript_id "')[1].split('"')[0]

def exons_by_tid(path, keep=None):
    ex = defaultdict(list); strand = {}; chrom = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = tid_of(f[8])
        if keep is not None and t not in keep:
            continue
        ex[t].append((int(f[3]), int(f[4])))
        strand[t] = f[6]; chrom[t] = f[0]
    return ex, strand, chrom

def introns(exlist):
    e = sorted(exlist)
    return tuple((e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1))  # 1-based inclusive

# simulated source transcripts
sim = set(l.rstrip('\n').split('\t')[1] for l in open(r2i))
# transcript ids present in the reduced annotation (= known); absent = genuine-novel
reduced_ids = set()
for line in open(reduced_gtf):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) >= 9 and f[2] == 'transcript':
        reduced_ids.add(tid_of(f[8]))

ex, strand, chrom = exons_by_tid(full_gtf, keep=sim)

sj = {}                 # (chrom, s, e, strandcode) -> 1
truth = {}              # chain_key -> label
n_gen = n_known = n_mono = 0
for t in sim:
    if t not in ex:
        continue
    intr = introns(ex[t])
    label = 'known' if t in reduced_ids else 'genuine_novel'
    if not intr:        # monoexonic: no intron chain
        n_mono += 1
        continue
    sc = 1 if strand[t] == '+' else 2
    key = f"{chrom[t]}|{strand[t]}|" + ",".join(f"{s}-{e}" for s, e in sorted(intr))
    # a chain shared by both a kept and a deleted isoform counts as known (it is annotated)
    if key not in truth or label == 'known':
        truth[key] = label
    for (s, e) in intr:
        sj[(chrom[t], s, e, sc)] = 1
    if label == 'genuine_novel':
        n_gen += 1
    else:
        n_known += 1

with open(out + '.SJ.tab', 'w') as o:
    for (c, s, e, sc) in sj:
        o.write(f"{c}\t{s}\t{e}\t{sc}\t1\t1\t50\t0\t30\n")  # canonical, annotated, n_uniq=50
with open(out + '.truth.tsv', 'w') as o:
    for k, v in truth.items():
        o.write(f"{k}\t{v}\n")

print(f"simulated={len(sim)} genuine_novel(multi-exon)={n_gen} known={n_known} "
      f"monoexonic_skipped={n_mono}")
print(f"SJ.tab junctions={len(sj)}  truth chains={len(truth)}")
