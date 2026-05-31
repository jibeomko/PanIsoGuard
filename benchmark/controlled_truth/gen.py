#!/usr/bin/env python3
"""Controlled-truth E1 generator (incomplete-reference protocol).

Builds a labelled benchmark for PanIsoGuard's adjudication from a reference GTF
(e.g. one chromosome of GENCODE):

  * hide ~20% of multi-exon transcripts -> they become GENUINE novel truth (real,
    but absent from the catalog handed to PanIsoGuard);
  * fabricate FALSE novel isoforms by shifting one internal exon boundary of a
    kept transcript (creates a junction that is neither annotated nor read-supported);
  * every REAL junction gets short-read support in the emitted SJ.tab; fabricated
    junctions do not.

A correct adjudicator should call genuine novels HIGH/MEDIUM and fabricated ones
LOW/ARTIFACT, and may honestly abstain (AMBIGUOUS) on genuine isoforms whose every
junction is already known (no novel junction to corroborate).

Usage:  python gen.py <reference.gtf> [outdir=.]
Tip:    awk -F'\\t' '$1=="chr22"' gencode.vNN.annotation.gtf > chr22.gtf
"""
import re, sys

GTF = sys.argv[1] if len(sys.argv) > 1 else "chr22.gtf"
OUT = (sys.argv[2] if len(sys.argv) > 2 else ".").rstrip("/")
tid_re = re.compile(r'transcript_id "([^"]+)"')
gid_re = re.compile(r'gene_id "([^"]+)"')

tx = {}
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = tid_re.search(f[8])
    if not m:
        continue
    g = gid_re.search(f[8])
    d = tx.setdefault(m.group(1), {"chrom": f[0], "strand": f[6],
                                   "gene": g.group(1) if g else "NA", "exons": []})
    d["exons"].append((int(f[3]), int(f[4])))

multi = {t: d for t, d in tx.items() if len(d["exons"]) >= 2}
for d in multi.values():
    d["exons"].sort()

tids = sorted(multi)
hidden = [t for i, t in enumerate(tids) if i % 5 == 0]
kept = [t for i, t in enumerate(tids) if i % 5 != 0]

def introns_1based(exons):
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]

with open(f"{OUT}/reduced_catalog.gtf", "w") as cat:
    for tid in kept:
        d = multi[tid]
        for (s, e) in d["exons"]:
            cat.write(f'{d["chrom"]}\tctrl\texon\t{s}\t{e}\t.\t{d["strand"]}\t.\t'
                      f'gene_id "{d["gene"]}"; transcript_id "{tid}";\n')

strand_code = {"+": 1, "-": 2}
seen = set()
with open(f"{OUT}/real.SJ.tab", "w") as sj:
    for d in multi.values():
        for (a, b) in introns_1based(d["exons"]):
            if b < a or (d["chrom"], a, b) in seen:
                continue
            seen.add((d["chrom"], a, b))
            sj.write(f'{d["chrom"]}\t{a}\t{b}\t{strand_code.get(d["strand"], 0)}\t1\t1\t20\t0\t30\n')

HEADER = ("isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
          "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
          "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n")
gtf = open(f"{OUT}/caller.gtf", "w")
cls = open(f"{OUT}/classification.tsv", "w"); cls.write(HEADER)
truth = open(f"{OUT}/truth.tsv", "w"); truth.write("isoform\ttruth\n")

def emit(iso, d, exons):
    for (s, e) in exons:
        gtf.write(f'{d["chrom"]}\tctrl\texon\t{s}\t{e}\t.\t{d["strand"]}\t.\t'
                  f'gene_id "{d["gene"]}"; transcript_id "{iso}";\n')
    cls.write(f"{iso}\t{d['chrom']}\t{d['strand']}\tnovel_not_in_catalog\t{d['gene']}\tnovel\t"
              f"multi-exon\tFALSE\tcanonical\t10.0\t0\tNA\tNA\tIsoform\n")

npos = nneg = 0
for tid in hidden:
    emit(f"POS_{tid}", multi[tid], multi[tid]["exons"])
    truth.write(f"POS_{tid}\tgenuine\n"); npos += 1

for i, tid in enumerate(kept):
    if i % 7 != 0:
        continue
    d = multi[tid]
    ex = [list(p) for p in d["exons"]]
    j = next((k for k in range(1, len(ex) - 1) if ex[k][1] - ex[k][0] > 12), None)
    if j is None:
        continue
    ex[j][0] += 6  # fabricated junction
    emit(f"NEG_{tid}", d, [tuple(p) for p in ex])
    truth.write(f"NEG_{tid}\tfalse\n"); nneg += 1

gtf.close(); cls.close(); truth.close()
print(f"multi-exon={len(multi)} hidden={len(hidden)} kept={len(kept)} "
      f"positives={npos} negatives={nneg} real_junctions={len(seen)}")
