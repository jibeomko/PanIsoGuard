#!/usr/bin/env python3
"""Build the adjudicate inputs for the GM12878 caller-novel chains.

From the wholegenome_multicaller combine matrix + the three caller GTFs, emit:
  - novel_chains.gtf : one transcript per novel chain (exons from a representative native id)
  - novel.cls        : a minimal SQANTI-style classification (all novel_not_in_catalog)
  - pig_ncallers.json: PIG id -> n_callers (for the consensus stratification)

Usage: build_inputs.py <WORK: gm12878_mc dir> <OUT dir>
"""
import json, sys
from collections import defaultdict

WORK, OUT = sys.argv[1], sys.argv[2]
GTFS = {"bambu": f"{WORK}/bambu.novel.gtf",
        "isoquant": f"{WORK}/isoquant_out/gm/gm.transcript_models.gtf",
        "espresso": f"{WORK}/espresso_out/samples_N2_R0_updated.gtf"}

def tidof(a):
    return a.split('transcript_id "')[1].split('"')[0] if 'transcript_id "' in a else None

EX = defaultdict(lambda: defaultdict(list)); META = {}
for caller, path in GTFS.items():
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = tidof(f[8])
        if t is None:
            continue
        EX[caller][t].append((int(f[3]), int(f[4]))); META[(caller, t)] = (f[0], f[6])

HDR = ("isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
       "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
       "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n")
g = open(f"{OUT}/novel_chains.gtf", "w"); c = open(f"{OUT}/novel.cls", "w"); c.write(HDR)
ncallers = {}; n = 0
for ln in open(f"{WORK}/matrix3.tsv"):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "pig_id" or len(f) < 9 or f[7] != "novel":
        continue
    pig = f[0]; nc = int(f[4]); rep = None
    for pair in f[8].split(";"):
        if "=" not in pair:
            continue
        caller, t = pair.split("=", 1)
        if t in EX.get(caller, {}):
            rep = (caller, t); break
    if rep is None:
        continue
    caller, t = rep; chrom, strand = META[rep]
    for s, e in sorted(EX[caller][t]):
        g.write(f'{chrom}\t{caller}\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{pig}"; transcript_id "{pig}";\n')
    c.write(f"{pig}\t{chrom}\t{strand}\tnovel_not_in_catalog\t{pig}\tnovel\tmulti-exon\tFALSE\tcanonical\t"
            f"10.0\t0\tNA\tNA\tIsoform\n")
    ncallers[pig] = nc; n += 1
g.close(); c.close()
json.dump(ncallers, open(f"{OUT}/pig_ncallers.json", "w"))
print(f"built {n} novel-chain transcripts + classification")
