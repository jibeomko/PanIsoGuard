#!/usr/bin/env python3
"""Multi-contig incomplete-reference prep for the SIRV spike-in control.

Hides ~20% of multi-exon SIRV transcripts (genuine-novel truth) and emits PBSIM3
sim input, a reduced annotation, contig-aware truth chains, and a short-read SJ.tab.
Usage: python prep.py <genome.fa(+.fai)> <annotation.gtf> <outdir>
"""
import re, sys

GEN, ANN, W = sys.argv[1], sys.argv[2], sys.argv[3].rstrip("/")

seq, cur = {}, None
for ln in open(GEN):
    if ln.startswith(">"):
        cur = ln[1:].split()[0]; seq[cur] = []
    else:
        seq[cur].append(ln.strip())
seq = {k: "".join(v).upper() for k, v in seq.items()}
comp = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(comp)[::-1]

tid = re.compile(r'transcript_id "([^"]+)"')
tx = {}
for ln in open(ANN):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = tid.search(f[8])
    if m:
        tx.setdefault(m.group(1), {"c": f[0], "st": f[6], "ex": []})["ex"].append((int(f[3]), int(f[4])))

multi = {t: d for t, d in tx.items() if len(d["ex"]) >= 2}
for d in multi.values():
    d["ex"].sort()
tids = sorted(multi); hidden = set(tids[::5])
tseq = lambda d: (lambda s: rc(s) if d["st"] == "-" else s)("".join(seq[d["c"]][a - 1:b] for a, b in d["ex"]))
intr = lambda d: [(d["ex"][i][1], d["ex"][i + 1][0] - 1) for i in range(len(d["ex"]) - 1)]

with open(f"{W}/sim.transcript", "w") as o:
    for t in tids:
        s = tseq(multi[t])
        if len(s) >= 300:
            o.write(f"{t}\t60\t0\t{s}\n")
with open(f"{W}/reduced.gtf", "w") as o:
    for ln in open(ANN):
        m = tid.search(ln)
        if m and m.group(1) in hidden:
            continue
        o.write(ln)
with open(f"{W}/truth_hidden_chains.tsv", "w") as o:
    for t in hidden:
        d = multi[t]
        o.write(f"{t}\t{d['c']}\t" + ";".join(f"{a}-{b}" for a, b in intr(d)) + "\n")
sc = {"+": 1, "-": 2}; seen = set()
with open(f"{W}/real.SJ.tab", "w") as o:
    for t in tids:
        d = multi[t]
        for a, b in intr(d):
            if b < a or (d["c"], a, b) in seen:
                continue
            seen.add((d["c"], a, b))
            o.write(f"{d['c']}\t{a+1}\t{b}\t{sc[d['st']]}\t1\t1\t30\t0\t30\n")
print(f"multi-exon SIRV tx={len(multi)} hidden={len(hidden)} real_junctions={len(seen)}")
