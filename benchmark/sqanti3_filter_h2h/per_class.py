#!/usr/bin/env python3
"""Per-confidence-class precision on the SQANTI-SIM truth, for every short-read setting of the
head-to-head run. Same labels as score_h2h.py (exact intron-chain match; multi-exon only).

Usage: per_class.py <work> <h2h_dir>   (work has flair.isoforms.gtf, truth.truth.tsv,
       sqanti_out/flair_classification.txt; h2h_dir has <setting>/pig*.adjudicated.tsv)
"""
import collections, os, sys

work, h2h = sys.argv[1:3]
LEGACY = {"HIGH_CONF_KNOWN": "KNOWN", "HIGH_CONF_NOVEL": "NOVEL_CONFIRMED", "MEDIUM_CONF_NOVEL": "NOVEL_PROBABLE",
          "LOW_CONF_PARTIAL": "UNCONFIRMED", "AMBIGUOUS": "UNRESOLVED", "PAN_REF_RESCUED_FALSE_NOVEL": "REFERENCE_BIAS"}
ex, strand, chrom = collections.defaultdict(list), {}, {}
for line in open(f"{work}/flair.isoforms.gtf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    t = f[8].split('transcript_id "')[1].split('"')[0]
    ex[t].append((int(f[3]), int(f[4]))); strand[t] = f[6]; chrom[t] = f[0]
truth = dict(l.rstrip("\n").split("\t") for l in open(f"{work}/truth.truth.tsv"))
label = {}
for t, e in ex.items():
    e = sorted(e)
    if len(e) < 2: continue
    key = f"{chrom[t]}|{strand[t]}|" + ",".join(f"{e[i][1] + 1}-{e[i + 1][0] - 1}" for i in range(len(e) - 1))
    label[t] = truth.get(key, "false_novel")
label = {t: l for t, l in label.items() if l in ("genuine_novel", "false_novel")}
pos_rate = sum(l == "genuine_novel" for l in label.values()) / len(label)
print(f"labelled multi-exon isoforms: {len(label)}  genuine={sum(l == 'genuine_novel' for l in label.values())}  positive rate={pos_rate:.3f}")
order = ["KNOWN", "NOVEL_CONFIRMED", "NOVEL_PROBABLE", "UNRESOLVED", "UNCONFIRMED", "ARTIFACT"]
for setting in ["none", "oracle", "sr_250000", "sr_1000000", "sr_4000000"]:
    for run in ["pig", "pigbam"]:
        p = f"{h2h}/{setting}/{run}.adjudicated.tsv"
        if not os.path.exists(p): continue
        it = (l.rstrip("\n").split("\t") for l in open(p)); hdr = next(it)
        cls = {r[0]: LEGACY.get(r[hdr.index("confidence_class")], r[hdr.index("confidence_class")]) for r in it}
        c = collections.Counter((cls[t], label[t]) for t in label if t in cls)
        print(f"\n=== {setting} / {run}")
        print(f"{'class':<20}{'genuine':>8}{'false':>7}{'precision':>11}")
        for k in order:
            g, f_ = c[(k, "genuine_novel")], c[(k, "false_novel")]
            if g + f_: print(f"{k:<20}{g:>8}{f_:>7}{g / (g + f_):>11.3f}")
