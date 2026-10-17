#!/usr/bin/env python3
"""Join LRGASP Supplementary Data 15 (WTC11 validation targets) with the WTC11 consolidated
models BED to get each target's intron chain (GRCh38, 1-based inclusive introns).

Usage: targets.py wtc11_validation_targets.tsv WTC11_consolidated.bed.gz outdir
Writes outdir/targets.tsv (target, transcript, type, category, detail, chain_key, n_introns, ...)
and outdir/targets.gtf (the target models, for inspection).
"""
import csv, gzip, sys

tsv, bed, out = sys.argv[1:4]
targets = list(csv.DictReader(open(tsv), delimiter="\t"))
want = {t["transcript"] for t in targets}
models = {}
with gzip.open(bed, "rt") as f:
    for l in f:
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        if c[3] not in want: continue
        chrom, start, strand = c[0], int(c[1]), c[5]
        sizes = [int(x) for x in c[10].rstrip(",").split(",")]; offs = [int(x) for x in c[11].rstrip(",").split(",")]
        exons = [(start + o, start + o + s) for o, s in zip(offs, sizes)]  # 0-based half-open
        introns = [(exons[i][1] + 1, exons[i + 1][0]) for i in range(len(exons) - 1)]  # 1-based inclusive
        models[c[3]] = (chrom, strand, exons, introns)
with open(f"{out}/targets.tsv", "w") as o, open(f"{out}/targets.gtf", "w") as g:
    w = csv.writer(o, delimiter="\t")
    w.writerow(["target", "transcript", "type", "category", "detail", "position", "chrom", "strand", "n_introns", "chain_key",
                "pb_chain", "ont_chain", "mediaCpm", "pbPlatformCount", "ontPlatformCount", "notes"])
    n_missing = 0
    for t in targets:
        m = models.get(t["transcript"])
        if m is None:
            n_missing += 1
            w.writerow([t["target"], t["transcript"], t["type"], t["category"], t["detail"], t["position"], "", "", "", "NA",
                        t["pb_chain"], t["ont_chain"], t["mediaCpm"], t["pbPlatformCount"], t["ontPlatformCount"], t["notes"]])
            continue
        chrom, strand, exons, introns = m
        key = f"{chrom}|{strand}|" + ",".join(f"{s}-{e}" for s, e in introns) if introns else "monoexonic"
        w.writerow([t["target"], t["transcript"], t["type"], t["category"], t["detail"], t["position"], chrom, strand, len(introns), key,
                    t["pb_chain"], t["ont_chain"], t["mediaCpm"], t["pbPlatformCount"], t["ontPlatformCount"], t["notes"]])
        for s, e in exons:
            g.write(f'{chrom}\tlrgasp\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tgene_id "{t["target"]}"; transcript_id "{t["transcript"]}";\n')
print(f"targets {len(targets)}, models found {len(targets) - n_missing}, missing {n_missing}", file=sys.stderr)
