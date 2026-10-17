#!/usr/bin/env python3
"""Join LRGASP Supplementary Data 15 with the LRGASP transcript-model BEDs to get each
validation target's intron chain (GRCh38, 1-based inclusive introns).

  sheet wtc11-support   -> targets.tsv          (NIC/NNC/FSM models from WTC11_consolidated.bed.gz)
  sheet gencode-support -> targets_gencode.tsv  (GENCODE-reviewed known/novel/rejected models from
                           LRGASP_manual_annotation.human.all.bed.gz + human_GENCODE_tmerge_transcripts.bed.gz)

Usage: targets.py <supp dir> <models dir> <outdir>
"""
import csv, gzip, os, sys

supp, models_dir, out = sys.argv[1:4]


def bed_models(path, want):
    found = {}
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("#"): continue
            c = l.rstrip("\n").split("\t")
            name = c[3].split(".")[0] if c[3].startswith("ENST") else c[3]
            if name not in want: continue
            chrom, start, strand = c[0], int(c[1]), c[5]
            sizes = [int(x) for x in c[10].rstrip(",").split(",")]; offs = [int(x) for x in c[11].rstrip(",").split(",")]
            exons = [(start + o, start + o + s) for o, s in zip(offs, sizes)]
            introns = [(exons[i][1] + 1, exons[i + 1][0]) for i in range(len(exons) - 1)]
            found[name] = (chrom, strand, exons, introns)
    return found


def write(targets, models, type_of, tsv_out, gtf_out):
    n_missing = 0
    with open(tsv_out, "w") as o, open(gtf_out, "w") as g:
        w = csv.writer(o, delimiter="\t")
        w.writerow(["target", "transcript", "type", "category", "detail", "position", "chrom", "strand", "n_introns", "chain_key",
                    "pb_chain", "ont_chain", "mediaCpm", "pbPlatformCount", "ontPlatformCount", "notes"])
        for t in targets:
            m = models.get(t["transcript"])
            row = [t["target"], t["transcript"], type_of(t), t["category"], t.get("detail", ""), t["position"]]
            if m is None:
                n_missing += 1
                w.writerow(row + ["", "", "", "NA", t.get("pb_chain", ""), t.get("ont_chain", ""), t.get("mediaCpm", ""), t.get("pbPlatformCount", ""), t.get("ontPlatformCount", ""), t.get("notes", "")])
                continue
            chrom, strand, exons, introns = m
            key = f"{chrom}|{strand}|" + ",".join(f"{s}-{e}" for s, e in introns) if introns else "monoexonic"
            w.writerow(row + [chrom, strand, len(introns), key, t.get("pb_chain", ""), t.get("ont_chain", ""), t.get("mediaCpm", ""), t.get("pbPlatformCount", ""), t.get("ontPlatformCount", ""), t.get("notes", "")])
            for s, e in exons:
                g.write(f'{chrom}\tlrgasp\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tgene_id "{t["target"]}"; transcript_id "{t["transcript"]}";\n')
    return n_missing


os.makedirs(out, exist_ok=True)
wtc = list(csv.DictReader(open(f"{supp}/wtc11_validation_targets.tsv"), delimiter="\t"))
m = bed_models(f"{models_dir}/WTC11_consolidated.bed.gz", {t["transcript"] for t in wtc})
miss = write(wtc, m, lambda t: t["type"], f"{out}/targets.tsv", f"{out}/targets.gtf")
print(f"wtc11-support: {len(wtc)} targets, {len(wtc) - miss} chains found", file=sys.stderr)

gpath = f"{supp}/gencode_validation_targets.tsv"
if os.path.exists(gpath):
    gen = list(csv.DictReader(open(gpath), delimiter="\t"))
    want = {t["transcript"] for t in gen}
    m = bed_models(f"{models_dir}/LRGASP_manual_annotation.human.all.bed.gz", want)
    m.update(bed_models(f"{models_dir}/human_GENCODE_tmerge_transcripts.bed.gz", want))
    miss = write(gen, m, lambda t: t["gencode"], f"{out}/targets_gencode.tsv", f"{out}/targets_gencode.gtf")
    print(f"gencode-support: {len(gen)} targets, {len(gen) - miss} chains found", file=sys.stderr)
