#!/usr/bin/env python3
"""Find real splice-motif-altering variants and build a labelled variant-axis test.

Scans a reference annotation's introns for splice donor/acceptor dinucleotides
that DIFFER between the linear reference and a personalized (SNV-consensus)
haplotype FASTA, classifies each as CREATED (non-canonical on reference, canonical
on the haplotype = reference bias), DISRUPTED (canonical -> non-canonical), or
CONTROL (canonical on both), and emits PanIsoGuard adjudicate inputs so the variant
axis can be validated on REAL variants.

Usage: python scan_variant_axis.py <ref.fa> <hap.fa> <chrom.gtf> <chrom> <outdir>
"""
import re, sys

ref_fa, hap_fa, gtf, CHROM, OUT = sys.argv[1:6]
OUT = OUT.rstrip("/")

def load(p):
    s = []
    for ln in open(p):
        if not ln.startswith(">"):
            s.append(ln.strip())
    return "".join(s).upper()

ref, hap = load(ref_fa), load(hap_fa)
assert len(ref) == len(hap), "ref and hap must be length-matched (use SNV-only consensus)"
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}

def canonical(seq, a, b, strand):
    d, ac = seq[a:a + 2], seq[b - 2:b]
    if strand == "-":
        d, ac = rc(ac), rc(d)
    return (d, ac) in CANON

tid = re.compile(r'transcript_id "([^"]+)"')
tx = {}
for ln in open(gtf):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = tid.search(f[8])
    if m:
        tx.setdefault(m.group(1), {"st": f[6], "ex": []})["ex"].append((int(f[3]), int(f[4])))

introns = set()
for d in tx.values():
    e = sorted(d["ex"])
    for i in range(len(e) - 1):
        a, b = e[i][1], e[i + 1][0] - 1   # 0-based half-open intron [a, b)
        if b > a:
            introns.add((a, b, d["st"]))

created, disrupted, control = [], [], []
for (a, b, st) in introns:
    if b - a < 4 or a < 100 or b + 100 > len(ref):
        continue
    cr, ch = canonical(ref, a, b, st), canonical(hap, a, b, st)
    if not cr and ch:
        created.append((a, b, st))
    elif cr and not ch:
        disrupted.append((a, b, st))
    elif cr and ch and len(control) < 5:
        control.append((a, b, st))

HDR = ("isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
       "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
       "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n")
g = open(f"{OUT}/vt.gtf", "w"); c = open(f"{OUT}/vt.cls", "w"); c.write(HDR); t = open(f"{OUT}/vt.truth", "w")
for label, lst in [("CREATED", created), ("DISRUPTED", disrupted), ("CONTROL_canon", control)]:
    for k, (a, b, st) in enumerate(lst):
        iso = f"{label}_{k}_{a}"
        g.write(f'{CHROM}\tt\texon\t{a-99}\t{a}\t.\t{st}\t.\tgene_id "{iso}"; transcript_id "{iso}";\n')
        g.write(f'{CHROM}\tt\texon\t{b+1}\t{b+100}\t.\t{st}\t.\tgene_id "{iso}"; transcript_id "{iso}";\n')
        canon = "canonical" if canonical(ref, a, b, st) else "non_canonical"
        c.write(f"{iso}\t{CHROM}\t{st}\tnovel_not_in_catalog\t{iso}\tnovel\tmulti-exon\tFALSE\t{canon}\t"
                f"10.0\t0\tNA\tNA\tIsoform\n")
        t.write(f"{iso}\t{label}\n")
open(f"{OUT}/vt.catalog.gtf", "w").write(
    f'{CHROM}\tt\texon\t100\t200\t.\t+\t.\tgene_id "D"; transcript_id "D";\n'
    f'{CHROM}\tt\texon\t400\t500\t.\t+\t.\tgene_id "D"; transcript_id "D";\n')
g.close(); c.close(); t.close()
print(f"introns scanned: {len(introns)}  CREATED(reference-bias)={len(created)} "
      f"DISRUPTED={len(disrupted)} CONTROL={len(control)}")
