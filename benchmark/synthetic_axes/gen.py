#!/usr/bin/env python3
"""Synthetic-genome controlled truth for the BAM (mapping) and variant axes.

Builds a synthetic contig with four labelled isoform categories so each evidence
axis can be validated end-to-end through `panisoguard adjudicate` against known
truth (no real genome / simulator needed):

  A genuine            canonical ref motif + short-read support  -> HIGH_CONF_NOVEL
  B noncanon-artifact  non-canonical motif + no SR support       -> ARTIFACT (noncanonical)
  C mapping-artifact   canonical + no SR + low-MAPQ spanning BAM -> ARTIFACT (mapping)   [needs --bam]
  D reference-bias     non-canonical on ref, canonical on hap    -> PAN_REF_RESCUED      [needs --reference-haplotype]

Each isoform is a 2-exon transcript with one test intron in its own 1 kb window.

Usage: python gen.py [outdir=.] [n_per_category=40]
"""
import sys

OUT = (sys.argv[1] if len(sys.argv) > 1 else ".").rstrip("/")
NPC = int(sys.argv[2]) if len(sys.argv) > 2 else 40
CATS = ["A", "B", "C", "D"]
N = NPC * len(CATS)
LEN = 1000 * (N + 2)

# window/intron geometry for isoform i (0-based genome coords)
def geom(i):
    b = 1000 * (i + 1)
    return b, (b + 100, b + 300)  # exon1=[b,b+100) ; intron=[b+100,b+300) ; exon2=[b+300,b+400)

ref = bytearray(b"A" * LEN)
hap = bytearray(b"A" * LEN)

def put(buf, pos, s):
    buf[pos:pos + len(s)] = s.encode()

cat_of = {}
for i in range(N):
    cat = CATS[i % len(CATS)]
    cat_of[i] = cat
    _, (ds, de) = geom(i)
    donor_ref = "GT" if cat in ("A", "C") else "GG"   # canonical vs non-canonical on reference
    put(ref, ds, donor_ref); put(ref, de - 2, "AG")
    donor_hap = "GT" if cat == "D" else donor_ref       # variant flips D to canonical on haplotype
    put(hap, ds, donor_hap); put(hap, de - 2, "AG")

def write_fa(path, buf):
    with open(path, "w") as f:
        f.write(">chrS\n")
        s = buf.decode()
        for k in range(0, len(s), 70):
            f.write(s[k:k + 70] + "\n")

write_fa(f"{OUT}/ref.fa", ref)
write_fa(f"{OUT}/hap.fa", hap)

# decoy "known" transcript so the catalog is non-empty (no test intron is in it)
with open(f"{OUT}/catalog.gtf", "w") as cat:
    s = LEN - 600
    for (a, b) in [(s, s + 100), (s + 300, s + 400)]:
        cat.write(f'chrS\tsyn\texon\t{a+1}\t{b}\t.\t+\t.\tgene_id "DECOY"; transcript_id "DECOY1";\n')

HEADER = ("isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
          "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
          "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n")
gtf = open(f"{OUT}/caller.gtf", "w")
cls = open(f"{OUT}/classification.tsv", "w"); cls.write(HEADER)
sj = open(f"{OUT}/real.SJ.tab", "w")
truth = open(f"{OUT}/truth.tsv", "w"); truth.write("isoform\tcategory\ttruth\n")
sam = open(f"{OUT}/syn.sam", "w")
sam.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrS\tLN:%d\n" % LEN)

truth_label = {"A": "genuine", "B": "false", "C": "false", "D": "reference_bias"}

def reads(iso, b, mapq, flag, n):
    for r in range(n):
        # 100M200N100M starting at 0-based b -> N op spans the intron [b+100,b+300)
        sam.write(f"{iso}_r{r}\t{flag}\tchrS\t{b+1}\t{mapq}\t100M200N100M\t*\t0\t0\t*\t*\n")

for i in range(N):
    cat = cat_of[i]
    b, (ds, de) = geom(i)
    iso = f"{cat}_{i:04d}"
    # caller GTF: exon1 [b,b+100), exon2 [b+300,b+400)  (1-based GTF)
    gtf.write(f'chrS\tsyn\texon\t{b+1}\t{b+100}\t.\t+\t.\tgene_id "{iso}"; transcript_id "{iso}";\n')
    gtf.write(f'chrS\tsyn\texon\t{b+301}\t{b+400}\t.\t+\t.\tgene_id "{iso}"; transcript_id "{iso}";\n')
    canon = "canonical" if cat in ("A", "C") else "non_canonical"
    cls.write(f"{iso}\tchrS\t+\tnovel_not_in_catalog\t{iso}\tnovel\tmulti-exon\tFALSE\t{canon}\t"
              f"10.0\t0\tNA\tNA\tIsoform\n")
    truth.write(f"{iso}\t{cat}\t{truth_label[cat]}\n")
    # short-read support: category A only (1-based intron [ds+1, de])
    if cat == "A":
        sj.write(f"chrS\t{ds+1}\t{de}\t1\t1\t1\t20\t0\t30\n")
    # BAM coverage
    if cat == "A":
        reads(iso, b, 60, 0, 3)         # clean
    elif cat == "C":
        reads(iso, b, 5, 0, 3)          # low-MAPQ -> mapping artifact
    else:                                # B, D: a couple of clean reads (axis evaluable, not flagged)
        reads(iso, b, 60, 0, 2)

for fh in (gtf, cls, sj, truth, sam):
    fh.close()
print(f"contig chrS len={LEN}  isoforms={N} ({NPC}/category)")
