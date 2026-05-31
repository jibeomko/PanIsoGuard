# HG002 variant-axis validation (real, non-circular)

The acceptance-critical validation of the headline `PAN_REF_RESCUED_FALSE_NOVEL`
mechanism on **real, independent (WGS-derived) variants** — GIAB HG002 — so the
rescue evidence is genuinely independent of the RNA being adjudicated (the
circularity firewall's `wgs` provenance).

## Method

1. Download the GIAB HG002 v4.2.1 small-variant benchmark VCF (public, NIST FTP).
2. Build a personalized haplotype FASTA with `bcftools consensus` using **strict
   1 bp biallelic SNVs only**, so it stays coordinate-aligned to GRCh38 (indels are
   excluded — the documented coordinate-shift caveat).
3. Scan the reference annotation's introns for splice donor/acceptor dinucleotides
   that differ between GRCh38 and the HG002 haplotype:
   - **CREATED** — non-canonical on GRCh38, canonical on HG002 = reference bias
   - **DISRUPTED** — canonical on GRCh38, non-canonical on HG002
   - **CONTROL** — canonical on both
4. Build adjudicate inputs from these loci and run PanIsoGuard with
   `--reference-haplotype hap.fa --haplotype-provenance wgs`.

## Run

```bash
# edit GENOME / GTF / tool paths at the top of run.sh, then:
bash run.sh hg002_work
```
Needs `bcftools` (e.g. conda env with bcftools ≥1.18), `samtools`, and `panisoguard`.

## Result (GENCODE v49 chr22, GIAB HG002 v4.2.1)

Of 13,496 unique chr22 introns, exactly **1 CREATED** (reference-bias) and **2
DISRUPTED** splice motifs are produced by real HG002 SNVs (annotated introns are
overwhelmingly canonical on GRCh38, so motif-altering SNVs are rare — as expected).

```
label           verdict
CREATED         PAN_REF_RESCUED_FALSE_NOVEL     <- chr22:23688305-23688721  GT-GG (ref) -> GT-AG (HG002)
DISRUPTED ×2    AMBIGUOUS                       (not falsely rescued)
CONTROL  ×5     AMBIGUOUS                       (not falsely rescued)
```

PanIsoGuard's variant axis, fed the real HG002 WGS-derived haplotype (non-circular),
correctly rescued the single genuine reference-bias junction — a junction that looks
non-canonical/suspect on GRCh38 but is a clean GT-AG on HG002's own genome — and
produced **zero false rescues** on the disrupted/control loci. This validates the
headline mechanism on real human-genome variation.

(For a controlled, high-N version of the same rescue logic see
[../synthetic_axes](../synthetic_axes); for the short-read axis see
[../controlled_truth](../controlled_truth).)
