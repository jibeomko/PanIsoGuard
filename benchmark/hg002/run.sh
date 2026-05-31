#!/bin/bash
# Validate the variant (reference-bias) axis on REAL, independent (WGS-derived)
# variants: GIAB HG002 -> SNV-only personalized haplotype -> find splice motifs that
# the variant creates/disrupts -> PanIsoGuard adjudicate (provenance=wgs, non-circular).
#
# This is the acceptance-critical, non-circular validation of the headline
# PAN_REF_RESCUED_FALSE_NOVEL mechanism (reference bias != new splicing).
set -euo pipefail

# ---- config (edit) ----------------------------------------------------------
GENOME=/path/to/GRCh38.primary_assembly.genome.fa
GTF=/path/to/gencode.vNN.annotation.gtf
CHROM=chr22
BCFTOOLS=$HOME/miniconda3/envs/pigval/bin/bcftools     # bcftools >=1.18
SAMTOOLS=samtools
PIG=$(cd "$(dirname "$0")/../../build" && pwd)/panisoguard
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-./hg002_work}; mkdir -p "$OUT"; cd "$OUT"
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

# ---- 1) GIAB HG002 benchmark variants ---------------------------------------
[ -f hg002.vcf.gz ] || curl -fsSL -o hg002.vcf.gz "$GIAB"
[ -f hg002.vcf.gz.tbi ] || curl -fsSL -o hg002.vcf.gz.tbi "$GIAB.tbi" || "$BCFTOOLS" index -t hg002.vcf.gz

# ---- 2) reference contig + strict-SNV personalized haplotype ----------------
"$SAMTOOLS" faidx "$GENOME" "$CHROM" > ref.fa; "$SAMTOOLS" faidx ref.fa
# strict 1bp biallelic SNVs only so the haplotype stays coordinate-aligned to ref
"$BCFTOOLS" view -r "$CHROM" -v snps -m2 -M2 hg002.vcf.gz \
  | "$BCFTOOLS" view -e 'strlen(REF)!=1 || strlen(ALT)!=1' -Oz -o snps.vcf.gz
"$BCFTOOLS" index -t snps.vcf.gz
"$BCFTOOLS" consensus -f ref.fa -H 1 snps.vcf.gz > hap.fa
"$SAMTOOLS" faidx hap.fa
[ "$(cut -f2 ref.fa.fai)" = "$(cut -f2 hap.fa.fai)" ] || { echo "ERROR: ref/hap length mismatch"; exit 1; }

# ---- 3) chromosome annotation + scan/build labelled variant-axis test -------
awk -F'\t' -v c="$CHROM" '$1==c' "$GTF" > chrom.gtf
python "$HERE/scan_variant_axis.py" ref.fa hap.fa chrom.gtf "$CHROM" .

# ---- 4) PanIsoGuard adjudicate (real HG002 haplotype, WGS provenance) --------
"$PIG" adjudicate --classification vt.cls --isoforms-gtf vt.gtf --ref-gtf vt.catalog.gtf \
  --reference ref.fa --reference-haplotype hap.fa --haplotype-provenance wgs --out-prefix vt

# ---- 5) report verdict by label ---------------------------------------------
echo "=== verdict by label (expect CREATED -> PAN_REF_RESCUED_FALSE_NOVEL) ==="
join -t$'\t' <(sort vt.truth) \
     <(awk -F'\t' 'NR>1{print $1"\t"$7}' vt.adjudicated.tsv | sort) \
  | awk -F'\t' '{print "  "$2"\t"$1"\t"$3}' | sort
