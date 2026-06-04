#!/usr/bin/env bash
# Cohort yield of the reference-bias rescue across FOUR GIAB individuals
# (HG001/NA12878, HG002/NA24385, HG003/NA24149, HG004/NA24143). Runs the same
# whole-genome annotation+variant scan as ../hg002/run_wholegenome.sh per individual and
# totals the rescues. Annotation+variant only (no RNA-seq / caller runs). Public GIAB data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
SCAN="$HERE/../hg002/scan_variant_axis.py"

# --- placeholders (fill in) -------------------------------------------------
GENOME=/path/to/GRCh38.primary_assembly.genome.fa
GBC=/path/to/gencode_by_chr            # GENCODE vNN pre-split per chromosome (chrN.gtf)
BCFTOOLS=bcftools; SAMTOOLS=samtools
# GIAB v4.2.1 GRCh38 benchmark VCFs (download from ftp-trace.ncbi.nlm.nih.gov/.../giab/release):
declare -A VCF=(
  [HG001]=/path/to/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG002]=/path/to/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG003]=/path/to/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  [HG004]=/path/to/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz )
# ----------------------------------------------------------------------------
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT"
printf 'sample\tcreated\trescued\tfalse_rescues\theld_circular\n' > "$OUT/cohort_summary.tsv"
for s in HG001 HG002 HG003 HG004; do
  C=0 R=0 F=0 H=0
  for chrom in chr{1..22}; do
    d="$OUT/$s/$chrom"; mkdir -p "$d"
    "$SAMTOOLS" faidx "$GENOME" "$chrom" > "$d/ref.fa"; "$SAMTOOLS" faidx "$d/ref.fa"
    "$BCFTOOLS" view -r "$chrom" -v snps -m2 -M2 "${VCF[$s]}" \
      | "$BCFTOOLS" view -e 'strlen(REF)!=1 || strlen(ALT)!=1' -Oz -o "$d/snps.vcf.gz"
    "$BCFTOOLS" index -t "$d/snps.vcf.gz"
    "$BCFTOOLS" consensus -f "$d/ref.fa" -H 1 "$d/snps.vcf.gz" > "$d/hap.fa"; "$SAMTOOLS" faidx "$d/hap.fa"
    python3 "$SCAN" "$d/ref.fa" "$d/hap.fa" "$GBC/$chrom.gtf" "$chrom" "$d" >/dev/null
    "$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" --ref-gtf "$d/vt.catalog.gtf" \
      --reference "$d/ref.fa" --reference-haplotype "$d/hap.fa" --haplotype-provenance wgs --out-prefix "$d/vt" >/dev/null
    "$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" --ref-gtf "$d/vt.catalog.gtf" \
      --reference "$d/ref.fa" --reference-haplotype "$d/hap.fa" --haplotype-provenance unknown --out-prefix "$d/vtc" >/dev/null
    pk(){ join -t$'\t' <(sort "$d/vt.truth") <(awk -F'\t' 'NR>1{print $1"\t"$7}' "$1"|sort); }
    C=$((C+$(awk -F'\t' '$2=="CREATED"' "$d/vt.truth"|wc -l)))
    R=$((R+$(pk "$d/vt.adjudicated.tsv"|awk -F'\t' '$2=="CREATED"&&$3=="PAN_REF_RESCUED_FALSE_NOVEL"'|wc -l)))
    F=$((F+$(pk "$d/vt.adjudicated.tsv"|awk -F'\t' '$2!="CREATED"&&$3=="PAN_REF_RESCUED_FALSE_NOVEL"'|wc -l)))
    H=$((H+$(pk "$d/vtc.adjudicated.tsv"|awk -F'\t' '$2=="CREATED"&&$3=="AMBIGUOUS"'|wc -l)))
  done
  printf '%s\t%s\t%s\t%s\t%s\n' "$s" "$C" "$R" "$F" "$H" | tee -a "$OUT/cohort_summary.tsv"
done
