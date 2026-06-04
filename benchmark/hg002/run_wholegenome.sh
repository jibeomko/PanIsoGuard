#!/usr/bin/env bash
# Whole-genome YIELD of the variant reference-bias rescue on REAL GIAB HG002 v4.2.1
# variants (extends run.sh from chr22 N=1 to chr1-22). Per chromosome: build an
# SNV-consensus haplotype, scan every GENCODE intron for a reference-bias splice site
# (CREATED = non-canonical on the linear reference, canonical on the HG002 haplotype),
# and adjudicate the variant axis. Annotation+variant scan only -- no RNA-seq / caller
# runs. Answers "how often does the rescue actually fire, and does it stay specific at
# scale?" Result envelope: ../results/hg002_wholegenome/metrics.json. No private data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# --- placeholders (fill in) -------------------------------------------------
GENOME=/path/to/GRCh38.primary_assembly.genome.fa     # + .fai
GENCODE=/path/to/gencode.vNN.annotation.gtf
VCF=/path/to/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz # GIAB (+ .tbi)
BCFTOOLS=bcftools                                      # >= 1.18
SAMTOOLS=samtools
# ----------------------------------------------------------------------------
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT/work" "$OUT/gencode_by_chr"; cd "$OUT"

# split GENCODE by main chromosome once
[ "$(ls gencode_by_chr | wc -l)" -ge 22 ] || \
  awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print >> "gencode_by_chr/"$1".gtf"}' "$GENCODE"

printf 'chrom\tintrons\tcreated\tdisrupted\tcontrol\trescued\tfalse_rescues\theld_circular\n' > summary.tsv
for chrom in chr{1..22}; do
  d="work/$chrom"; mkdir -p "$d"
  "$SAMTOOLS" faidx "$GENOME" "$chrom" > "$d/ref.fa"; "$SAMTOOLS" faidx "$d/ref.fa"
  "$BCFTOOLS" view -r "$chrom" -v snps -m2 -M2 "$VCF" \
    | "$BCFTOOLS" view -e 'strlen(REF)!=1 || strlen(ALT)!=1' -Oz -o "$d/snps.vcf.gz"
  "$BCFTOOLS" index -t "$d/snps.vcf.gz"
  "$BCFTOOLS" consensus -f "$d/ref.fa" -H 1 "$d/snps.vcf.gz" > "$d/hap.fa"; "$SAMTOOLS" faidx "$d/hap.fa"
  scan=$(python3 "$HERE/scan_variant_axis.py" "$d/ref.fa" "$d/hap.fa" "gencode_by_chr/$chrom.gtf" "$chrom" "$d")
  # variant axis: wgs provenance = independent (rescue should fire)
  "$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" --ref-gtf "$d/vt.catalog.gtf" \
    --reference "$d/ref.fa" --reference-haplotype "$d/hap.fa" --haplotype-provenance wgs --out-prefix "$d/vt" >/dev/null
  # circular provenance = firewall must HOLD (not promote)
  "$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" --ref-gtf "$d/vt.catalog.gtf" \
    --reference "$d/ref.fa" --reference-haplotype "$d/hap.fa" --haplotype-provenance unknown --out-prefix "$d/vt_circ" >/dev/null
  pick() { join -t$'\t' <(sort "$d/vt.truth") <(awk -F'\t' 'NR>1{print $1"\t"$7}' "$1"|sort); }
  introns=$(echo "$scan"|grep -oE 'introns scanned: [0-9]+'|grep -oE '[0-9]+')
  created=$(echo "$scan"|grep -oE 'CREATED\(reference-bias\)=[0-9]+'|grep -oE '[0-9]+')
  disr=$(echo "$scan"|grep -oE 'DISRUPTED=[0-9]+'|grep -oE '[0-9]+')
  ctrl=$(echo "$scan"|grep -oE 'CONTROL=[0-9]+'|grep -oE '[0-9]+')
  rescued=$(pick "$d/vt.adjudicated.tsv"|awk -F'\t' '$2=="CREATED"&&$3=="PAN_REF_RESCUED_FALSE_NOVEL"'|wc -l)
  false_r=$(pick "$d/vt.adjudicated.tsv"|awk -F'\t' '$2!="CREATED"&&$3=="PAN_REF_RESCUED_FALSE_NOVEL"'|wc -l)
  held=$(pick "$d/vt_circ.adjudicated.tsv"|awk -F'\t' '$2=="CREATED"&&$3=="AMBIGUOUS"'|wc -l)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$chrom" "$introns" "$created" "$disr" "$ctrl" "$rescued" "$false_r" "$held" | tee -a summary.tsv
done
awk -F'\t' 'NR>1{i+=$2;c+=$3;r+=$6;f+=$7;h+=$8} END{print "TOTAL introns="i" CREATED="c" rescued="r" false="f" held="h}' summary.tsv
