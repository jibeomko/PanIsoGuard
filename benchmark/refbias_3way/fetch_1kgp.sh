#!/usr/bin/env bash
# The 1000 Genomes high-coverage phased call set for the given samples: SNVs, indels and SVs
# called from 30x Illumina reads and phased with the pedigree (Byrska-Bishop et al., Cell 2022;
# 3,202 samples on GRCh38, chr1-22 and chrX). The panel is 30 GB and the server serves about
# 0.1 MB/s per connection, so only the regions of a BED file are read, through the remote tabix
# index, and only the samples' own columns are kept: one VCF per sample with the records where
# that sample carries a non-reference allele. INFO is dropped except END / SVLEN / SVTYPE of
# symbolic alleles.
#
# Usage: fetch_1kgp.sh OUTDIR REGIONS.bed SAMPLE [SAMPLE ...]     (env: JOBS, default 24 connections)
#   REGIONS.bed: 0-based, sorted, not overlapping (run.sh C_1kgp writes one: the span of every
#   isoform with a novel junction, +- 1 kb). Several samples can share one fetch: merge their BEDs.
#   writes OUTDIR/<SAMPLE>.kgp.vcf.gz(.tbi) and, with the phase removed from every genotype,
#   OUTDIR/<SAMPLE>.kgp.unphased.vcf.gz(.tbi)
set -euo pipefail
OUT=$(mkdir -p "$1" && cd "$1" && pwd); BED=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2; SAMPLES=("$@"); JOBS=${JOBS:-24}
U=${KGP_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV}
T="$OUT/tmp"; mkdir -p "$T"; cd "$T"
name() { [ "$1" = chrX ] && echo "1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" || echo "1kGP_high_coverage_Illumina.$1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"; }
CHROMS=$(printf '%s\n' chr{1..22} chrX)

# the header of one file (the sample order is the same in all of them): the first 3 MB hold it
[ -s header.txt ] || { curl -s -r 0-3000000 -o head.vcf.gz "$U/$(name chr22)"
  { bgzip -dc head.vcf.gz 2>/dev/null || true; } | awk '/^##/ { print > "header.txt"; next } /^#CHROM/ { n = split($0, a, "\t"); for (i = 1; i <= n; i++) print a[i]; exit }' > columns.txt || true; }
[ -s columns.txt ] || { echo "cannot read the header of $(name chr22)" >&2; exit 1; }
pairs=""
for s in "${SAMPLES[@]}"; do
  c=$(grep -n -x "$s" columns.txt | cut -d: -f1); [ -n "$c" ] || { echo "$s is not in the call set" >&2; exit 1; }
  pairs="$pairs$c $s"$'\n'
done
# cut writes the columns in the file's order, whatever order they are asked in: name them in that order
cols=$(printf '%s' "$pairs" | sort -n | awk '{ printf ",%s", $1 }')
NAMES=$(printf '%s' "$pairs" | sort -n | awk '{ printf "%s%s", (NR > 1 ? " " : ""), $2 }')

# the indexes first (tabix would fetch each on its own, once per process), then the regions in
# chunks of about 12 Mb, one connection per chunk
for chrom in $CHROMS; do [ -s "$(name "$chrom").tbi" ] || curl -s --retry 3 -o "$(name "$chrom").tbi" "$U/$(name "$chrom").tbi"; done
: > chunks.txt
for chrom in $CHROMS; do
  awk -F'\t' -v c="$chrom" -v dir="$T" '$1 == c { if (n == 0 || bp > 12000000) { n++; bp = 0 } print > (dir "/" c "." n ".bed"); bp += $3 - $2 } END { for (i = 1; i <= n; i++) print c "." i }' "$BED" >> chunks.txt
done
one() {   # <chrom.n>: the chunk's regions, the samples' columns, split per sample
  local id=$1 chrom=${1%%.*} try
  [ -e "$id.done" ] && return 0
  for try in 1 2 3 4; do
    rm -f "$id".*.part
    if tabix "$U/$(name "$chrom")" -R "$id.bed" 2> "$id.err" | cut -f "1-9$cols" | awk -F'\t' -v OFS='\t' -v id="$id" -v names="$NAMES" '
        BEGIN { n = split(names, nm, " ") }
        { info = "."
          if ($5 ~ /</) { info = ""; k = split($8, a, ";")
            for (j = 1; j <= k; j++) if (a[j] ~ /^(END|SVLEN|SVTYPE)=/) info = info (info == "" ? "" : ";") a[j]
            if (info == "") info = "." }
          for (i = 1; i <= n; i++) if ($(9 + i) ~ /[1-9]/) print $1, $2, $3, $4, $5, ".", ".", info, "GT", $(9 + i) > (id "." nm[i] ".part") }' \
       && [ ! -s "$id.err" ]; then
      : > "$id.done"; return 0
    fi
    sleep 5
  done
  echo "failed: $id ($(head -c 200 "$id.err"))" >&2; return 1
}
export -f one name; export U cols NAMES; set -o pipefail
xargs -P "$JOBS" -I{} bash -c 'set -o pipefail; one {}' < chunks.txt
[ "$(ls ./*.done | wc -l)" -eq "$(wc -l < chunks.txt)" ] || { echo "some chunks are missing; run again to resume" >&2; exit 1; }

for s in "${SAMPLES[@]}"; do
  { echo '##fileformat=VCFv4.2'
    grep '^##contig' header.txt
    echo '##ALT=<ID=DEL,Description="Deletion">'
    echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">'
    echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">'
    echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV Type">'
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$s"
    for chrom in $CHROMS; do   # a record that reaches into two regions is returned twice
      cat "$chrom".*."$s".part 2> /dev/null | sort -s -t$'\t' -k2,2n | awk '!seen[$0]++' || true
    done
  } | bgzip -c > "$OUT/$s.kgp.vcf.gz"
  tabix -f -p vcf "$OUT/$s.kgp.vcf.gz"
  # the same calls as an unphased file: what a call set without phasing gives
  bgzip -dc "$OUT/$s.kgp.vcf.gz" | awk -F'\t' -v OFS='\t' '/^#/ { print; next } { gsub(/\|/, "/", $10); print }' | bgzip -c > "$OUT/$s.kgp.unphased.vcf.gz"
  tabix -f -p vcf "$OUT/$s.kgp.unphased.vcf.gz"
  echo "$s: $(bgzip -dc "$OUT/$s.kgp.vcf.gz" | grep -vc '^#') records"
done
cp "$BED" "$OUT/regions.bed"
