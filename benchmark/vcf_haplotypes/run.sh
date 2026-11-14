#!/usr/bin/env bash
# The variant axis from a VCF (--haplotype-vcf) against the same axis from consensus FASTAs
# (--reference-haplotype) on the runs that already exist on disk: HG002 chr22 (benchmark/hg002; GIAB
# biallelic SNVs, hap.fa = bcftools consensus -H 1) and every per-chromosome scan of HG03516 and HG02717
# (benchmark/hg03516_refbias; paftools SNV calls per parental assembly, consensus per haplotype). The
# FASTA-route outputs are the ones those protocols wrote (vt.*, adj_wgs.*); only the VCF route is run
# here, into vtvcf.* / adjvcf_wgs.*. compare_routes.py then counts isoforms with the same verdict and
# junctions with the same variant call and the same haplotype motif.
# Usage: bash run.sh      (env: B = the benchmark data root)
set -euo pipefail
B=${B:-/mnt/Data/PanIsoGuard_bench}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# HG002 chr22: one haplotype FASTA (-H 1) vs the VCF's two haplotypes; the first is the FASTA's
d=$B/pig_hg002
"$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" --ref-gtf "$d/vt.catalog.gtf" --reference "$d/chr22.ref.fa" \
  --haplotype-vcf "$d/chr22.snps.vcf.gz" --haplotype-provenance wgs --out-prefix "$d/vtvcf" > "$d/vtvcf.log" 2>&1
pairs=("HG002:$d/vt:$d/vtvcf")

# HG03516 / HG02717: per chromosome and parental haplotype, consensus FASTA vs the haplotype's SNV VCF
for s in hg03516 hg02717; do
  W=$B/pig_$s
  for chrom in chr{1..22} chrX; do for H in pat mat; do
    sd="$W/scan/$chrom.$H"; [ -f "$sd/vt.cls" ] && [ -s "$sd/adj_wgs.adjudicated.tsv" ] || continue
    "$PIG" adjudicate --classification "$sd/vt.cls" --isoforms-gtf "$sd/vt.gtf" --ref-gtf "$sd/vt.catalog.gtf" --reference "$W/scan/$chrom.ref.fa" \
      --haplotype-vcf "$W/$H.snv.vcf.gz" --haplotype-provenance wgs --out-prefix "$sd/adjvcf_wgs" > "$sd/adjvcf_wgs.log" 2>&1
    pairs+=("$(echo $s | tr a-z A-Z):$sd/adj_wgs:$sd/adjvcf_wgs")
  done; done
  echo "[$(date +%T)] $s: VCF route done"
done
python3 "$HERE/compare_routes.py" "${pairs[@]}" --emit-metrics "$HERE/../results/vcf_haplotypes/metrics.json"
