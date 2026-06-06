#!/usr/bin/env bash
# Reference-bias rescue on real divergent-individual long-read RNA (HG03516, HPRC R2, ESN/AFR).
# All inputs are public (HPRC open S3; GENCODE v49; GRCh38). Heavy: ~14GB download +
# whole-genome asm5 + 11.3M-read splice alignment + IsoQuant (multi-hour).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/hg03516}"; mkdir -p "$WORK"; cd "$WORK"
GRCH38="${GRCH38:-/path/to/GRCh38.primary_assembly.genome.fa}"   # + .fai
REFGTF="${REFGTF:-/path/to/gencode.v49.annotation.gtf}"
SCAN="$HERE/../hg002/scan_variant_axis.py"
ISOQUANT="${ISOQUANT:-isoquant}"; BCF="${BCF:-bcftools}"; TBX="${TBX:-tabix}"
S3=s3://human-pangenomics/working/HPRC/HG03516

# 1) download Iso-Seq flnc + phased assembly (public, no-sign-request)
aws s3 cp --no-sign-request "$S3/raw_data/PacBio_Kinnex/HG03516.lymph.m84081_240703_203142_s1-m84081_240709_185528_s2.flnc.bam" flnc.bam
aws s3 cp --no-sign-request "$S3/assemblies/release2/HG03516_pat_hprc_r2_v1.1.0.fa.gz" pat.fa.gz
aws s3 cp --no-sign-request "$S3/assemblies/release2/HG03516_mat_hprc_r2_v1.1.0.fa.gz" mat.fa.gz

# 2) personal SNV-consensus haplotypes (assembly -> GRCh38 variants -> length-matched consensus)
for H in pat mat; do
  minimap2 -cx asm5 --cs -t 16 "$GRCH38" $H.fa.gz | sort -k6,6 -k8,8n -S4G | \
    paftools.js call -f "$GRCH38" -s $H - > $H.vcf
  "$BCF" view -v snps $H.vcf | "$BCF" sort -Oz -o $H.snv.vcf.gz; "$TBX" -p vcf $H.snv.vcf.gz
  "$BCF" consensus -f "$GRCH38" $H.snv.vcf.gz > ${H}_consensus.fa; samtools faidx ${H}_consensus.fa
done

# 3) RNA: flnc -> splice alignment -> IsoQuant novel transcript models
samtools fastq -@8 flnc.bam | minimap2 -ax splice:hq -uf -t16 "$GRCH38" - | \
  samtools sort -@8 -m2G -o rna.aln.bam - && samtools index rna.aln.bam
"$ISOQUANT" --reference "$GRCH38" --genedb "$REFGTF" --complete_genedb \
  --bam rna.aln.bam --data_type pacbio_ccs -o isoquant_out --threads 32 --prefix iq

# 4) reference-bias scan of the novel junctions vs BOTH personal haplotypes + adjudicate
python3 - "$WORK" <<'PY'
import sys, collections
W=sys.argv[1]; src=f"{W}/isoquant_out/iq/iq.transcript_models.gtf"
def tid(a): return a.split('transcript_id "')[1].split('"')[0] if 'transcript_id "' in a else ""
by=collections.defaultdict(list)
for ln in open(src):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)<9 or f[2]!='exon': continue
    if '.nic' in tid(f[8]) or '.nnic' in tid(f[8]): by[f[0]].append(ln)
import os; os.makedirs(f"{W}/scan",exist_ok=True)
for ch,ls in by.items(): open(f"{W}/scan/{ch}.novel.gtf","w").writelines(ls)
PY
for chrom in chr{1..22} chrX; do
  g="$WORK/scan/$chrom.novel.gtf"; [ -f "$g" ] || continue
  samtools faidx "$GRCH38" $chrom > "$WORK/scan/$chrom.ref.fa"
  for H in pat mat; do
    samtools faidx ${H}_consensus.fa $chrom > "$WORK/scan/$chrom.$H.fa"
    d="$WORK/scan/$chrom.$H"; mkdir -p "$d"
    python3 "$SCAN" "$WORK/scan/$chrom.ref.fa" "$WORK/scan/$chrom.$H.fa" "$g" $chrom "$d" >/dev/null
    for prov in wgs unknown; do
      tag=$([ "$prov" = wgs ] && echo wgs || echo unk)
      "$PIG" adjudicate --classification "$d/vt.cls" --isoforms-gtf "$d/vt.gtf" \
        --ref-gtf "$d/vt.catalog.gtf" --reference "$WORK/scan/$chrom.ref.fa" \
        --reference-haplotype "$WORK/scan/$chrom.$H.fa" --haplotype-provenance $prov \
        --out-prefix "$d/adj_$tag" >/dev/null
    done
  done
done

# 5) dedup CREATED hits across pat/mat -> refbias_unique.tsv, then emit metrics
python3 "$HERE/dedup_hits.py" "$WORK"
SAMTOOLS=samtools python3 "$HERE/analyze.py" "$WORK" "$REFGTF" "$WORK/rna.aln.bam" \
  --emit-metrics "$HERE/../results/hg03516_refbias/metrics.json"
