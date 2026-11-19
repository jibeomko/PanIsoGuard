#!/usr/bin/env bash
# Reference bias three ways (see README.md): the same Iso-Seq reads through the GRCh38
# pipeline (A), an SNV-personalized diploid (B1), the HPRC phased assembly (B2), and
# PanIsoGuard post hoc on A: with the two SNV-consensus FASTAs (C), with every call of the
# assembly, insertions and deletions included, as one haploid VCF per haplotype (C_vcf), and
# with a call set that does not come from the assembly, the 1000 Genomes short-read calls of
# the same person, phased and with the phase removed (C_1kgp).
# Builds on the hg03516_refbias outputs in $DATA.
#
# Usage: run.sh <stage>...   or   run.sh all
# Stages are idempotent (each returns early when its final output exists) and each is
# timed with /usr/bin/time -v into $WORK/time.<stage>.txt for the cost table.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
SAMPLE="${SAMPLE:-hg03516}"
DATA="${DATA:-/mnt/Data/PanIsoGuard_bench/pig_$SAMPLE}"            # dl.sh / geno.sh / rna.sh / isoquant.sh outputs
WORK="${WORK:-$DATA/refbias_3way}"
GRCH38="${GRCH38:-/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa}"
REFGTF="${REFGTF:-/mnt/Data/LR-RNAseq/reference/gencode.v49.annotation.gtf}"
THREADS="${THREADS:-16}"
CONDA="${CONDA:-$HOME/miniconda3/envs}"
MM2="${MM2:-$CONDA/flair300/bin/minimap2}"
SAMTOOLS="${SAMTOOLS:-$CONDA/sqanti3/bin/samtools}"
ISOQUANT="${ISOQUANT:-$CONDA/isoquant/bin/isoquant}"
LIFTOFF="${LIFTOFF:-$CONDA/liftoff/bin/liftoff}"
SQ3BIN="${SQ3BIN:-$CONDA/sqanti3/bin}"
BCFTOOLS="${BCFTOOLS:-$CONDA/pigval/bin/bcftools}"
KGP="${KGP:-$(dirname "$DATA")/kgp_1000g}"                          # fetch_1kgp.sh output: <ID>.kgp.vcf.gz, <ID>.kgp.unphased.vcf.gz
export PATH="$SQ3BIN:$PATH"   # sqanti3_qc.py's `#!/usr/bin/env python` must find its own env
ALL="A_sqanti B1_build B1_align B1_isoquant B1_sqanti B2_paf B2_liftoff B2_build B2_align B2_isoquant B2_sqanti C_adjudicate C_vcf C_1kgp"

mkdir -p "$WORK"/{armA,armB1,armB2,armC,armC_vcf,armC_1kgp,paf}; cd "$WORK"
log() { echo "[$(date '+%F %T')] $*"; }

# Outer invocation: run each stage in a child under /usr/bin/time so wall/RSS are recorded.
if [ -z "${_INNER:-}" ]; then
  [ "${1:-}" = all ] && set -- $ALL
  for s in "$@"; do
    log "== $s"
    if [ -x /usr/bin/time ]; then _INNER=1 /usr/bin/time -v -o "$WORK/time.$s.txt" bash "$HERE/run.sh" "$s"
    else _INNER=1 bash "$HERE/run.sh" "$s"; fi
  done
  exit 0
fi

# ---- helpers -------------------------------------------------------------------------
align() {  # <genome.fa> <out.bam>: flnc -> minimap2 splice:hq (one index part) -> sorted BAM
  $SAMTOOLS fastq -@4 "$DATA/flnc.bam" 2>/dev/null \
    | $MM2 -ax splice:hq -uf -t "$THREADS" -I 8G "$1" - 2> "$2.mm2.log" \
    | $SAMTOOLS sort -@4 -m 2G -o "$2" - && $SAMTOOLS index "$2"
}
isoquant() {  # <genome.fa> <genedb.gtf> <in.bam> <outdir>
  $ISOQUANT --reference "$1" --genedb "$2" --complete_genedb --bam "$3" --data_type pacbio_ccs \
    -o "$4" --threads 32 --prefix iq > "$4.log" 2>&1
}
sqanti() {  # <isoforms.gtf> <ref.gtf> <genome.fa> <outdir>: QC (no ORF, no report) + rules filter
  sqanti3_qc.py --isoforms "$1" --refGTF "$2" --refFasta "$3" -o iq -d "$4" --report skip -t 4 -n 4 \
    > "$4.log" 2>&1
  sqanti3_filter.py rules --sqanti_class "$4/iq_classification.txt" -o rules -d "$4/rules" \
    --skip_report > "$4.filter.log" 2>&1
}
suffix_ids() {  # <H>: append _<H> to gene_id / transcript_id (and to the chromosome if $2 = chrom)
  awk -v H="$1" -v C="${2:-}" 'BEGIN{FS=OFS="\t"} /^#/{next}
    { if (C == "chrom") $1 = $1 "_" H
      sub(/gene_id "[^"]+/, "&_" H); sub(/transcript_id "[^"]+/, "&_" H); print }'
}

# ---- A: GRCh38 pipeline (alignment + IsoQuant already in $DATA) ------------------------
A_sqanti() {
  [ -s armA/sq/rules/rules_RulesFilter_classification.txt ] && return
  sqanti "$DATA/isoquant_out/iq/iq.transcript_models.gtf" "$REFGTF" "$GRCH38" armA/sq
}

# ---- B1: SNV-personalized diploid, GRCh38 coordinates ---------------------------------
B1_build() {
  [ -s armB1/diploid.gtf ] && [ -s armB1/diploid.fa.fai ] && return
  for H in mat pat; do sed "s/^>\(\S*\).*/>\1_$H/" "$DATA/${H}_consensus.fa"; done > armB1/diploid.fa
  $SAMTOOLS faidx armB1/diploid.fa
  for H in mat pat; do suffix_ids $H chrom < "$REFGTF"; done > armB1/diploid.gtf
}
B1_align()    { [ -s armB1/rna.aln.bam.bai ] || align armB1/diploid.fa armB1/rna.aln.bam; }
B1_isoquant() { [ -s armB1/isoquant/iq/iq.transcript_models.gtf ] || isoquant armB1/diploid.fa armB1/diploid.gtf armB1/rna.aln.bam armB1/isoquant; }
B1_sqanti()   { [ -s armB1/sq/rules/rules_RulesFilter_classification.txt ] || sqanti armB1/isoquant/iq/iq.transcript_models.gtf armB1/diploid.gtf armB1/diploid.fa armB1/sq; }

# ---- B2: HPRC phased assembly diploid ------------------------------------------------
B2_paf() {  # assembly (query) -> GRCh38 (target), kept for paftools.js liftover
  for H in mat pat; do
    [ -s paf/$H.asm5.paf ] || $MM2 -cx asm5 --cs -t "$THREADS" "$GRCH38" "$DATA/$H.fa.gz" > paf/$H.asm5.paf 2> paf/$H.mm2.log
  done
}
B2_liftoff() {
  [ -s armB2/gencode_pat.gtf ] && return
  [ -s armB2/gencode.v49.gtf ] || cp "$REFGTF" armB2/gencode.v49.gtf   # liftoff writes its <gtf>_db next to the GTF
  for H in mat pat; do
    [ -s armB2/$H.fa ] || zcat "$DATA/$H.fa.gz" > armB2/$H.fa
    [ -s armB2/gencode_$H.gtf ] || $LIFTOFF -g armB2/gencode.v49.gtf -o armB2/gencode_$H.gtf -u armB2/unmapped_$H.txt \
      -p "$THREADS" -dir armB2/liftoff_$H armB2/$H.fa "$GRCH38" > armB2/liftoff_$H.log 2>&1
  done
}
B2_build() {
  [ -s armB2/diploid.gtf ] && [ -s armB2/diploid.fa.fai ] && return
  cat armB2/mat.fa armB2/pat.fa > armB2/diploid.fa      # contig names are already HG03516#1#/#2#-prefixed
  $SAMTOOLS faidx armB2/diploid.fa
  for H in mat pat; do suffix_ids $H < armB2/gencode_$H.gtf; done > armB2/diploid.gtf
}
B2_align()    { [ -s armB2/rna.aln.bam.bai ] || align armB2/diploid.fa armB2/rna.aln.bam; }
B2_isoquant() { [ -s armB2/isoquant/iq/iq.transcript_models.gtf ] || isoquant armB2/diploid.fa armB2/diploid.gtf armB2/rna.aln.bam armB2/isoquant; }
B2_sqanti() {  # liftoff writes CDS lines without frames and SQANTI3's gtfToGenePred rejects them: feed gene/transcript/exon only
  [ -s armB2/sq/rules/rules_RulesFilter_classification.txt ] && return
  awk -F'\t' '$3=="gene" || $3=="transcript" || $3=="exon"' armB2/diploid.gtf > armB2/diploid.sqanti.gtf
  sqanti armB2/isoquant/iq/iq.transcript_models.gtf armB2/diploid.sqanti.gtf armB2/diploid.fa armB2/sq
}

# ---- C: PanIsoGuard post hoc on arm A --------------------------------------------------
C_adjudicate() {
  [ -s armC/unk.adjudicated.tsv ] && return
  common=(--classification armA/sq/iq_classification.txt --isoforms-gtf "$DATA/isoquant_out/iq/iq.transcript_models.gtf"
          --ref-gtf "$REFGTF" --bam "$DATA/rna.aln.bam" --reference "$GRCH38")
  haps=(--reference-haplotype "$DATA/mat_consensus.fa" --reference-haplotype "$DATA/pat_consensus.fa")
  "$PIG" adjudicate "${common[@]}" --out-prefix armC/nohap > armC/nohap.log 2>&1
  "$PIG" adjudicate "${common[@]}" "${haps[@]}" --haplotype-provenance wgs     --out-prefix armC/wgs > armC/wgs.log 2>&1
  "$PIG" adjudicate "${common[@]}" "${haps[@]}" --haplotype-provenance unknown --out-prefix armC/unk > armC/unk.log 2>&1
}

# ---- C with the assembly's own calls: substitutions, insertions and deletions ----------
# $DATA/<H>.vcf is `paftools.js call` of one assembly haplotype against GRCh38 (geno.sh); it
# writes every genotype as 1/1, and one haplotype per file is GT 1.
C_vcf() {
  [ -s armC_vcf/wgs.adjudicated.tsv ] && return
  for H in mat pat; do
    [ -s "$DATA/$H.hap.vcf.gz.tbi" ] && continue
    "$BCFTOOLS" sort "$DATA/$H.vcf" 2> /dev/null | sed 's/\t1\/1$/\t1/' | bgzip -c > "$DATA/$H.hap.vcf.gz"
    tabix -f -p vcf "$DATA/$H.hap.vcf.gz"
  done
  "$PIG" adjudicate --classification armA/sq/iq_classification.txt --isoforms-gtf "$DATA/isoquant_out/iq/iq.transcript_models.gtf" \
    --ref-gtf "$REFGTF" --bam "$DATA/rna.aln.bam" --reference "$GRCH38" \
    --haplotype-vcf "$DATA/mat.hap.vcf.gz" --haplotype-vcf "$DATA/pat.hap.vcf.gz" --haplotype-provenance wgs \
    --out-prefix armC_vcf/wgs > armC_vcf/wgs.log 2>&1
}

# ---- C with a call set independent of the assembly: 1000 Genomes high-coverage short reads --
# The sample's own genotypes from the 3,202-sample phased panel (SNVs, indels, SVs), read in the
# regions of the isoforms that have a novel junction, and the same file with every genotype
# unphased. One diploid VCF, so the two haplotypes are the two alleles of each genotype.
C_1kgp() {
  [ -s armC_1kgp/unphased.adjudicated.tsv ] && return
  ID=$(echo "$SAMPLE" | tr a-z A-Z)
  if [ ! -s "$KGP/$ID.kgp.vcf.gz.tbi" ]; then
    python3 "$HERE/novel_isoform_regions.py" armC/wgs.junctions.tsv "$DATA/isoquant_out/iq/iq.transcript_models.gtf" > armC_1kgp/regions.bed
    bash "$HERE/fetch_1kgp.sh" "$KGP" armC_1kgp/regions.bed "$ID"
  fi
  for v in phased unphased; do
    f="$KGP/$ID.kgp.vcf.gz"; [ $v = unphased ] && f="$KGP/$ID.kgp.unphased.vcf.gz"
    "$PIG" adjudicate --classification armA/sq/iq_classification.txt --isoforms-gtf "$DATA/isoquant_out/iq/iq.transcript_models.gtf" \
      --ref-gtf "$REFGTF" --bam "$DATA/rna.aln.bam" --reference "$GRCH38" \
      --haplotype-vcf "$f" --haplotype-provenance wgs --out-prefix armC_1kgp/$v > armC_1kgp/$v.log 2>&1
  done
}

for s in "$@"; do "$s"; done
