#!/usr/bin/env bash
# LRGASP WTC11: PanIsoGuard against experimentally validated isoforms on real PacBio cDNA,
# ONT cDNA and Illumina data (ENCODE, public). Stages are idempotent; each is timed with
# /usr/bin/time -v into $WORK/time.<stage>.txt.
#   reads   ENCSR507JOF (PacBio cDNA, 3 reps) ENCSR539ZXJ (ONT cDNA, rep 1) ENCSR673UKZ (Illumina, rep 1)
#   truth   LRGASP Supplementary Data 15 (WTC11 validation analysis) + the WTC11 consolidated models BED
# Usage: run.sh <stage>...  or  run.sh all   (after raw/dl.sh has downloaded the reads)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/mnt/Data/pig_lrgasp}"
GRCH38="${GRCH38:-/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa}"
REFGTF="${REFGTF:-/mnt/Data/LR-RNAseq/reference/gencode.v49.annotation.gtf}"
STAR_IDX="${STAR_IDX:-/mnt/Data/LR-RNAseq/reference/STAR_index_v49}"
THREADS="${THREADS:-16}"
CONDA="${CONDA:-$HOME/miniconda3/envs}"
MM2="${MM2:-$CONDA/flair300/bin/minimap2}"
SAMTOOLS="${SAMTOOLS:-$CONDA/sqanti3/bin/samtools}"
STAR="${STAR:-$CONDA/star274/bin/STAR}"   # must be the STAR version that built $STAR_IDX (a newer index has parameters 2.7.4a rejects)
ISOQUANT="${ISOQUANT:-$CONDA/isoquant/bin/isoquant}"
SQ3BIN="${SQ3BIN:-$CONDA/sqanti3/bin}"
export PATH="$SQ3BIN:$PATH"
ALL="targets align_pb align_ont star isoquant_pb isoquant_ont sqanti_pb sqanti_ont adjudicate score"
mkdir -p "$WORK"/{pb,ont,star,targets}; cd "$WORK"
log() { echo "[$(date '+%F %T')] $*"; }

if [ -z "${_INNER:-}" ]; then
  [ "${1:-}" = all ] && set -- $ALL
  for s in "$@"; do
    log "== $s"
    if [ -x /usr/bin/time ]; then _INNER=1 /usr/bin/time -v -o "$WORK/time.$s.txt" bash "$HERE/run.sh" "$s"
    else _INNER=1 bash "$HERE/run.sh" "$s"; fi
  done
  exit 0
fi

targets() { [ -s targets/targets.tsv ] || python3 "$HERE/targets.py" supp/wtc11_validation_targets.tsv models/WTC11_consolidated.bed.gz targets; }
align_pb() {
  [ -s pb/pb.bam.bai ] && return
  zcat raw/pb_cdna_rep1.fastq.gz raw/pb_cdna_rep2.fastq.gz raw/pb_cdna_rep3.fastq.gz \
    | $MM2 -ax splice:hq -uf --MD -t "$THREADS" "$GRCH38" - 2> pb/mm2.log \
    | $SAMTOOLS sort -@4 -m 2G -o pb/pb.bam - && $SAMTOOLS index pb/pb.bam
}
align_ont() {
  [ -s ont/ont.bam.bai ] && return
  zcat raw/ont_cdna_rep1.fastq.gz \
    | $MM2 -ax splice --MD -t "$THREADS" "$GRCH38" - 2> ont/mm2.log \
    | $SAMTOOLS sort -@4 -m 2G -o ont/ont.bam - && $SAMTOOLS index ont/ont.bam
}
star() {
  [ -s star/SJ.out.tab ] && return
  $STAR --genomeDir "$STAR_IDX" --readFilesIn raw/illumina_rep1_R1.fastq.gz raw/illumina_rep1_R2.fastq.gz \
    --readFilesCommand zcat --twopassMode Basic --outSAMtype None --runThreadN "$THREADS" --outFileNamePrefix star/ > star/star.log 2>&1
}
isoquant_pb()  { [ -s pb/isoquant/iq/iq.transcript_models.gtf ] || $ISOQUANT --reference "$GRCH38" --genedb "$REFGTF" --complete_genedb --bam pb/pb.bam --data_type pacbio_ccs -o pb/isoquant --threads 32 --prefix iq > pb/isoquant.log 2>&1; }
isoquant_ont() { [ -s ont/isoquant/iq/iq.transcript_models.gtf ] || $ISOQUANT --reference "$GRCH38" --genedb "$REFGTF" --complete_genedb --bam ont/ont.bam --data_type nanopore -o ont/isoquant --threads 32 --prefix iq > ont/isoquant.log 2>&1; }
sqanti() {  # <isoforms.gtf> <outdir>
  [ -s "$2/iq_classification.txt" ] && return
  sqanti3_qc.py --isoforms "$1" --refGTF "$REFGTF" --refFasta "$GRCH38" -o iq -d "$2" --report skip -t 4 -n 4 > "$2.log" 2>&1
}
sqanti_pb()  { sqanti pb/isoquant/iq/iq.transcript_models.gtf pb/sq; }
sqanti_ont() { sqanti ont/isoquant/iq/iq.transcript_models.gtf ont/sq; }
adjudicate() {
  for p in pb ont; do
    preset=$([ $p = pb ] && echo hifi || echo ont)
    common=(--classification $p/sq/iq_classification.txt --isoforms-gtf $p/isoquant/iq/iq.transcript_models.gtf --ref-gtf "$REFGTF"
            --bam $p/$p.bam --reference "$GRCH38" --platform $preset)
    [ -s $p/adj_sr.adjudicated.tsv ]   || "$PIG" adjudicate "${common[@]}" --sj-tab star/SJ.out.tab --out-prefix $p/adj_sr   > $p/adj_sr.log 2>&1
    [ -s $p/adj_nosr.adjudicated.tsv ] || "$PIG" adjudicate "${common[@]}"                       --out-prefix $p/adj_nosr > $p/adj_nosr.log 2>&1
  done
}
score() { python3 "$HERE/score_lrgasp.py" --work "$WORK" --emit-metrics "$HERE/../results/lrgasp_wtc11/metrics.json"; }

for s in "$@"; do "$s"; done
