#!/bin/bash
# End-to-end read-simulation E1: simulate -> align -> FLAIR -> SQANTI3 -> PanIsoGuard,
# scored against intron-chain truth (hidden GENCODE transcripts = genuine novels).
#
# This is the realistic, messy complement to the clean controlled_truth / synthetic_axes
# checks: it validates that PanIsoGuard behaves correctly on REAL caller (FLAIR) +
# SQANTI3 output, not hand-constructed inputs.
#
# Tools (each typically in its own conda env -- adjust the paths below):
#   pbsim3, minimap2, samtools, FLAIR 3.x, SQANTI3 5.x, and the built panisoguard.
# NOTE: FLAIR 3.0.0's `bam2Bed12` entry point is broken on bioconda
#   (ModuleNotFoundError: flair.bam2Bed12); bam2bed12.py here replaces it.
set -euo pipefail

# ---- config (edit, or override via environment) -----------------------------
GENOME=${GENOME:-/path/to/GRCh38.primary_assembly.genome.fa}     # + .fai
GTF=${GTF:-/path/to/gencode.vNN.annotation.gtf}
CHR=chr22
PBSIM=$HOME/miniconda3/envs/pigval/bin/pbsim
PBSIM_MODEL=$HOME/miniconda3/envs/pigval/data/ERRHMM-SEQUEL.model
FLAIR_BIN=$HOME/miniconda3/envs/flair/bin            # has flair (bam2Bed12 is broken; we bypass)
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin        # has sqanti3_qc.py
PIG=$(cd "$(dirname "$0")/../../build" && pwd)/panisoguard
OUT=${1:-./e2e_work}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"; cd "$OUT"

# ---- 1) reference subset ----------------------------------------------------
awk -F'\t' -v c="$CHR" '$1==c' "$GTF" > chr.gtf
samtools faidx "$GENOME" "$CHR" > chr.fa; samtools faidx chr.fa

# ---- 2) labelled truth + sim input (hide ~20% transcripts) ------------------
python "$HERE/prep.py" chr.gtf chr.fa .

# ---- 3) simulate reads (PBSIM3 transcriptome) -------------------------------
"$PBSIM" --strategy trans --method errhmm --errhmm "$PBSIM_MODEL" \
         --transcript sim.transcript --prefix sim --seed 7 >/dev/null 2>&1
cat sim*.fq.gz > reads.fq.gz 2>/dev/null || { cat sim*.fastq | gzip > reads.fq.gz; }

# ---- 4) align (minimap2 splice) ---------------------------------------------
minimap2 -ax splice:hq -uf -t 8 chr.fa reads.fq.gz 2>/dev/null | samtools sort -@8 -o aln.bam -
samtools index aln.bam

# ---- 5) FLAIR collapse (BAM->bed12 via our converter, then collapse) --------
python "$HERE/bam2bed12.py" aln.bam reads.bed12
PATH="$FLAIR_BIN:$PATH" flair collapse -g chr.fa -q reads.bed12 -r reads.fq.gz \
     --gtf reduced.gtf -o flair --threads 8 >collapse.log 2>&1

# ---- 6) SQANTI3 QC vs the reduced annotation --------------------------------
PATH="$SQANTI3_BIN:$PATH" sqanti3_qc.py --isoforms flair.isoforms.gtf \
     --refGTF reduced.gtf --refFasta chr.fa -o e2e -d sqout --cpus 8 --report skip >sq.log 2>&1

# ---- 7) PanIsoGuard adjudicate ----------------------------------------------
"$PIG" adjudicate --classification sqout/e2e_classification.txt \
     --isoforms-gtf flair.isoforms.gtf --ref-gtf reduced.gtf \
     --sj-tab real.SJ.tab --bam aln.bam --reference chr.fa --out-prefix e2e_adj

# ---- 8) score vs intron-chain truth -----------------------------------------
python "$HERE/score.py"
