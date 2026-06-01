#!/usr/bin/env bash
# SQANTI-SIM truth-based benchmark for PanIsoGuard (per-class P/R + AUPRC).
# Edit the paths/envs below, then: bash run.sh
# Tools: SQANTI-SIM (lean env, no R) + FLAIR 3.x + SQANTI3 5.x + the built panisoguard.
set -euo pipefail

# --- edit these ---------------------------------------------------------------
SQSIM=$HOME/SQANTI-SIM                 # SQANTI-SIM checkout
SQSIM_BIN=$HOME/miniconda3/envs/sqsim/bin   # lean env (biopython bcbio-gff pysam pbsim3 pbccs gffread minimap2 samtools scikit-learn pybedtools bx-python)
FLAIR_BIN=$HOME/miniconda3/envs/flair/bin
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin
PIG=$HOME/PanIsoGuard/build/panisoguard
HERE=$(cd "$(dirname "$0")" && pwd)
GENOME=/path/to/GRCh38.primary_assembly.genome.fa
FULL_GTF=/path/to/gencode.vNN.annotation.gtf
W=work; mkdir -p "$W"; cd "$W"
PB_MODEL=$SQSIM_BIN/../data/QSHMM-RSII.model

# --- chr22 inputs -------------------------------------------------------------
awk -F'\t' '$1=="chr22"' "$FULL_GTF" > chr22.gtf
samtools faidx "$GENOME" chr22 > chr22.fa && samtools faidx chr22.fa

# --- SQANTI-SIM: classify -> design (delete novel truth) -> simulate ----------
export PATH="$SQSIM_BIN:$PATH"
python "$SQSIM/sqanti-sim.py" classif --gtf chr22.gtf -o chr22 -d . --cores 8
python "$SQSIM/sqanti-sim.py" design equal -i chr22_index.tsv --gtf chr22.gtf \
    -o chr22 -d . -nt 1500 --NIC 300 --NNC 300 --ISM 200 --read_count 100000 -s 1 -k 8
python "$SQSIM/sqanti-sim.py" sim --gtf chr22.gtf --genome chr22.fa -i chr22_index.tsv \
    -d . --pb --pbsim --pbsim_model "$PB_MODEL" --long_count 100000 -k 8

# --- align -> FLAIR collapse (vs the REDUCED annotation) ----------------------
minimap2 -ax splice:hq -t 8 chr22.fa PBSIM3_simulated.fasta | samtools sort -o aln.bam - && samtools index aln.bam
python "$HERE/../end2end/bam2bed12.py" aln.bam reads.bed12
PATH="$FLAIR_BIN:$SQSIM_BIN:$PATH" flair collapse -g chr22.fa -q reads.bed12 \
    -r PBSIM3_simulated.fasta --gtf chr22_modified.gtf -o flair --threads 8

# --- SQANTI3 qc (vs reduced) --------------------------------------------------
PATH="$SQANTI3_BIN:$PATH" sqanti3_qc.py --isoforms flair.isoforms.gtf \
    --refGTF chr22_modified.gtf --refFasta chr22.fa -o flair -d sqanti_out --report skip --cpus 8

# --- truth SJ.tab + adjudicate + score ---------------------------------------
python "$HERE/prep_truth_sj.py" chr22.gtf chr22_modified.gtf PBSIM3_simulated.read_to_isoform.tsv truth
"$PIG" adjudicate --classification sqanti_out/flair_classification.txt \
    --isoforms-gtf flair.isoforms.gtf --ref-gtf chr22_modified.gtf --sj-tab truth.SJ.tab --out-prefix sqsim
python "$HERE/score.py" flair.isoforms.gtf truth.truth.tsv sqsim.adjudicated.tsv
