#!/usr/bin/env bash
# LRGASP WTC11 validation targets as would-be isoforms: every amplicon-tested intron chain of
# Supplementary Data 15 (both target sets of ../lrgasp_wtc11/targets.py), classified by SQANTI3 against
# GENCODE v49 and adjudicated by PanIsoGuard with the PacBio and the ONT library, with and without the
# Illumina junctions, independently of whether a caller reconstructed the chain. Negatives = chains the
# amplicons did not support; positives = the supported ones. Also SQANTI3's default rules filter.
# Usage: bash run.sh [/mnt/Data/pig_lrgasp]
set -euo pipefail
L=${1:-/mnt/Data/pig_lrgasp}; W=$L/targets_run; mkdir -p "$W"; cd "$W"
GRCH38=${GRCH38:-/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa}
REFGTF=${REFGTF:-/mnt/Data/LR-RNAseq/reference/gencode.v49.annotation.gtf}
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin; PIG=$HOME/PanIsoGuard/build/panisoguard
echo "[$(date +%T)] targets GTF"
[ -s targets_all.gtf ] || cat "$L/targets/targets.gtf" "$L/targets/targets_gencode.gtf" > targets_all.gtf
echo "[$(date +%T)] SQANTI3 qc vs GENCODE v49"
[ -s sq/tg_classification.txt ] || PATH="$SQANTI3_BIN:$PATH" sqanti3_qc.py --isoforms targets_all.gtf --refGTF "$REFGTF" --refFasta "$GRCH38" -o tg -d sq --report skip -t 4 -n 4 > sq.log 2>&1
[ -s rules/rules_RulesFilter_classification.txt ] || PATH="$SQANTI3_BIN:$PATH" sqanti3_filter.py rules --sqanti_class sq/tg_classification.txt -o rules -d rules --skip_report > rules.log 2>&1
for p in pb ont; do
  for s in nosr sr; do
    [ -s adj_${p}_$s.adjudicated.tsv ] && continue
    echo "[$(date +%T)] adjudicate $p $s"
    args=(adjudicate --classification sq/tg_classification.txt --isoforms-gtf targets_all.gtf --ref-gtf "$REFGTF" --bam "$L/$p/$p.bam" --reference "$GRCH38")
    [ $p = ont ] && args+=(--platform ont)
    [ $s = sr ] && args+=(--sj-tab "$L/star/SJ.out.tab")
    "$PIG" "${args[@]}" --out-prefix adj_${p}_$s > adj_${p}_$s.log 2>&1
  done
done
echo "[$(date +%T)] done: $W"
