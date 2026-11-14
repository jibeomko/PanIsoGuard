#!/usr/bin/env bash
# --no-sqanti3 against SQANTI3 on the chr22 SQANTI-SIM five-caller set: every caller's output is
# adjudicated twice, from its SQANTI3 classification and from the computed one, in the two
# settings of the main text (no short reads + BAM; 1 M pairs + BAM). score_no_sqanti3.py compares
# the categories, the QC priors and the verdicts, and scores both modes against the truth.
# Env: WORK (sqanti_sim work dir), MC (multicaller data dir), GENOME (chr22.fa), SJ (1 M SJ.out.tab), OUT.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/sqanti_sim/work}"
MC="${MC:-/path/to/multicaller/data}"
GENOME="${GENOME:-/path/to/chr22.fa}"
SJ="${SJ:-/path/to/h2h_sr/sr_1000000/SJ.out.tab}"
OUT="${OUT:-/tmp/no_sqanti3}"; mkdir -p "$OUT"
declare -A GTF=([flair]="$WORK/flair.isoforms.gtf" [isoquant]="$MC/isoquant_out/iq/iq.transcript_models.gtf"
                [bambu]="$MC/bambu.novel.gtf" [espresso]="$MC/espresso_out/espresso_samples_N2_R0_updated.gtf"
                [talon]="$MC/talon_out/talon_filt_talon.gtf")
declare -A CLS=([flair]="$WORK/sqanti_out/flair_classification.txt")
for c in isoquant bambu espresso talon; do CLS[$c]="$MC/sq_$c/${c}_classification.txt"; done
for c in flair isoquant bambu espresso talon; do
  common=(--isoforms-gtf "${GTF[$c]}" --ref-gtf "$WORK/chr22_modified.gtf" --bam "$WORK/aln.bam" --reference "$GENOME")
  for s in nosr sr1M; do
    sj=(); [ $s = sr1M ] && sj=(--sj-tab "$SJ")
    "$PIG" adjudicate --classification "${CLS[$c]}" "${common[@]}" "${sj[@]}" --out-prefix "$OUT/$c.$s.sqanti3" > /dev/null 2>&1
    "$PIG" adjudicate --no-sqanti3 "${common[@]}" "${sj[@]}" --out-prefix "$OUT/$c.$s.computed" > /dev/null 2>&1
  done
done
python3 "$HERE/score_no_sqanti3.py" --out "$OUT" --truth "$WORK/truth.truth.tsv" \
  --gtf flair="${GTF[flair]}" --gtf isoquant="${GTF[isoquant]}" --gtf bambu="${GTF[bambu]}" --gtf espresso="${GTF[espresso]}" --gtf talon="${GTF[talon]}" \
  --cls flair="${CLS[flair]}" --cls isoquant="${CLS[isoquant]}" --cls bambu="${CLS[bambu]}" --cls espresso="${CLS[espresso]}" --cls talon="${CLS[talon]}" \
  --tool-version "$("$PIG" --version | head -1 | awk '{print $2}')" --emit-metrics "$HERE/../results/no_sqanti3/metrics.json"
