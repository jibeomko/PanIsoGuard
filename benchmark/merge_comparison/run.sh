#!/usr/bin/env bash
# Merge-tool head-to-head: PanIsoGuard `combine` (exact intron-chain) vs gffcompare -i
# (incumbent N-way), TAMA merge (fuzzy), and a controlled exact->fuzzy wobble sweep, on
# the SAME 5 chr22 caller GTFs as the multicaller protocol. Answers two reviewer
# objections: (1) does combine just re-derive gffcompare? (2) is the consensus precision
# an artifact of exact matching? See README.md. All inputs public/simulated; no private data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# --- placeholders (fill in; reuse the multicaller protocol's outputs) -------
WORK=/path/to/sqanti_sim/work            # truth.truth.tsv, chr22_modified.gtf
MULTI=/path/to/pig_multicaller           # matrix5.tsv + the 5 caller GTFs
MERGECMP_BIN=$HOME/miniconda3/envs/mergecmp/bin     # gffcompare + ucsc-gtftogenepred + ucsc-genepredtobed
TAMA=/path/to/tama                       # github.com/GenomeRIK/tama (python2 + biopython)
PY2=python2                              # a python2 with numpy + biopython
# ----------------------------------------------------------------------------
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT/bed"; cd "$OUT"

# caller -> GTF (same order everywhere)
cat > callers.tsv <<EOF
flair	$WORK/flair.isoforms.gtf
isoquant	$MULTI/isoquant_out/iq/iq.transcript_models.gtf
bambu	$MULTI/bambu.novel.gtf
espresso	$MULTI/espresso_out/espresso_samples_N2_R0_updated.gtf
talon	$MULTI/talon_out/talon_filt_talon.gtf
EOF
cut -f2 callers.tsv > gtflist.txt

# 1) gffcompare -i (incumbent N-way comparison, no reference)
"$MERGECMP_BIN/gffcompare" -i gtflist.txt -o gffcmp

# 2) TAMA merge (fuzzy, junction wobble 10bp). Needs BED12 with "gene_id;transcript_id" col4.
: > filelist.txt
while IFS=$'\t' read -r name gtf; do
  "$MERGECMP_BIN/gtfToGenePred" -ignoreGroupsWithoutExons "$gtf" "bed/$name.gp"
  "$MERGECMP_BIN/genePredToBed" "bed/$name.gp" "bed/$name.raw.bed"
  awk -F'\t' 'BEGIN{OFS="\t"} {$4=$4";"$4; print}' "bed/$name.raw.bed" > "bed/$name.bed"
  printf '%s\tno_cap\t1,1,1\t%s\n' "$OUT/bed/$name.bed" "$name" >> filelist.txt
done < callers.tsv
"$PY2" "$TAMA/tama_merge.py" -f filelist.txt -p tama -a 50 -m 10 -z 50 -d merge_dup

# 3) Score all matchers on one footing (+ the wobble sweep, computed in-script).
python3 "$HERE/score_merge.py" --callers callers.tsv \
  --truth "$WORK/truth.truth.tsv" --ref "$WORK/chr22_modified.gtf" \
  --paniso "$MULTI/matrix5.tsv" --gffcompare gffcmp.tracking --tama tama_trans_report.txt \
  --wobble 0,5,10,20 --tool-version "$("$PIG" --version | head -1 | awk '{print $2}')" \
  --emit-metrics "$HERE/../results/merge_comparison/metrics.json"
