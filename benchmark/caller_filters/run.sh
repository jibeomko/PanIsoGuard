#!/usr/bin/env bash
# Post-caller filters on the same caller output (none / SQANTI3 rules / PanIsoGuard + BAM), five
# callers on the chr22 SQANTI-SIM alignment, and each caller's own knob (Bambu NDR, TALON's
# filter) against PanIsoGuard on the permissive output. Reuses the PanIsoGuard runs of
# ../chain_support (OUT there) and the caller outputs of ../multicaller (MC).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/sqanti_sim/work}"            # aln.bam, chr22_modified.gtf, truth.truth.tsv, flair.isoforms.gtf, sqanti_out/
MC="${MC:-/path/to/multicaller/data}"               # the other callers' GTFs + sq_<caller>/ SQANTI3 runs (../multicaller/run.sh)
CHAIN="${CHAIN:-/tmp/chain_support}"                # ../chain_support/run.sh OUT: <caller>.{nosr,sr1M}.on.adjudicated.tsv
H2H="${H2H:-/path/to/h2h_work}"                     # ../sqanti3_filter_h2h/run.sh OUT: none/rules/, sr_1000000/SJ.out.tab
GENOME="${GENOME:-/path/to/chr22.fa}"
OUT="${OUT:-/tmp/caller_filters}"; mkdir -p "$OUT"
export PATH="${SQANTI3_BIN:-$HOME/miniconda3/envs/sqanti3/bin}:$PATH"
RSCRIPT="${RSCRIPT:-$HOME/miniconda3/envs/bambu/bin/Rscript}"

# 1. SQANTI3's rules filter on each caller's classification (no coverage column: perc_A + RTS rules)
for c in isoquant bambu espresso talon; do
  sqanti3_filter.py rules --sqanti_class "$MC/sq_$c/${c}_classification.txt" -o rules -d "$OUT/rules_$c" --skip_report > "$OUT/rules_$c.log" 2>&1
done
# 2. TALON before talon_filter_transcripts (the observed-only GTF): SQANTI3, its rules filter, PanIsoGuard
sqanti3_qc.py --isoforms "$MC/talon_out/talon_talon_observedOnly.gtf" --refGTF "$WORK/chr22_modified.gtf" --refFasta "$GENOME" \
  -o talon_unfilt -d "$OUT/talon_unfilt" --report skip --cpus 8 > "$OUT/talon_unfilt.log" 2>&1
sqanti3_filter.py rules --sqanti_class "$OUT/talon_unfilt/talon_unfilt_classification.txt" -o rules -d "$OUT/rules_talon_unfilt" --skip_report > /dev/null 2>&1
common=(--classification "$OUT/talon_unfilt/talon_unfilt_classification.txt" --isoforms-gtf "$MC/talon_out/talon_talon_observedOnly.gtf"
        --ref-gtf "$WORK/chr22_modified.gtf" --bam "$WORK/aln.bam" --reference "$GENOME")
"$PIG" adjudicate "${common[@]}" --out-prefix "$CHAIN/talon_unfilt.nosr.on" > /dev/null 2>&1
"$PIG" adjudicate "${common[@]}" --sj-tab "$H2H/sr_1000000/SJ.out.tab" --out-prefix "$CHAIN/talon_unfilt.sr1M.on" > /dev/null 2>&1
# 3. Bambu at stricter NDR (the tracked set is NDR 0.5; ../multicaller/run_bambu_ndr.R)
for ndr in 0.3 0.2 0.1; do
  d="$OUT/bambu_ndr$ndr"; mkdir -p "$d"
  cat > "$d/run.R" <<R
suppressMessages(library(bambu))
ann <- prepareAnnotations("$WORK/chr22_modified.gtf")
se <- bambu(reads = "$WORK/aln.bam", annotations = ann, genome = "$GENOME", ncore = 4, NDR = $ndr, quant = FALSE)
writeToGTF(se, "$d/extended_annotations.gtf")
R
  "$RSCRIPT" "$d/run.R" > "$d/run.log" 2>&1
  awk -F'\t' '$9 ~ /transcript_id "BambuTx/' "$d/extended_annotations.gtf" > "$d/bambu.novel.gtf"
  sqanti3_qc.py --isoforms "$d/bambu.novel.gtf" --refGTF "$WORK/chr22_modified.gtf" --refFasta "$GENOME" -o bambu -d "$d/sq" --report skip --cpus 8 > "$d/sq.log" 2>&1
done
# 4. score
"${PY:-python3}" "$HERE/score_filters.py" --work "$WORK" --pig-dir "$CHAIN" \
  --caller "flair=$WORK/flair.isoforms.gtf,$WORK/sqanti_out/flair_classification.txt,$H2H/none/rules/rules_RulesFilter_classification.txt" \
  --caller "isoquant=$MC/isoquant_out/iq/iq.transcript_models.gtf,$MC/sq_isoquant/isoquant_classification.txt,$OUT/rules_isoquant/rules_RulesFilter_classification.txt" \
  --caller "bambu=$MC/bambu.novel.gtf,$MC/sq_bambu/bambu_classification.txt,$OUT/rules_bambu/rules_RulesFilter_classification.txt" \
  --caller "espresso=$MC/espresso_out/espresso_samples_N2_R0_updated.gtf,$MC/sq_espresso/espresso_classification.txt,$OUT/rules_espresso/rules_RulesFilter_classification.txt" \
  --caller "talon=$MC/talon_out/talon_filt_talon.gtf,$MC/sq_talon/talon_classification.txt,$OUT/rules_talon/rules_RulesFilter_classification.txt" \
  --family "bambu=bambu=$MC/bambu.novel.gtf,$MC/sq_bambu/bambu_classification.txt,$OUT/rules_bambu/rules_RulesFilter_classification.txt" \
  --family "bambu=bambu_ndr0.3=$OUT/bambu_ndr0.3/bambu.novel.gtf,$OUT/bambu_ndr0.3/sq/bambu_classification.txt" \
  --family "bambu=bambu_ndr0.2=$OUT/bambu_ndr0.2/bambu.novel.gtf,$OUT/bambu_ndr0.2/sq/bambu_classification.txt" \
  --family "bambu=bambu_ndr0.1=$OUT/bambu_ndr0.1/bambu.novel.gtf,$OUT/bambu_ndr0.1/sq/bambu_classification.txt" \
  --family "talon=talon_unfilt=$MC/talon_out/talon_talon_observedOnly.gtf,$OUT/talon_unfilt/talon_unfilt_classification.txt,$OUT/rules_talon_unfilt/rules_RulesFilter_classification.txt" \
  --family "talon=talon=$MC/talon_out/talon_filt_talon.gtf,$MC/sq_talon/talon_classification.txt,$OUT/rules_talon/rules_RulesFilter_classification.txt" \
  --emit-metrics "$HERE/../results/caller_filters/metrics.json"
