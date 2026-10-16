#!/usr/bin/env bash
# Whole-chain long-read support (ruleset 0.0.4): does counting the reads that carry an
# isoform's entire intron chain separate genuine from false novel isoforms, and for which
# callers? Five callers on the same chr22 SQANTI-SIM alignment (the benchmark/multicaller
# set), each adjudicated with the chain axis ON (shipped default, min_full_reads = 3) and
# OFF (min_full_reads raised beyond any read count), with 1 M simulated short-read pairs
# and without short reads. Scored against the SQANTI-SIM truth by score_chain.py.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/sqanti_sim/work}"          # aln.bam, chr22_modified.gtf, truth.truth.tsv, flair.isoforms.gtf
MC="${MC:-/path/to/multicaller/data}"             # the other callers' GTFs + sq_<caller>/ SQANTI3 runs (see ../multicaller/run.sh)
SJ="${SJ:-/path/to/h2h_sr/sr_1000000/SJ.out.tab}" # 1 M wgsim pairs -> STAR (see ../sqanti3_filter_h2h/run.sh)
GENOME="${GENOME:-/path/to/chr22.fa}"
CFG="${CFG:-$HERE/../../config/rules.default.toml}"
OUT="${OUT:-/tmp/chain_support}"; mkdir -p "$OUT"
sed 's/^min_full_reads.*/min_full_reads = 1000000/' "$CFG" > "$OUT/cfg_off.toml"   # the chain axis can never fire

declare -A GTF=([flair]="$WORK/flair.isoforms.gtf" [isoquant]="$MC/isoquant_out/iq/iq.transcript_models.gtf"
                [bambu]="$MC/bambu.novel.gtf" [espresso]="$MC/espresso_out/espresso_samples_N2_R0_updated.gtf"
                [talon]="$MC/talon_out/talon_filt_talon.gtf")
declare -A CLS=([flair]="$WORK/sqanti_out/flair_classification.txt")
for c in isoquant bambu espresso talon; do CLS[$c]="$MC/sq_$c/${c}_classification.txt"; done

for c in flair isoquant bambu espresso talon; do
  common=(--classification "${CLS[$c]}" --isoforms-gtf "${GTF[$c]}" --ref-gtf "$WORK/chr22_modified.gtf"
          --bam "$WORK/aln.bam" --reference "$GENOME")
  "$PIG" adjudicate --config "$CFG"           "${common[@]}" --sj-tab "$SJ" --out-prefix "$OUT/$c.sr1M.on"   > /dev/null 2>&1
  "$PIG" adjudicate --config "$OUT/cfg_off.toml" "${common[@]}" --sj-tab "$SJ" --out-prefix "$OUT/$c.sr1M.off" > /dev/null 2>&1
  "$PIG" adjudicate --config "$CFG"           "${common[@]}"               --out-prefix "$OUT/$c.nosr.on"   > /dev/null 2>&1
  "$PIG" adjudicate --config "$OUT/cfg_off.toml" "${common[@]}"               --out-prefix "$OUT/$c.nosr.off" > /dev/null 2>&1
  echo "adjudicated $c"
done
"${PY:-python3}" "$HERE/score_chain.py" --work "$WORK" --out "$OUT" \
  $(for c in flair isoquant bambu espresso talon; do echo "--gtf $c:${GTF[$c]} --cls $c:${CLS[$c]}"; done) \
  --emit-metrics "$HERE/../results/chain_support/metrics.json"
