#!/usr/bin/env bash
# BAM mapping axis: wire + validate the indel-near / soft-clip signals.
#
# These read-level fractions are extracted by the BAM reader at every junction but,
# before v0.0.4, never fed the verdict (only low-MAPQ + supplementary did). This
# protocol adjudicates the SAME chr22 SQANTI-SIM inputs twice with the current binary
# -- ON (shipped default config) vs OFF (indel/softclip thresholds raised > 1.0 so they
# can never fire) -- and scores both against truth to measure the discrimination and the
# before/after impact. Fully reproducible; public/simulated data only.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/sqanti_sim/work}"   # the shared sqanti_sim work dir (FLAIR + SQANTI3 + aln.bam)
CFG="${CFG:-$HERE/../../config/rules.default.toml}"
OUT="${OUT:-/tmp/bam_axis}"; mkdir -p "$OUT"

cd "$WORK"

common=(--classification sqanti_out/flair_classification.txt --isoforms-gtf flair.isoforms.gtf
        --ref-gtf chr22_modified.gtf --sj-tab truth.SJ.tab --bam aln.bam)

# ON: the shipped default config (indel/softclip gates at 0.5).
"$PIG" adjudicate --config "$CFG" "${common[@]}" --out-prefix "$OUT/on"

# OFF: disable just the two newly-wired triggers by raising them above 1.0.
sed 's/^max_indel_near_frac.*/max_indel_near_frac = 1.1/; s/^max_softclip_frac.*/max_softclip_frac = 1.1/' \
    "$CFG" > "$OUT/cfg_off.toml"
"$PIG" adjudicate --config "$OUT/cfg_off.toml" "${common[@]}" --out-prefix "$OUT/off"

python3 "$HERE/score_bam_axis.py" \
    --flair flair.isoforms.gtf --truth truth.truth.tsv \
    --off "$OUT/off.attribution.jsonl" --on "$OUT/on.attribution.jsonl" \
    --emit-metrics "$HERE/../results/bam_axis/metrics.json"
