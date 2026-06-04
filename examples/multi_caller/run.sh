#!/usr/bin/env bash
# Self-contained multi-caller quickstart: combine 3 callers -> adjudicate with the
# consensus axis -> (optional) PDF report. Tiny committed fixtures, runs in <1s offline.
# Shows the headline: a novel chain recovered by >= 2 callers is promoted to
# MEDIUM_CONF_NOVEL, a single-caller novel is held AMBIGUOUS.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd); cd "$SCRIPT_DIR"

if   [[ -n "${PANISOGUARD:-}" ]]; then PIG="$PANISOGUARD"
elif command -v panisoguard >/dev/null 2>&1; then PIG=panisoguard
elif [[ -x ../../build/panisoguard ]]; then PIG=../../build/panisoguard
else echo "error: panisoguard not found (set PANISOGUARD=… or build ../../build/panisoguard)" >&2; exit 1; fi

OUT=${1:-output}; mkdir -p "$OUT"

# 1) integrate the 3 callers by intron-chain fingerprint -> caller-support matrix
"$PIG" combine \
  --gtf flair:data/flair.gtf \
  --gtf isoquant:data/isoquant.gtf \
  --gtf bambu:data/bambu.gtf \
  --ref-gtf data/reference.gtf \
  --out "$OUT/caller_support.tsv"

# 2) adjudicate FLAIR's isoforms, feeding cross-caller agreement (no short-read here,
#    so the consensus axis decides). >=2 callers -> MEDIUM_CONF_NOVEL; 1 caller -> AMBIGUOUS.
"$PIG" adjudicate \
  --classification data/flair.classification.tsv \
  --isoforms-gtf data/flair.gtf \
  --ref-gtf data/reference.gtf \
  --caller-support "$OUT/caller_support.tsv" \
  --out-prefix "$OUT/sample"

echo
echo "=== verdict per isoform (isoA: 3 callers, isoB: 2, isoC: 1) ==="
awk -F'\t' 'NR==1 || $1 ~ /^iso/ {printf "  %-6s %-22s %s\n", $1, $7, $5}' "$OUT/sample.adjudicated.tsv"

# 3) (optional) a SQANTI3-style PDF report, if panisoguard-report is installed
if command -v panisoguard-report >/dev/null 2>&1; then
  panisoguard-report --prefix "$OUT/sample" && echo "wrote $OUT/sample.report.pdf"
else
  echo "(install ./python for a PDF report: pip install ../../python && panisoguard-report --prefix $OUT/sample)"
fi

# compare against the committed expected verdicts
if [[ -f expected/sample.adjudicated.tsv ]]; then
  diff -u expected/sample.adjudicated.tsv "$OUT/sample.adjudicated.tsv" && echo "output matches expected"
fi
