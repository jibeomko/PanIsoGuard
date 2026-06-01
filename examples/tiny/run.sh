#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

if [[ -n "${PANISOGUARD:-}" ]]; then
  PIG="$PANISOGUARD"
elif command -v panisoguard >/dev/null 2>&1; then
  PIG=panisoguard
elif [[ -x ../../build/panisoguard ]]; then
  PIG=../../build/panisoguard
else
  echo "error: panisoguard not found. Set PANISOGUARD=/path/to/panisoguard or build ../../build/panisoguard" >&2
  exit 1
fi

OUT_DIR=${1:-output}
mkdir -p "$OUT_DIR"

"$PIG" adjudicate \
  --classification data/classification.tsv \
  --isoforms-gtf data/caller.gtf \
  --ref-gtf data/reference.gtf \
  --sj-tab data/short_read.SJ.tab \
  --out-prefix "$OUT_DIR/sample"

if [[ -f expected/sample.adjudicated.tsv && -f expected/sample.attribution.jsonl && -f expected/sample.provenance.log ]]; then
  diff -u expected/sample.adjudicated.tsv "$OUT_DIR/sample.adjudicated.tsv"
  diff -u expected/sample.attribution.jsonl "$OUT_DIR/sample.attribution.jsonl"
  diff -u expected/sample.provenance.log "$OUT_DIR/sample.provenance.log"
  echo "example output matches expected files"
else
  echo "expected files not found; generated output under $OUT_DIR/"
fi
