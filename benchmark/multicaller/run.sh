#!/usr/bin/env bash
# Multi-caller consensus benchmark: three long-read isoform callers on one shared
# alignment, integrated by PanIsoGuard `combine`, scored against SQANTI-SIM truth.
#
# Demonstrates the caller-agnostic value proposition: a single caller's novel calls
# carry that caller's (often poor) false-discovery rate; multi-caller agreement,
# surfaced by `combine` and consumed by the `adjudicate` consensus axis, separates
# genuine novel isoforms from per-caller artifacts.
#
# Inputs are gitignored and live off-repo (see the data-locations note). Paths below
# are the placeholders to fill in; the committed result is benchmark/results/
# multicaller/metrics.json.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# --- placeholders (fill in) -------------------------------------------------
WORK=/path/to/sqanti_sim/work            # has aln.bam, chr22_modified.gtf, truth.truth.tsv, flair.isoforms.gtf
GENOME=/path/to/chr22.fa                 # + .fai
ISOQUANT_BIN=isoquant                    # conda env: isoquant
BAMBU_RSCRIPT=Rscript                    # conda env with bioconductor-bambu (+ r-biocmanager)
# ----------------------------------------------------------------------------

REF="$WORK/chr22_modified.gtf"
ALN="$WORK/aln.bam"
FLAIR_GTF="$WORK/flair.isoforms.gtf"     # FLAIR collapse already run by the sqanti_sim protocol
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT"   # data/ is gitignored (heavy intermediates)

# 1) IsoQuant on the shared alignment (PacBio-CCS-like; reads were minimap2 splice:hq).
"$ISOQUANT_BIN" --reference "$GENOME" --genedb "$REF" --complete_genedb \
  --bam "$ALN" --data_type pacbio_ccs -o "$OUT/isoquant_out" --threads 8 --prefix iq
IQ_GTF="$OUT/isoquant_out/iq/iq.transcript_models.gtf"

# 2) Bambu on the shared alignment; keep only its NOVEL (BambuTx) transcripts so the
#    set is comparable to the other callers' discoveries (bambu also echoes the full
#    reference annotation). NDR fixed at 0.5 to avoid the online NDR-recommendation step.
"$BAMBU_RSCRIPT" "$HERE/run_bambu_ndr.R" "$ALN" "$REF" "$GENOME" "$OUT/bambu_out" 4
awk -F'\t' '$9 ~ /transcript_id "BambuTx/' "$OUT/bambu_out/extended_annotations.gtf" > "$OUT/bambu.novel.gtf"

# 3) Integrate the three callers caller-agnostically by intron-chain fingerprint.
"$PIG" combine --gtf flair:"$FLAIR_GTF" --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" \
  --ref-gtf "$REF" --out "$OUT/matrix.tsv"

# 4) Engine demo: adjudicate one caller (FLAIR) long-read-only, with vs without the
#    cross-caller consensus matrix. Baseline holds every novel AMBIGUOUS; consensus
#    promotes the cross-caller-agreed novels to MEDIUM_CONF_NOVEL.
FLAIR_CLS="$WORK/sqanti_out/flair_classification.txt"
"$PIG" adjudicate --classification "$FLAIR_CLS" --isoforms-gtf "$FLAIR_GTF" --ref-gtf "$REF" \
  --out-prefix "$OUT/flair_base"
"$PIG" adjudicate --classification "$FLAIR_CLS" --isoforms-gtf "$FLAIR_GTF" --ref-gtf "$REF" \
  --caller-support "$OUT/matrix.tsv" --out-prefix "$OUT/flair_cons"

# 5) Score: per-caller vs consensus stratification vs truth, and emit the envelope.
python3 "$HERE/score_multicaller.py" --matrix "$OUT/matrix.tsv" \
  --truth "$WORK/truth.truth.tsv" \
  --gtf flair:"$FLAIR_GTF" --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" \
  --min-callers 2 --tool-version "$("$PIG" --version | head -1 | awk '{print $2}')" \
  --ruleset-version builtin-0.0.1 \
  --emit-metrics "$HERE/../results/multicaller/metrics.json"
