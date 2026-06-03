#!/usr/bin/env bash
# Multi-caller consensus benchmark: five long-read isoform callers on one shared
# alignment, integrated by PanIsoGuard `combine`, scored against SQANTI-SIM truth.
#
# Demonstrates the caller-agnostic value proposition: a single caller's novel calls
# carry that caller's (often poor) false-discovery rate; multi-caller agreement,
# surfaced by `combine` and consumed by the `adjudicate` consensus axis, separates
# genuine novel isoforms from per-caller artifacts. score_multicaller.py also emits a
# precision/recall curve over the consensus threshold (the F1-optimal threshold scales
# with the number of callers).
#
# Inputs are gitignored and live off-repo. Paths below are placeholders to fill in; the
# committed result is benchmark/results/multicaller/metrics.json.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# --- placeholders (fill in) -------------------------------------------------
WORK=/path/to/sqanti_sim/work        # aln.bam, chr22_modified.gtf, truth.truth.tsv, flair.isoforms.gtf, sqanti_out/
GENOME=/path/to/chr22.fa             # + .fai
ISOQUANT=isoquant                    # conda env: isoquant
BAMBU_RS=Rscript                     # conda env: bambu (bioconductor-bambu + r-biocmanager)
ESPRESSO_BIN=$HOME/miniconda3/envs/espresso/bin
SAMTOOLS=samtools
# TALON 6.0 on py3.7 needs a pyranges/importlib.metadata shim (see notes at bottom).
TALON_BIN=$HOME/miniconda3/envs/talon/bin
# ----------------------------------------------------------------------------

REF="$WORK/chr22_modified.gtf"; ALN="$WORK/aln.bam"
FLAIR_GTF="$WORK/flair.isoforms.gtf"          # FLAIR collapse from the sqanti_sim protocol
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT"     # data/ is gitignored

# 1) IsoQuant (PacBio-CCS-like; reads were minimap2 splice:hq)
"$ISOQUANT" --reference "$GENOME" --genedb "$REF" --complete_genedb \
  --bam "$ALN" --data_type pacbio_ccs -o "$OUT/isoquant_out" --threads 8 --prefix iq
IQ_GTF="$OUT/isoquant_out/iq/iq.transcript_models.gtf"

# 2) Bambu; keep only NOVEL (BambuTx) transcripts (it also echoes the full reference).
"$BAMBU_RS" "$HERE/run_bambu_ndr.R" "$ALN" "$REF" "$GENOME" "$OUT/bambu_out" 4
awk -F'\t' '$9 ~ /transcript_id "BambuTx/' "$OUT/bambu_out/extended_annotations.gtf" > "$OUT/bambu.novel.gtf"

# 3) ESPRESSO (S correct/detect, C per-read correction, Q quantify+GTF); novel rows have
#    source "novel_isoform" -- combine annotates novelty, so the full GTF is fine.
printf '%s\tsim1\n' "$ALN" > "$OUT/espresso_samples.tsv"
perl "$ESPRESSO_BIN/ESPRESSO_S.pl" -L "$OUT/espresso_samples.tsv" -F "$GENOME" -A "$REF" -O "$OUT/espresso_out" -T 8
perl "$ESPRESSO_BIN/ESPRESSO_C.pl" -I "$OUT/espresso_out" -F "$GENOME" -X 0 -T 8
perl "$ESPRESSO_BIN/ESPRESSO_Q.pl" -A "$REF" -L "$OUT/espresso_out/espresso_samples.tsv.updated" -O "$OUT/espresso_out" -T 8
ESP_GTF="$OUT/espresso_out/espresso_samples_N2_R0_updated.gtf"

# 4) TALON: needs MD-tagged SAM. Use the STANDARD filtered whitelist (minCount 5) so the
#    novel set is comparable -- talon --observed alone is one-transcript-per-read noise.
"$SAMTOOLS" calmd -@ 8 "$ALN" "$GENOME" 2>/dev/null > "$OUT/aln.md.sam"
"$TALON_BIN/talon_label_reads" --f "$OUT/aln.md.sam" --g "$GENOME" --t 8 --o "$OUT/talon/sim1" --deleteTmp
"$TALON_BIN/talon_initialize_database" --f "$REF" --g chr22 --a gencode_chr22 --o "$OUT/talon/talon"
printf 'sim1,sim1,SequelII,%s\n' "$OUT/talon/sim1_labeled.sam" > "$OUT/talon/config.csv"
"$TALON_BIN/talon" --f "$OUT/talon/config.csv" --db "$OUT/talon/talon.db" --build chr22 --threads 8 --o "$OUT/talon/run"
"$TALON_BIN/talon_filter_transcripts" --db "$OUT/talon/talon.db" -a gencode_chr22 \
  --datasets sim1 --maxFracA 0.5 --minCount 5 --minDatasets 1 --o "$OUT/talon/whitelist.csv"
"$TALON_BIN/talon_create_GTF" --db "$OUT/talon/talon.db" --build chr22 -a gencode_chr22 \
  --whitelist "$OUT/talon/whitelist.csv" --o "$OUT/talon/talon_filt"
TALON_GTF="$OUT/talon/talon_filt_talon.gtf"

# 5) Integrate all five caller-agnostically by intron-chain fingerprint.
"$PIG" combine --gtf flair:"$FLAIR_GTF" --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" \
  --gtf espresso:"$ESP_GTF" --gtf talon:"$TALON_GTF" --ref-gtf "$REF" --out "$OUT/matrix.tsv"

# 6) Engine demo: adjudicate FLAIR long-read-only, with vs without the consensus matrix.
FLAIR_CLS="$WORK/sqanti_out/flair_classification.txt"
"$PIG" adjudicate --classification "$FLAIR_CLS" --isoforms-gtf "$FLAIR_GTF" --ref-gtf "$REF" --out-prefix "$OUT/flair_base"
"$PIG" adjudicate --classification "$FLAIR_CLS" --isoforms-gtf "$FLAIR_GTF" --ref-gtf "$REF" \
  --caller-support "$OUT/matrix.tsv" --out-prefix "$OUT/flair_cons"

# 7) Score: per-caller, consensus stratification, PR curve, engine demo -> envelope.
python3 "$HERE/score_multicaller.py" --matrix "$OUT/matrix.tsv" --truth "$WORK/truth.truth.tsv" \
  --gtf flair:"$FLAIR_GTF" --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" \
  --gtf espresso:"$ESP_GTF" --gtf talon:"$TALON_GTF" \
  --min-callers 2 --tool-version "$("$PIG" --version | head -1 | awk '{print $2}')" \
  --ruleset-version builtin-0.0.1 --emit-metrics "$HERE/../results/multicaller/metrics.json"

# TALON py3.7 shim (one-time, if `import pyranges` fails on importlib.metadata):
#   pip install 'importlib-metadata'; add sitecustomize.py to site-packages:
#     import sys, importlib, importlib_metadata
#     sys.modules['importlib.metadata'] = importlib_metadata
#     importlib.metadata = importlib_metadata
