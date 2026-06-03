#!/usr/bin/env bash
# Whole-genome multi-caller consensus on REAL public data (GM12878 ONT dRNA).
# Three callers on one shared genome alignment, integrated by `combine`, validated by
# splice-motif canonicality (no ground truth on real data). Honest, small-N result --
# see README.md. All inputs public (ENCODE ENCSR368UNC; GRCh38/GENCODE v49); no private data.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"

# --- placeholders (fill in) -------------------------------------------------
BAM=/path/to/gm12878.sorted.bam        # ENCODE GM12878 ONT dRNA, minimap2 to GRCh38
GENOME=/path/to/GRCh38.primary_assembly.genome.fa   # + .fai
REFGTF=/path/to/gencode.v49.annotation.gtf
ISOQUANT=isoquant; BAMBU_RS=Rscript; ESPRESSO_BIN=$HOME/miniconda3/envs/espresso/bin
SAMTOOLS=samtools
# ----------------------------------------------------------------------------
OUT="${OUT:-$HERE/data}"; mkdir -p "$OUT"   # data/ is gitignored

# 0) restrict to main chromosomes so ESPRESSO does not choke on alt/random contigs
MAIN=$(for i in $(seq 1 22) X Y M; do echo -n "chr$i "; done)
"$SAMTOOLS" view -@ 8 -b "$BAM" $MAIN > "$OUT/gm.main.bam"; "$SAMTOOLS" index "$OUT/gm.main.bam"

# 1) IsoQuant (ONT dRNA)
"$ISOQUANT" --reference "$GENOME" --genedb "$REFGTF" --complete_genedb \
  --bam "$OUT/gm.main.bam" --data_type nanopore -o "$OUT/isoquant_out" --threads 8 --prefix gm
IQ_GTF="$OUT/isoquant_out/gm/gm.transcript_models.gtf"

# 2) Bambu; keep NOVEL (BambuTx) transcripts only
"$BAMBU_RS" "$HERE/../multicaller/run_bambu_ndr.R" "$OUT/gm.main.bam" "$REFGTF" "$GENOME" "$OUT/bambu_out" 6
awk -F'\t' '$9 ~ /transcript_id "BambuTx/' "$OUT/bambu_out/extended_annotations.gtf" > "$OUT/bambu.novel.gtf"

# 3) ESPRESSO (S/C/Q)
printf '%s\tgm1\n' "$OUT/gm.main.bam" > "$OUT/espresso_samples.tsv"
perl "$ESPRESSO_BIN/ESPRESSO_S.pl" -L "$OUT/espresso_samples.tsv" -F "$GENOME" -A "$REFGTF" -O "$OUT/espresso_out" -T 8
perl "$ESPRESSO_BIN/ESPRESSO_C.pl" -I "$OUT/espresso_out" -F "$GENOME" -X 0 -T 8
perl "$ESPRESSO_BIN/ESPRESSO_Q.pl" -A "$REFGTF" -L "$OUT/espresso_out/espresso_samples.tsv.updated" -O "$OUT/espresso_out" -T 8
ESP_GTF=$(ls "$OUT"/espresso_out/*_updated.gtf | head -1)

# 4) integrate + score (canonical-motif stratification)
"$PIG" combine --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" --gtf espresso:"$ESP_GTF" \
  --ref-gtf "$REFGTF" --out "$OUT/matrix.tsv"
python3 "$HERE/score_wholegenome.py" --matrix "$OUT/matrix.tsv" --genome "$GENOME" \
  --gtf isoquant:"$IQ_GTF" --gtf bambu:"$OUT/bambu.novel.gtf" --gtf espresso:"$ESP_GTF" \
  --sample gm12878 --tool-version "$("$PIG" --version | head -1 | awk '{print $2}')" \
  --ruleset-version builtin-0.0.1 --emit-metrics "$HERE/../results/wholegenome_multicaller/metrics.json"
