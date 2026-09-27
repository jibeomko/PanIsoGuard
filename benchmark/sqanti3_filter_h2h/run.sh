#!/usr/bin/env bash
# Head-to-head: PanIsoGuard vs the SQANTI3 rules / ML filter on the SQANTI-SIM truth set.
# Reuses the ../sqanti_sim work dir (run ../sqanti_sim/run.sh first). Two short-read settings:
#   oracle    the truth-derived SJ.tab from ../sqanti_sim (every simulated junction, 50 reads)
#   sr_<N>    realistic: N wgsim read pairs from the simulated transcripts (weighted by their
#             long-read counts) -> STAR 2-pass vs the reduced annotation -> SJ.out.tab
# Edit the paths below, then: bash run.sh
set -euo pipefail

# --- edit these ---------------------------------------------------------------
SQSIM_WORK=/path/to/sqanti_sim/work    # has flair.isoforms.gtf, chr22_modified.gtf, truth.*, aln.bam
FULL_GTF=/path/to/chr22.gtf            # unreduced chr22 annotation (transcript sequences)
GENOME=/path/to/chr22.fa
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin
STAR=$HOME/miniconda3/envs/star274/bin/STAR
PY=$HOME/miniconda3/envs/sqsim/bin/python   # needs scikit-learn
PIG=$HOME/PanIsoGuard/build/panisoguard
DEPTHS="250000 1000000 4000000"        # read pairs per realistic run
OUT=${OUT:-h2h_work}
# ------------------------------------------------------------------------------
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT" && cd "$OUT"
W=$SQSIM_WORK
export PATH=$SQANTI3_BIN:$PATH

# --- realistic short reads ----------------------------------------------------
if [ ! -s tx.weighted.fa ]; then
  gffread -w tx.fa -g "$GENOME" "$FULL_GTF"
  # one FASTA copy per simulated long-read molecule: wgsim samples pairs in proportion to
  # sequence length, so copies x length = expression x length, as in real RNA-seq.
  "$PY" - "$W/PBSIM3_simulated.read_to_isoform.tsv" tx.fa > tx.weighted.fa <<'EOF'
import sys
from collections import Counter
n = Counter(l.split('\t')[1].strip() for l in open(sys.argv[1]))
seq, name = {}, None
for l in open(sys.argv[2]):
    if l.startswith('>'):
        name = l[1:].split()[0]; seq[name] = []
    else:
        seq[name].append(l.strip())
for t, k in n.items():
    s = ''.join(seq[t])
    for i in range(k):
        print(f">{t}_{i}\n{s}")
EOF
fi
if [ ! -s star_idx/SA ]; then
  mkdir -p star_idx
  "$STAR" --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles "$GENOME" \
      --sjdbGTFfile "$W/chr22_modified.gtf" --sjdbOverhang 99 --genomeSAindexNbases 11 --runThreadN 16
fi
for N in $DEPTHS; do
  d=sr_$N; mkdir -p $d
  [ -s $d/SJ.out.tab ] || {
    wgsim -N "$N" -1 100 -2 100 -d 250 -s 30 -e 0.002 -r 0 -R 0 -S 11 tx.weighted.fa $d/r1.fq $d/r2.fq > /dev/null 2>&1
    "$STAR" --genomeDir star_idx --readFilesIn $d/r1.fq $d/r2.fq --twopassMode Basic \
        --outSAMtype None --runThreadN 16 --outFileNamePrefix $d/ > /dev/null
    rm -f $d/r1.fq $d/r2.fq
  }
done

# --- per short-read setting: SQANTI3 QC with coverage, filters, PanIsoGuard -----
cat > strict.json <<'EOF'
{
    "full-splice_match": [ { "perc_A_downstream_TTS":[0,59] } ],
    "rest": [ { "perc_A_downstream_TTS":[0,59], "RTS_stage":"FALSE", "min_cov":3 } ]
}
EOF
common=(--isoforms-gtf "$W/flair.isoforms.gtf" --ref-gtf "$W/chr22_modified.gtf")
settings="none oracle"; for N in $DEPTHS; do settings="$settings sr_$N"; done
for s in $settings; do
  mkdir -p $s
  case $s in none) sj="";; oracle) sj=$W/truth.SJ.tab;; *) sj=$PWD/$s/SJ.out.tab;; esac
  cls=$W/sqanti_out/flair_classification.txt
  if [ -n "$sj" ]; then
    [ -s $s/sq/flair_classification.txt ] || sqanti3_qc.py --isoforms "$W/flair.isoforms.gtf" \
        --refGTF "$W/chr22_modified.gtf" --refFasta "$GENOME" -c "$sj" -o flair -d $PWD/$s/sq \
        --report skip --cpus 8 > $s/sq.log 2>&1
    cls=$PWD/$s/sq/flair_classification.txt
    sqanti3_filter.py rules --sqanti_class "$cls" -j strict.json -o strict -d $PWD/$s/strict --skip_report > /dev/null 2>&1
  fi
  sqanti3_filter.py rules --sqanti_class "$cls" -o rules -d $PWD/$s/rules --skip_report > /dev/null 2>&1
  sqanti3_filter.py ml --sqanti_class "$cls" -o ml -d $PWD/$s/ml --skip_report > /dev/null 2>&1
  pig=(adjudicate --classification "$W/sqanti_out/flair_classification.txt" "${common[@]}")
  [ -n "$sj" ] && pig+=(--sj-tab "$sj")
  "$PIG" "${pig[@]}" --out-prefix $s/pig > /dev/null
  "$PIG" "${pig[@]}" --bam "$W/aln.bam" --reference "$GENOME" --out-prefix $s/pigbam > /dev/null

  m=("SQ3rules=rules:$s/rules/rules_RulesFilter_classification.txt"
     "SQ3ml=ml:$s/ml/ml_ML_classification.txt")
  [ -n "$sj" ] && m+=("SQ3strict=rules:$s/strict/strict_RulesFilter_classification.txt"
                      "allSJ3=mincov:$cls")
  m+=("PIG=pig:$s/pig.adjudicated.tsv" "PIG+bam=pig:$s/pigbam.adjudicated.tsv")
  echo "##### short reads: $s"
  "$PY" "$HERE/score_h2h.py" "$W/flair.isoforms.gtf" "$W/truth.truth.tsv" \
      "$W/sqanti_out/flair_classification.txt" "${m[@]}"
done | tee h2h.txt
