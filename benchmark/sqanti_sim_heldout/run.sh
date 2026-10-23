#!/usr/bin/env bash
# Held-out chromosome: the chr22 protocol (../sqanti_sim + ../sqanti3_filter_h2h) repeated unchanged on
# another chromosome, with PanIsoGuard's thresholds left exactly as set on chr22. Two short-read
# settings, the two of the main text: none, and 1 M wgsim pairs -> STAR 2-pass.
# Usage: bash run.sh chr20 [/mnt/Data/pig_sqantisim_chr20]      then emit_heldout.py over the runs
set -euo pipefail
CHROM=$1; ROOT=${2:-/mnt/Data/pig_sqantisim_$CHROM}
SQSIM=$HOME/SQANTI-SIM; SQSIM_BIN=$HOME/miniconda3/envs/sqsim/bin; FLAIR_BIN=$HOME/miniconda3/envs/flair/bin
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin; STAR=$HOME/miniconda3/envs/star274/bin/STAR; PY=$SQSIM_BIN/python
PIG=$HOME/PanIsoGuard/build/panisoguard
GENOME=${GENOME:-/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa}
FULL_GTF=${FULL_GTF:-/mnt/Data/LR-RNAseq/reference/gencode.v49.annotation.gtf}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$ROOT/work; mkdir -p "$W"; cd "$W"
export PATH="$SQSIM_BIN:$PATH"
echo "[$(date +%T)] $CHROM: inputs"
[ -s $CHROM.gtf ] || awk -F'\t' -v c=$CHROM '$1==c' "$FULL_GTF" > $CHROM.gtf
[ -s $CHROM.fa.fai ] || { samtools faidx "$GENOME" $CHROM > $CHROM.fa && samtools faidx $CHROM.fa; }
echo "[$(date +%T)] SQANTI-SIM classif / design / sim (same parameters as chr22)"
[ -s ${CHROM}_index.tsv ] || python "$SQSIM/sqanti-sim.py" classif --gtf $CHROM.gtf -o $CHROM -d "$W" --cores 8 > classif.log 2>&1
[ -s ${CHROM}_modified.gtf ] || python "$SQSIM/sqanti-sim.py" design equal -i ${CHROM}_index.tsv --gtf $CHROM.gtf \
    -o $CHROM -d "$W" -nt 1500 --NIC 300 --NNC 300 --ISM 200 --read_count 100000 -s 1 -k 8 > design.log 2>&1
[ -s PBSIM3_simulated.fasta ] || python "$SQSIM/sqanti-sim.py" sim --gtf $CHROM.gtf --genome $CHROM.fa -i ${CHROM}_index.tsv \
    -d "$W" --pb --pbsim --pbsim_model "$SQSIM_BIN/../data/QSHMM-RSII.model" --long_count 100000 -k 8 > sim.log 2>&1
echo "[$(date +%T)] minimap2 / FLAIR collapse vs the reduced annotation"
[ -s aln.bam.bai ] || { minimap2 -ax splice:hq -t 8 $CHROM.fa PBSIM3_simulated.fasta 2> mm2.log | samtools sort -o aln.bam - && samtools index aln.bam; }
[ -s reads.bed12 ] || python "$HERE/../end2end/bam2bed12.py" aln.bam reads.bed12
[ -s flair.isoforms.gtf ] || PATH="$FLAIR_BIN:$SQSIM_BIN:$PATH" flair collapse -g $CHROM.fa -q reads.bed12 \
    -r PBSIM3_simulated.fasta --gtf ${CHROM}_modified.gtf -o flair --threads 8 > collapse.log 2>&1
echo "[$(date +%T)] SQANTI3 qc vs the reduced annotation; truth"
[ -s sqanti_out/flair_classification.txt ] || PATH="$SQANTI3_BIN:$PATH" sqanti3_qc.py --isoforms flair.isoforms.gtf \
    --refGTF ${CHROM}_modified.gtf --refFasta $CHROM.fa -o flair -d sqanti_out --report skip --cpus 8 > sqanti.log 2>&1
[ -s truth.truth.tsv ] || python "$HERE/../sqanti_sim/prep_truth_sj.py" $CHROM.gtf ${CHROM}_modified.gtf PBSIM3_simulated.read_to_isoform.tsv truth

# --- head-to-head (as ../sqanti3_filter_h2h/run.sh, settings none + sr_1000000) ----------------
H=$ROOT/h2h; mkdir -p "$H"; cd "$H"
export PATH="$SQANTI3_BIN:$PATH"
echo "[$(date +%T)] short reads: wgsim 1 M pairs -> STAR 2-pass"
if [ ! -s tx.weighted.fa ]; then
  gffread -w tx.fa -g "$W/$CHROM.fa" "$W/$CHROM.gtf"
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
[ -s star_idx/SA ] || { mkdir -p star_idx; "$STAR" --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles "$W/$CHROM.fa" \
    --sjdbGTFfile "$W/${CHROM}_modified.gtf" --sjdbOverhang 99 --genomeSAindexNbases 11 --runThreadN 16 > star_idx.log 2>&1; }
N=1000000; d=sr_$N; mkdir -p $d
[ -s $d/SJ.out.tab ] || {
  wgsim -N "$N" -1 100 -2 100 -d 250 -s 30 -e 0.002 -r 0 -R 0 -S 11 tx.weighted.fa $d/r1.fq $d/r2.fq > /dev/null 2>&1
  "$STAR" --genomeDir star_idx --readFilesIn $d/r1.fq $d/r2.fq --twopassMode Basic --outSAMtype None --runThreadN 16 --outFileNamePrefix $d/ > /dev/null
  rm -f $d/r1.fq $d/r2.fq; }
cat > strict.json <<'EOF'
{
    "full-splice_match": [ { "perc_A_downstream_TTS":[0,59] } ],
    "rest": [ { "perc_A_downstream_TTS":[0,59], "RTS_stage":"FALSE", "min_cov":3 } ]
}
EOF
common=(--isoforms-gtf "$W/flair.isoforms.gtf" --ref-gtf "$W/${CHROM}_modified.gtf")
for s in none sr_$N; do
  echo "[$(date +%T)] setting $s: SQANTI3 filters + PanIsoGuard"
  mkdir -p $s
  case $s in none) sj="";; *) sj=$PWD/$s/SJ.out.tab;; esac
  cls=$W/sqanti_out/flair_classification.txt
  if [ -n "$sj" ]; then
    [ -s $s/sq/flair_classification.txt ] || sqanti3_qc.py --isoforms "$W/flair.isoforms.gtf" --refGTF "$W/${CHROM}_modified.gtf" \
        --refFasta "$W/$CHROM.fa" -c "$sj" -o flair -d $PWD/$s/sq --report skip --cpus 8 > $s/sq.log 2>&1
    cls=$PWD/$s/sq/flair_classification.txt
    [ -s $s/strict/strict_RulesFilter_classification.txt ] || sqanti3_filter.py rules --sqanti_class "$cls" -j strict.json -o strict -d $PWD/$s/strict --skip_report > /dev/null 2>&1
  fi
  [ -s $s/rules/rules_RulesFilter_classification.txt ] || sqanti3_filter.py rules --sqanti_class "$cls" -o rules -d $PWD/$s/rules --skip_report > /dev/null 2>&1
  [ -s $s/ml/ml_ML_classification.txt ] || sqanti3_filter.py ml --sqanti_class "$cls" -o ml -d $PWD/$s/ml --skip_report > /dev/null 2>&1
  pig=(adjudicate --classification "$W/sqanti_out/flair_classification.txt" "${common[@]}")
  [ -n "$sj" ] && pig+=(--sj-tab "$sj")
  "$PIG" "${pig[@]}" --out-prefix $s/pig > /dev/null
  "$PIG" "${pig[@]}" --bam "$W/aln.bam" --reference "$W/$CHROM.fa" --out-prefix $s/pigbam > /dev/null
done
echo "[$(date +%T)] $CHROM done: $ROOT"
