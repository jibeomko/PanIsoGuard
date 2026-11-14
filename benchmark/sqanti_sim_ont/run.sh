#!/usr/bin/env bash
# SQANTI-SIM chr22 with ONT cDNA reads. The design of the PacBio run in ../sqanti_sim (the reduced
# annotation, the hidden transcripts and their requested counts) is reused unchanged, so the truth set is
# the same and only the chemistry differs: NanoSim (SQANTI-SIM's pre-trained human WTC11 cDNA guppy
# model) instead of PBSIM3, minimap2 -ax splice instead of splice:hq. Then FLAIR, SQANTI3, 1 M wgsim
# pairs -> STAR 2-pass and the SQANTI3 filters as in ../sqanti3_filter_h2h, and PanIsoGuard without the
# BAM, with the BAM (default hifi preset) and with --platform ont. Thresholds untouched.
# Usage: bash run.sh [ROOT]   (env: PB = the ../sqanti_sim work root, THREADS)
set -euo pipefail
ROOT=${1:-/mnt/Data/PanIsoGuard_bench/pig_sqantisim_ont}
PB=${PB:-/mnt/Data/PanIsoGuard_bench/pig_sqantisim}     # the PacBio run: work/{chr22.gtf,chr22.fa,chr22_index.tsv,chr22_modified.gtf}
SQSIM=$HOME/SQANTI-SIM; SQSIM_BIN=$HOME/miniconda3/envs/sqsim/bin; FLAIR_BIN=$HOME/miniconda3/envs/flair/bin
SQANTI3_BIN=$HOME/miniconda3/envs/sqanti3/bin; STAR=$HOME/miniconda3/envs/star274/bin/STAR; PY=$SQSIM_BIN/python
PIG=$HOME/PanIsoGuard/build/panisoguard
HERE=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-16}
W=$ROOT/work; mkdir -p "$W"; cd "$W"
export PATH="$SQSIM_BIN:$PATH"

echo "[$(date +%T)] inputs: the chr22 design of $PB"
for f in chr22.gtf chr22.fa chr22.fa.fai chr22_modified.gtf chr22_index.tsv; do   # sim rewrites sim_counts in the index: work on copies
  [ -s $f ] || { [ -s "$PB/work/$f" ] && cp "$PB/work/$f" . || cp "$PB/$f" .; }
done

echo "[$(date +%T)] NanoSim: SQANTI-SIM's WTC11 cDNA model, 100,000 reads, seed 1 (NanoSim runs single-threaded)"
# nanosim_compat/sitecustomize.py lets the scikit-learn-1.2 model pickles load under a newer scikit-learn
[ -s ONT_simulated.fastq ] || PYTHONPATH="$HERE/nanosim_compat${PYTHONPATH:+:$PYTHONPATH}" python "$SQSIM/sqanti-sim.py" sim --gtf chr22.gtf --genome chr22.fa -i chr22_index.tsv \
    -d "$W" --ont --read_type cDNA --long_count 100000 -s 1 -k 8 > sim.log 2>&1

echo "[$(date +%T)] minimap2 -ax splice / FLAIR collapse vs the reduced annotation"
[ -s aln.bam.bai ] || { minimap2 -ax splice -t "$THREADS" chr22.fa ONT_simulated.fastq 2> mm2.log | samtools sort -@4 -o aln.bam - && samtools index aln.bam; }
[ -s reads.bed12 ] || python "$HERE/../end2end/bam2bed12.py" aln.bam reads.bed12
[ -s flair.isoforms.gtf ] || PATH="$FLAIR_BIN:$SQSIM_BIN:$PATH" flair collapse -g chr22.fa -q reads.bed12 \
    -r ONT_simulated.fastq --gtf chr22_modified.gtf -o flair --threads "$THREADS" > collapse.log 2>&1

echo "[$(date +%T)] SQANTI3 qc vs the reduced annotation; truth"
[ -s sqanti_out/flair_classification.txt ] || PATH="$SQANTI3_BIN:$PATH" sqanti3_qc.py --isoforms flair.isoforms.gtf \
    --refGTF chr22_modified.gtf --refFasta chr22.fa -o flair -d sqanti_out --report skip --cpus 8 > sqanti.log 2>&1
[ -s truth.truth.tsv ] || python "$HERE/../sqanti_sim/prep_truth_sj.py" chr22.gtf chr22_modified.gtf ONT_simulated.read_to_isoform.tsv truth

H=$ROOT/h2h; mkdir -p "$H"; cd "$H"
export PATH="$SQANTI3_BIN:$PATH"
echo "[$(date +%T)] short reads: wgsim 1 M pairs -> STAR 2-pass"
if [ ! -s tx.weighted.fa ]; then
  gffread -w tx.fa -g "$W/chr22.fa" "$W/chr22.gtf"
  # one FASTA copy per simulated long read, so wgsim samples pairs in proportion to expression x length
  "$PY" - "$W/ONT_simulated.read_to_isoform.tsv" tx.fa > tx.weighted.fa <<'EOF'
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
[ -s star_idx/SA ] || { mkdir -p star_idx; "$STAR" --runMode genomeGenerate --genomeDir star_idx --genomeFastaFiles "$W/chr22.fa" \
    --sjdbGTFfile "$W/chr22_modified.gtf" --sjdbOverhang 99 --genomeSAindexNbases 11 --runThreadN "$THREADS" > star_idx.log 2>&1; }
N=1000000; d=sr_$N; mkdir -p $d
[ -s $d/SJ.out.tab ] || {
  wgsim -N "$N" -1 100 -2 100 -d 250 -s 30 -e 0.002 -r 0 -R 0 -S 11 tx.weighted.fa $d/r1.fq $d/r2.fq > /dev/null 2>&1
  "$STAR" --genomeDir star_idx --readFilesIn $d/r1.fq $d/r2.fq --twopassMode Basic --outSAMtype None --runThreadN "$THREADS" --outFileNamePrefix $d/ > /dev/null
  rm -f $d/r1.fq $d/r2.fq; }
cat > strict.json <<'EOF'
{
    "full-splice_match": [ { "perc_A_downstream_TTS":[0,59] } ],
    "rest": [ { "perc_A_downstream_TTS":[0,59], "RTS_stage":"FALSE", "min_cov":3 } ]
}
EOF
common=(--isoforms-gtf "$W/flair.isoforms.gtf" --ref-gtf "$W/chr22_modified.gtf")
for s in none sr_$N; do
  echo "[$(date +%T)] setting $s: SQANTI3 filters + PanIsoGuard (no BAM / BAM, hifi preset / BAM, ont preset)"
  mkdir -p $s
  case $s in none) sj="";; *) sj=$PWD/$s/SJ.out.tab;; esac
  cls=$W/sqanti_out/flair_classification.txt
  if [ -n "$sj" ]; then
    [ -s $s/sq/flair_classification.txt ] || sqanti3_qc.py --isoforms "$W/flair.isoforms.gtf" --refGTF "$W/chr22_modified.gtf" \
        --refFasta "$W/chr22.fa" -c "$sj" -o flair -d $PWD/$s/sq --report skip --cpus 8 > $s/sq.log 2>&1
    cls=$PWD/$s/sq/flair_classification.txt
    [ -s $s/strict/strict_RulesFilter_classification.txt ] || sqanti3_filter.py rules --sqanti_class "$cls" -j strict.json -o strict -d $PWD/$s/strict --skip_report > /dev/null 2>&1
  fi
  [ -s $s/rules/rules_RulesFilter_classification.txt ] || sqanti3_filter.py rules --sqanti_class "$cls" -o rules -d $PWD/$s/rules --skip_report > /dev/null 2>&1
  [ -s $s/ml/ml_ML_classification.txt ] || sqanti3_filter.py ml --sqanti_class "$cls" -o ml -d $PWD/$s/ml --skip_report > /dev/null 2>&1
  pig=(adjudicate --classification "$W/sqanti_out/flair_classification.txt" "${common[@]}")
  [ -n "$sj" ] && pig+=(--sj-tab "$sj")
  "$PIG" "${pig[@]}" --out-prefix $s/pig > /dev/null
  "$PIG" "${pig[@]}" --bam "$W/aln.bam" --reference "$W/chr22.fa" --out-prefix $s/pigbam > /dev/null
  "$PIG" "${pig[@]}" --bam "$W/aln.bam" --reference "$W/chr22.fa" --platform ont --out-prefix $s/pigbam_ont > /dev/null
done
"$PY" "$HERE/emit_ont.py" "$ROOT" --emit-metrics "$HERE/../results/sqanti_sim_ont/metrics.json"
echo "[$(date +%T)] done: $ROOT"
