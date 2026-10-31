#!/usr/bin/env bash
# PanIsoGuard on a real published GM12878 ONT direct-RNA dataset (robustness).
# Adjudicates the wholegenome_multicaller caller-novel chains WITH the real ONT BAM (so the
# mapping axis sees real reads) and scans every novel junction for reference bias against
# NA12878's (= HG001) own SNV-consensus haplotype. Public data only.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIG="${PIG:-$HERE/../../build/panisoguard}"
WORK="${WORK:-/path/to/gm12878_mc}"        # the wholegenome_multicaller data dir: matrix.tsv + caller GTFs + ONT BAM
BAM="${BAM:-$WORK/gm.main.bam}"             # as wholegenome_multicaller/run.sh writes it
MATRIX="${MATRIX:-$WORK/matrix.tsv}"
HAP="${HAP:-/path/to/HG001/work}"          # per-chrom ref.fa/hap.fa (see benchmark/giab_cohort_rescue)
REFGTF="${REFGTF:-/path/to/gencode.v49.annotation.gtf}"
OUT="${OUT:-/tmp/gm12878_realdata}"; mkdir -p "$OUT/scan"
SCAN="$HERE/../hg002/scan_variant_axis.py"

# 1. build a GTF of the novel chains + a minimal SQANTI classification (keyed by PIG id)
python3 "$HERE/build_inputs.py" "$WORK" "$OUT" "$MATRIX"

# 2. adjudicate the novel chains WITH the real ONT BAM (mapping axis on real reads): the HiFi
#    defaults (rb) and the ONT preset (rb_ont; indel-near gate off), both with the 0.0.7
#    whole-chain read support
"$PIG" adjudicate --classification "$OUT/novel.cls" --isoforms-gtf "$OUT/novel_chains.gtf" \
    --ref-gtf "$REFGTF" --bam "$BAM" --out-prefix "$OUT/rb"
"$PIG" adjudicate --classification "$OUT/novel.cls" --isoforms-gtf "$OUT/novel_chains.gtf" \
    --ref-gtf "$REFGTF" --bam "$BAM" --platform ont --out-prefix "$OUT/rb_ont"
# 2b. how exact are ONT splice sites? chain reads at 0..10 bp tolerance (needs pysam)
"${PY:-python3}" "$HERE/chain_wobble.py" "$OUT/novel_chains.gtf" "$BAM" "$OUT/pig_ncallers.json" \
    --out "$OUT/chain_wobble.json"

# 3. reference-bias scan of every novel junction vs NA12878's own haplotype (chr1-22)
created=0; scanned=0
python3 -c "
import collections
by=collections.defaultdict(list)
for ln in open('$OUT/novel_chains.gtf'): by[ln.split('\t')[0]].append(ln)
for ch,ls in by.items(): open(f'$OUT/scan/{ch}.gtf','w').writelines(ls)
"
for chrom in chr{1..22}; do
  d="$HAP/$chrom"; g="$OUT/scan/$chrom.gtf"
  [ -f "$g" ] && [ -f "$d/hap.fa" ] || continue
  mkdir -p "$OUT/scan/$chrom"
  r=$(python3 "$SCAN" "$d/ref.fa" "$d/hap.fa" "$g" "$chrom" "$OUT/scan/$chrom" 2>&1 | tail -1)
  created=$((created + $(echo "$r" | grep -oE 'CREATED\(reference-bias\)=[0-9]+' | grep -oE '[0-9]+')))
  scanned=$((scanned + $(echo "$r" | grep -oE 'introns scanned: [0-9]+' | grep -oE '[0-9]+')))
done

# 4. emit the metrics envelope
python3 "$HERE/analyze.py" "$OUT/rb.attribution.jsonl" "$OUT/pig_ncallers.json" \
    "$created" "$scanned" --ont "$OUT/rb_ont.attribution.jsonl" --wobble "$OUT/chain_wobble.json" \
    --emit-metrics "$HERE/../results/gm12878_realdata/metrics.json"
