# Multi-caller `combine` on real caller output

Validates `panisoguard combine` (intron-chain integration of multiple callers'
differently-named isoforms) on **genuine multi-caller output** — without the
Synapse-gated LRGASP pre-run GTFs — by running a second caller (IsoQuant) on the
same simulated alignment used by the [end-to-end](../end2end) benchmark and
integrating it with FLAIR.

## Run (after the end-to-end FLAIR run)

```bash
# same simulated BAM (aln.bam) + reduced annotation from ../end2end
isoquant --reference chr.fa --genedb reduced.gtf --bam aln.bam \
         --data_type pacbio_ccs -o isoquant -t 8
panisoguard combine \
  --gtf flair:flair.isoforms.gtf \
  --gtf isoquant:isoquant/OUT/OUT.transcript_models.gtf \
  --ref-gtf reduced.gtf --out combine_matrix.tsv
```

## Result (GENCODE v49 chr22 simulation)

FLAIR produced 397 isoforms, IsoQuant 157; `combine` integrated them into 426 unique
intron chains:

```
n_callers=2 (agreed by BOTH callers): 128
n_callers=1 (single caller):          298
known vs reduced annotation: 143   novel: 283
```

The two callers use completely different native ID schemes, yet identical splice
chains are merged onto one stable `PIG.NNNNNN` id with both native IDs preserved:

```
pig_id      n_callers  callers          novelty  native_ids
PIG.000003  2          flair,isoquant   novel    flair=S_5022; isoquant=transcript19.chr22.nic
PIG.000005  2          flair,isoquant   known    flair=ENST00000215794.8; isoquant=ENST00000215794.8
```

This confirms the caller-running-ID recognition + integration device on real,
independently-named multi-caller output (complements the synthetic cross-caller unit
test).
