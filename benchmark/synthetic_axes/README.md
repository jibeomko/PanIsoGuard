# Synthetic-axes controlled truth (BAM + variant axes)

Validates every evidence axis end-to-end through `panisoguard adjudicate` against
known truth, on a fully synthetic contig (no real genome, simulator, or caller).
Four labelled isoform categories, one test intron each:

| category | reference motif | short-read | BAM | haplotype | expected verdict |
|----------|-----------------|------------|-----|-----------|------------------|
| A genuine | canonical | supported | clean | = ref | `HIGH_CONF_NOVEL` |
| B noncanon-artifact | non-canonical | none | clean | = ref | `ARTIFACT` (noncanonical) |
| C mapping-artifact | canonical | none | low-MAPQ spanning | = ref | `ARTIFACT` (mapping) — needs `--bam` |
| D reference-bias | non-canonical | none | clean | canonical (SNV) | `PAN_REF_RESCUED_FALSE_NOVEL` — needs `--reference-haplotype` |

## Run

```bash
python gen.py work/ 40
samtools faidx work/ref.fa && samtools faidx work/hap.fa
samtools sort -O bam -o work/syn.bam work/syn.sam && samtools index work/syn.bam
PIG=../../build/panisoguard
C="--classification work/classification.tsv --isoforms-gtf work/caller.gtf --ref-gtf work/catalog.gtf --sj-tab work/real.SJ.tab"
$PIG adjudicate $C --out-prefix work/base
$PIG adjudicate $C --bam work/syn.bam --reference work/ref.fa --out-prefix work/bam
$PIG adjudicate $C --reference work/ref.fa --reference-haplotype work/hap.fa --haplotype-provenance external --out-prefix work/var
$PIG adjudicate $C --bam work/syn.bam --reference work/ref.fa --reference-haplotype work/hap.fa --haplotype-provenance external --out-prefix work/full
python score_axes.py work/truth.tsv work/full.adjudicated.tsv
```

## Result (40 isoforms/category)

| category (truth) | base (SR) | +BAM | +variant | full |
|------------------|-----------|------|----------|------|
| A genuine | HIGH_CONF_NOVEL | HIGH | HIGH | HIGH |
| B noncanon-artifact | ARTIFACT | ARTIFACT | ARTIFACT | ARTIFACT |
| C mapping-artifact | LOW_CONF_PARTIAL | **ARTIFACT (mapping)** | LOW | **ARTIFACT** |
| D reference-bias | ARTIFACT | ARTIFACT | **PAN_REF_RESCUED** | **PAN_REF_RESCUED** |

All 160 isoforms are classified per truth in the full configuration (0 errors).
The per-config deltas isolate each axis: the BAM axis alone moves C from
`LOW_CONF_PARTIAL` to `ARTIFACT` (mapping mechanism); the variant axis alone moves
D from `ARTIFACT` to `PAN_REF_RESCUED_FALSE_NOVEL`; neither axis perturbs the other
categories. This is a truth-based, ablation-style correctness proof for the BAM
mapping and variant/reference-bias axes, complementing the short-read-axis check in
[../controlled_truth](../controlled_truth).
