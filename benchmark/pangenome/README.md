# Pangenome reference-bias rescue — GATE-1 (real HPRC graph)

Validates the **file-based pangenome reference-bias axis** on the real
[HPRC v1.1 Minigraph-Cactus GRCh38 graph](https://github.com/human-pangenomics/hpp_pangenome_resources),
satisfying **GATE-1** in [../../docs/decision_engine.md](../../docs/decision_engine.md).

## Idea

A population **deletion** carried on a pangenome haplotype is the genomic form of a
reference-bias "novel intron": an RNA read from a haplotype that carries the deletion,
aligned to the linear reference, spans the deleted span as an apparent novel intron,
though no new splice site exists. So a novel junction that exactly matches a population
deletion is reference bias, not new splicing — and only **independent** (non-circular)
evidence if it comes from haplotypes *other than the sample being adjudicated*.

## Pipeline

```
HPRC v1.1 MC GRCh38 graph (hprc-v1.1-mc-grch38.gbz)
  -> vg chunk (chr22)  ->  vg deconstruct (GRCh38 chr22 path)  ->  chr22.graph.vcf
  -> derive_junctions.py  (population deletions >= 50 bp -> graph-supported junctions,
                           both strands, n_haplotypes = AC)
  -> panisoguard adjudicate --pangenome-junctions ... --pangenome-provenance population
```

**Non-circular by construction:** the adjudicated sample is **out of the graph** —
`vg paths --list-samples` on the HPRC graph lists 47 assembly samples and does **not**
include HG002/NA24385, so any graph support comes from *other* population haplotypes.

## Validation (GENCODE v49 chr22; HPRC v1.1 chr22 → 1,520 graph junctions)

| Axis | Setup | Result |
|------|-------|--------|
| **Specificity** | real FLAIR chr22 isoforms (1,579; 730 genuine-novel) run through the HPRC junctions, provenance `population` | **0 false rescues** — no genuine novel junction coincides with a population deletion (0 rescued overall) |
| **Sensitivity** | constructed isoforms whose novel intron **is** a real HPRC deletion ([gate1_check.py](gate1_check.py)) | **6/6 rescued** `PAN_REF_RESCUED_FALSE_NOVEL` under `population` provenance |
| **Circularity firewall** | same isoforms, default `unknown` provenance | **6/6 held `AMBIGUOUS`** with `circularity_flag=true` (never promoted) |
| **Exactness** | near-miss (5 bp off) and off-graph control introns | **0 false rescues** (exact junction match required) |

So the pangenome axis rescues an apparent novel junction **only** when it is a real
population deletion realizable on the graph, **only** when the provenance is declared
independent, and **never** spuriously on genuine novel junctions.

## Reproduce

```bash
# 1) graph -> deconstruct VCF (heavy; needs vg + the HPRC gbz)
vg chunk -x hprc-v1.1-mc-grch38.gbz -p GRCh38#0#chr22 ... > chr22.vg
vg deconstruct -p GRCh38#0#chr22 chr22.vg > chr22.graph.vcf
# 2) derive graph-supported junctions
python3 derive_junctions.py chr22.graph.vcf chr22.pangenome_junctions.tsv 50
# 3a) specificity on a real chr22 FLAIR set (see ../sqanti_sim)
panisoguard adjudicate --classification flair_classification.txt --isoforms-gtf flair.isoforms.gtf \
  --ref-gtf chr22_modified.gtf --pangenome-junctions chr22.pangenome_junctions.tsv \
  --pangenome-provenance population --out-prefix spec
# 3b) sensitivity + firewall (self-contained on the committed sample)
python3 gate1_check.py /path/to/panisoguard testdata/sample_junctions.tsv
```

`testdata/sample_junctions.tsv` is a 50-junction slice of the real HPRC-derived set
(genomic coordinates only) so `gate1_check.py` runs self-contained (it is also the
`integration_pangenome_gate1` CTest).

## Whole-genome specificity on real public data (GM12878)

Extends the GATE-1 specificity from chr22 to **whole-genome** on a **real, fully public**
sample: GM12878 (NA12878) ONT direct-RNA from ENCODE (ENCSR368UNC), classified with
SQANTI3 vs GENCODE v49 → **13,017 novel isoforms**, adjudicated against the
**whole-genome** HPRC v1.1 deletion junctions (`vg deconstruct` of all GRCh38
chr1-22/X/Y paths → 71,435 junctions). NA12878 is out of the HPRC graph (non-circular).

```text
pangenome false rescues:                0 / 13,017 novel isoforms
novel junctions within +/-2 bp of a population deletion:  0   (genuinely no coincidence)
BAM mapping axis: novel isoforms flagged as mapping artifacts:  390
```

**Honest reading.** Zero false rescues confirms whole-genome specificity on real data.
A proximity scan shows the novel junctions are *not* near population deletions at all —
they are real splice sites, not reference bias. So **reference-bias-via-population-deletion
is rare in typical samples**, and the pangenome axis is best understood as a
**high-specificity guardrail** (it will not over-promote) rather than a high-yield
discovery axis. Recorded in
[../results/pangenome_gm12878](../results/pangenome_gm12878); fully reproducible from
public ENCODE/HPRC inputs.

## Scope / what remains

Validated for **specificity** at chr22 (constructed sensitivity/firewall) and
**whole-genome** scale (real public GM12878) on real HPRC v1.1. The derivation uses graph
**deletions** (the dominant reference-bias mechanism); insertion/inversion-based
junctions and an **in-process** graph traversal (`-DWITH_PANGENOME_LIB`, gbwtgraph/GBZ)
remain future work — today the axis consumes the pre-extracted junction file.
