# Reference-bias rescue on real long-read RNA — a divergent-individual cohort

**The positive real-data application, replicated across two independent divergent
individuals.** On real PacBio Kinnex Iso-Seq from **two unrelated West African individuals**
— HG03516 (ESN, Esan in Nigeria) and HG02717 (GWD, Gambian Mandinka), both HPRC Release 2
with same-individual Iso-Seq **and** a HiFi de-novo phased assembly — PanIsoGuard's
reference-bias rescue fires on real caller-reported novel splice junctions and is
**perfectly specific across both genomes**. The reference-grade European control GM12878
(= NA12878) yields **0** ([gm12878_realdata](../gm12878_realdata)).

Each individual is processed identically (see [hg03516_refbias](../hg03516_refbias) for the
full single-individual deep-dive + the 5-angle adversarial audit): `flnc → minimap2
splice:hq → IsoQuant` for the novel junctions, and `assembly → minimap2 asm5 + paftools.js
call → SNV-consensus haplotype` for the personal genome.

## Cohort result

| individual | population | novel junctions scanned | **reference-bias** | novel vs GENCODE | rescued (indep-DNA) | false | firewall-held |
|------------|-----------|------------------------:|-------------------:|-----------------:|--------------------:|------:|--------------:|
| HG03516 | ESN (Nigeria) | 98,672 | **40** | 35 | 40/40 | **0** | 40/40 |
| HG02717 | GWD (Gambia) | 88,622 | **46** | 42 | 46/46 | **0** | 46/46 |
| **cohort** | West African ×2 | 187,294 | **86** | **77** | **86/86** | **0** | **86/86** |
| GM12878 (control) | European, reference-grade | 1,364 | **0** | — | — | — | — |

**Reading.** Across two independent divergent genomes, **86** caller-reported novel splice
junctions are reference bias — **77** of which a reference-only pipeline would report as
novel discoveries. Each is non-canonical on GRCh38 but canonical on that individual's own
HiFi-assembly haplotype via a real personal SNV at the splice dinucleotide. PanIsoGuard
exonerates **all 86 with 0 false rescues**, and the circularity firewall holds every one
under circular-risk provenance. The rescue fires consistently (~40–46 per genome) and is
**perfectly specific in both** — the cohort evidence that it is a reliable safeguard whose
yield **scales with divergence** (0 on reference-grade GM12878; cf. the
[giab_cohort_rescue](../giab_cohort_rescue), 137 across 4 personalized genomes).

## Caveats (carried from the [hg03516_refbias](../hg03516_refbias) audit; apply to both)

- **A principled read-support filter** drops any candidate with **0 exact-coordinate
  spanning reads** (a displaced-coordinate caller-model artifact) — **1 per individual**.
  All 86 reported have ≥ 2 exact spanning Iso-Seq reads (medians 5 / 9).
- **The divergent-vs-GM12878(0) contrast is directional, NOT a controlled rate** —
  the arms differ in scan denominator (~137× vs GM12878's 1,364), chemistry (PacBio Kinnex
  vs ONT direct-RNA), caller set, and variant source (HiFi assembly genome-wide vs GIAB
  v4.2.1 high-confidence SNPs that exclude segdups). The cohort + `giab_cohort_rescue` are
  the apples-to-apples evidence that yield scales with divergence.
- **"Independent-DNA provenance" = HiFi de-novo ASSEMBLY variant calls** (paftools.js
  placeholder QC, SNV-only → indel-mediated reference bias invisible), not short-read WGS
  genotyping. The firewall's independence is **molecule** independence (DNA assembly vs RNA),
  which holds; both arms share the GRCh38 + minimap2 backbone.

## Reproduce

```bash
# per individual (HG03516, HG02717): downloads HPRC flnc + assembly from public S3
WORK=/path/to/hg03516 ../hg03516_refbias/run.sh
WORK=/path/to/hg02717 ../hg03516_refbias/run.sh   # edit the S3 sample id
# then the cohort envelope:
analyze.py /path/to/gencode.v49.annotation.gtf HG03516:/path/to/hg03516 HG02717:/path/to/hg02717 \
  --emit-metrics ../results/refbias_cohort/metrics.json
```

Public HPRC Release 2 + GENCODE v49 + GRCh38 only; no private data. Result envelope:
[../results/refbias_cohort/metrics.json](../results/refbias_cohort/metrics.json).
