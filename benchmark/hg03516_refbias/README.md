# Reference-bias rescue on real long-read RNA from a divergent individual (HG03516)

**The positive real-data application case for PanIsoGuard's differentiator.** On real
PacBio Kinnex Iso-Seq from **HG03516** — a divergent African individual (ESN, Esan in
Nigeria; HPRC Release 2) — the reference-bias rescue **fires on real caller-reported novel
splice junctions**: junctions a reference-only pipeline reports as novel discoveries that
are actually explained by the individual's *own* genome. On the reference-grade European
sample GM12878 (= NA12878) the same scan finds **0** ([gm12878_realdata](../gm12878_realdata)).

Everything is public: HPRC Release 2 pairs, for the **same** individual, a PacBio Kinnex
Iso-Seq `flnc.bam` and a HiFi de-novo **phased assembly** (no dbGaP/EGA).

## What was run

```
flnc Iso-Seq (11.34M reads) ── minimap2 splice:hq ─→ IsoQuant ─→ 37,568 novel models
HiFi assembly (pat + mat) ── minimap2 asm5 + paftools.js call vs GRCh38 ─→ SNV-consensus
                              personal haplotypes (length-matched)
  → scan every novel junction: non-canonical on GRCh38 AND canonical on HG03516's haplotype?
  → adjudicate the variant axis: --haplotype-provenance wgs (rescue) and unknown (firewall)
```

## Result (independently adversarially audited)

| | HG03516 (divergent, ESN) | GM12878 (reference-grade) |
|---|:---:|:---:|
| unique novel junctions scanned | 98,672 | 1,364 |
| **reference-bias junctions** | **40** | **0** |
| &nbsp;&nbsp;novel vs GENCODE v49 | **35** | — |
| &nbsp;&nbsp;trace to a real personal splice-site SNV | 40/40 | — |
| &nbsp;&nbsp;read-supported (≥2 exact spanning reads, median 5) | 40/40 | — |
| **rescued** `PAN_REF_RESCUED_FALSE_NOVEL` (independent-DNA prov.) | **40/40, 0 false** | — |
| **firewall-held** AMBIGUOUS (circular-risk prov.) | **40/40** | — |

**Reading.** On a divergent individual's real transcriptome, **40** novel splice junctions
are reference bias — **35** of which a reference-only pipeline would report as novel
discoveries. Each is non-canonical on GRCh38 but canonical on HG03516's own haplotype via a
real assembly SNV at the splice dinucleotide (17 homozygous, 23 heterozygous). PanIsoGuard
exonerates all 40 with **0 false rescues**, and the circularity firewall holds every one
when provenance is circular-risk. This is the rescue's value made concrete on real data, and
it **scales with divergence** — 0 on reference-grade GM12878, 40 here, consistent with the
[giab_cohort_rescue](../giab_cohort_rescue) (137 across 4 personalized genomes).

## Honest caveats (from a 5-angle independent adversarial audit)

1. **One hit was dropped.** `chr17:3681978-3685514` (P2RX5) is a **4 bp-displaced IsoQuant
   model coordinate** with **0** exact spanning reads — its 384 reads splice 4 bp away to a
   site already GT-AG canonical on GRCh38. It is excluded from the 40 (it was an annotated
   `non_canonical_polymorphism`, the most circular case). The 40 reported all re-derive
   cleanly with ≥ 2 exact spanning reads.
2. **The HG03516-vs-GM12878 "0" is directional, NOT a controlled rate.** The two arms differ
   ~72× in scan denominator (98,672 vs 1,364 junctions), plus chemistry (PacBio Kinnex vs
   ONT direct-RNA), caller set (IsoQuant vs multi-caller), and variant source (HiFi assembly,
   genome-wide incl. repeats, vs GIAB v4.2.1 high-confidence SNPs that exclude segdups). At
   HG03516's per-junction rate GM12878's 0 is expected (E = 0.55, P(0) = 0.58). So GM12878's
   0 **illustrates** the divergence effect but does not, alone, prove a low base rate; the
   **apples-to-apples** evidence is `giab_cohort_rescue` (137 across 4 genomes).
3. **"Independent-DNA provenance" = HiFi de-novo ASSEMBLY variant calls**, not short-read WGS
   genotyping. paftools.js writes placeholder genotype/quality (GT 1/1, QUAL 60) and the
   consensus is **SNV-only** (indel-mediated reference bias is invisible by construction). The
   firewall's independence is **molecule** independence (DNA assembly vs RNA), which holds —
   but both arms share the GRCh38 + minimap2 backbone. *Mitigations confirmed by the audit:*
   40/40 independent exact read support, and **0/40** hits in segdup / paralog / low-complexity
   regions, with each restoring SNV collinear on a single assembly contig.

## Reproduce

```bash
WORK=/path/to/hg03516 ./run.sh    # downloads HPRC HG03516 flnc + assemblies (public S3)
```

Public data only (HPRC Release 2 HG03516 Iso-Seq + assembly; GENCODE v49; GRCh38). Result
envelope: [../results/hg03516_refbias/metrics.json](../results/hg03516_refbias/metrics.json).
