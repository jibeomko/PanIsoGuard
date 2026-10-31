# Relationship to multi-caller merge tools

Multi-caller consensus is established practice in long-read transcriptomics (LRGASP;
Pardo-Palacios *et al.* *Nat Commun* 2024). PanIsoGuard's `combine` is not a new merge
algorithm: it matches isoforms across callers by exact intron chain and feeds the result to
the adjudicator as one evidence axis. This page says where it overlaps with existing tools
and where it differs.

## Overlap

| Tool | What it does | Relation to PanIsoGuard |
|------|--------------|-------------------------|
| **gffcompare `-i`** (Pertea & Pertea 2020) | N-way comparison of GTFs by intron chain; `.tracking` lists which inputs share each transfrag | `combine` does the same grouping. On the same five caller sets it reproduces gffcompare's n_callers distribution and precision-recall curve ([benchmark/merge_comparison](../benchmark/merge_comparison)). The difference is practical: one htslib binary, no GTF round trip, output in the form the adjudicator reads. |
| **TAMA merge** (Kuo *et al.* 2020) | Fuzzy merge with configurable 5′, junction and 3′ tolerances | `combine` matches junctions exactly. On the chr22 set, exact matching, TAMA at 10 bp, and a 0–20 bp wobble sweep give the same precision at ≥ 3 callers (0.975–0.980), and no genuine novel isoform moves across the single/multi-caller boundary under any tolerance. |
| **Bambu NDR** (Chen *et al.* 2023) | A trained per-transcript novel-discovery rate | Single-caller confidence from read evidence; PanIsoGuard's class combines callers and evidence axes. Bambu can be one of the callers `combine` integrates. |
| **IsoQuant, FLAIR, ESPRESSO collapse** | Each caller's own transcript construction | These are the inputs. PanIsoGuard sits above them and does not call isoforms. |
| **SQANTI3 rescue and filter** | Recovers reference transcripts a caller missed; rule and ML QC filters | The word "rescue" means the opposite thing in the two tools (below); the QC overlap is in [relationship_to_sqanti3.md](relationship_to_sqanti3.md). |
| **IsoformSwitchAnalyzeR** (Vitting-Seerup 2019) | Downstream isoform-switch analysis | A consumer of filtered annotations, not an alternative. |

## "Rescue" in SQANTI3 and in PanIsoGuard

SQANTI3 rescue adds back a **reference** transcript that a caller failed to report.
PanIsoGuard's reference-bias rescue (`REFERENCE_BIAS`) keeps a **sample's novel call**
that looks novel only because the sample's genome differs from the reference: a haplotype
or a pangenome path makes the junction's splice site canonical. A circularity firewall holds
the call unless the haplotype or junction set comes from DNA-level data independent of the
sample's RNA. (Also noted in [relationship_to_sqanti3.md](relationship_to_sqanti3.md).)

## What PanIsoGuard adds

The consensus axis is shared practice. Two things in PanIsoGuard are not in the merge or QC
tools above:

1. **Reference-bias rescue with a circularity firewall.** A novel call is reclassified as
   reference bias only on independent haplotype or pangenome evidence, never on the
   sample's own RNA-derived variants. The yield is low and depends on how far the sample
   is from the reference: across the four GIAB genomes it is 31–37 junctions each (137 in
   total, all rescued under independent provenance, none rescued under circular-risk
   provenance, no rescue of a control junction in that scan;
   [benchmark/giab_cohort_rescue](../benchmark/giab_cohort_rescue)); on real Iso-Seq from
   two West African individuals (HG03516, HG02717; HPRC Release 2) it is 40 and 46
   caller-reported junctions, 77 of them novel to GENCODE, with the same outcome
   ([benchmark/refbias_cohort](../benchmark/refbias_cohort)); on the reference-grade GM12878
   it is 0 ([benchmark/gm12878_realdata](../benchmark/gm12878_realdata)). These are
   constructed and observational scans, not controlled rates; see each protocol's caveats.
2. **One verdict with its reasons.** Consensus, short-read, long-read mapping, haplotype
   and pangenome evidence go into one confidence class with a `verdict_reason` and a
   `rule_trace`, instead of a `.tracking` file, an NDR score and a QC table that the user
   reconciles by hand.
