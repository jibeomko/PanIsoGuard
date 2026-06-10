# Relationship to multi-caller merge / consensus tools

Multi-caller consensus is **established practice** in long-read transcriptomics (LRGASP;
Pardo-Palacios *et al.* *Nat Commun* 2024). PanIsoGuard's `combine` is **not a novel
merge algorithm** — it operationalizes that practice as one evidence axis feeding the
adjudicator. This page maps the overlap honestly and says where PanIsoGuard is
genuinely differentiated vs merely duplicative.

## Overlap map

| Tool | What it does | Overlap with PanIsoGuard | Honest verdict |
|------|--------------|--------------------------|----------------|
| **gffcompare `-i`** (Pertea & Pertea 2020) | N-way comparison of multiple GTFs by intron chain; `.tracking` shows which inputs share each transfrag | **`combine` re-implements this.** Head-to-head, `combine` reproduces gffcompare's grouping **exactly** (same n_callers distribution + PR curve, see [benchmark/merge_comparison](../benchmark/merge_comparison)) | **Duplicative by design** — `combine` is a clean htslib-native re-implementation (no Perl, no GTF round-trip) that feeds the adjudicator directly. Not a new algorithm; we credit gffcompare. |
| **TAMA merge** (Kuo *et al.* 2020) | Fuzzy merge of transcript models with configurable 5′/junction/3′ wobble | `combine` uses **exact** junction matching, not wobble | **Matcher-robust, not a gap.** On the same data, exact vs TAMA(10 bp) vs a 0–20 bp wobble sweep give the same ≥3 consensus precision (0.975–0.980), and **0** genuine novels cross the single↔multi boundary under any fuzzy matcher. Exact matching does not lose the agreements the consensus gate depends on. |
| **Bambu NDR** (Chen *et al.* 2023) | A trained per-transcript novel-discovery-rate (probabilistic confidence) | Both assign confidence to novel transcripts | **Orthogonal & single-caller.** NDR scores one caller's models from read evidence; PanIsoGuard's confidence is **cross-caller + cross-axis**. Complementary — Bambu's NDR can be one of the callers PanIsoGuard integrates. |
| **IsoQuant / FLAIR / ESPRESSO built-in collapse** | Each caller's own transcript construction / quantification | These are the **inputs** PanIsoGuard consumes | **Upstream.** PanIsoGuard is caller-agnostic and sits above them. |
| **SQANTI3 rescue / filter** | Recovers lost **reference** transcripts; ML/rule QC filter | Name clash on "rescue" (see below); QC overlap covered in [relationship_to_sqanti3.md](relationship_to_sqanti3.md) | **Different direction.** See the disambiguation note. |
| **IsoformSwitchAnalyzeR** (Vitting-Seerup 2019) | Downstream isoform-switching / functional consequence analysis | none on the calling side | **Downstream consumer**, not a competitor. |

## "Rescue" means the opposite thing here vs SQANTI3

A frequent confusion: **SQANTI3 rescue** recovers a **reference** model that a caller
*failed to report* (it adds a known transcript back). **PanIsoGuard's reference-bias
rescue** (`PAN_REF_RESCUED_FALSE_NOVEL`) does the opposite — it **exonerates a sample's
novel call** that looks novel only because of *reference bias* (a haplotype/pangenome
path makes the apparent novel junction explainable without new splicing), and it is gated
by a **circularity firewall** so a rescue never rests on the sample's own RNA-derived
data. (Also noted in [relationship_to_sqanti3.md](relationship_to_sqanti3.md).)

## Where PanIsoGuard is actually differentiated

The consensus axis is **sound and validated** (above) but is **not** the differentiator —
it is shared practice done correctly. PanIsoGuard's genuine contributions, which no
surveyed merge/QC tool provides, are:

1. **Reference-bias rescue + circularity firewall** — reclassifying a *sample* novel call
   as reference bias on **independent** (variant / pangenome) evidence, never on circular
   RNA-derived data. No surveyed merge/QC tool does this. *Honestly scoped:* it is a
   **high-specificity, high-sensitivity guardrail with a low base rate**, not a high-yield
   discovery engine — reference bias at splice junctions is rare but *recurrent per genome*:
   across a **four-individual GIAB cohort** (HG001/2/3/4) it finds 31–37 genuine cases each
   (**137 total → 137/137 rescued, 0 false, 137/137 firewall-held**); population-deletion
   reference bias in typical cell lines is ~0. Its value is *correctness* (it never
   over-promotes, and catches every genuine case), and it rises for non-reference /
   personalized-genome samples. **On real long-read RNA from TWO divergent
   West African individuals (HG03516 ESN + HG02717 GWD, HPRC R2 — same-individual Iso-Seq +
   HiFi assembly each) the rescue fires on real caller calls and is perfectly specific in
   both: 40 + 46 = 86 IsoQuant novel junctions are reference bias (77 novel vs GENCODE —
   discoveries a reference-only pipeline would report), 86/86 rescued with 0 false and the
   firewall holding all 86; the reference-grade control GM12878 yields 0** (a directional
   contrast, not a controlled rate — see the protocol's caveats). See
   [benchmark/refbias_cohort](../benchmark/refbias_cohort) (the positive real-data case,
   replicated) and [benchmark/hg03516_refbias](../benchmark/hg03516_refbias) (its audited
   deep-dive), [benchmark/giab_cohort_rescue](../benchmark/giab_cohort_rescue),
   [benchmark/hg002](../benchmark/hg002), and [benchmark/pangenome](../benchmark/pangenome).
2. **One unified, auditable verdict** — consensus, short-read, long-read-mapping, variant,
   and pangenome evidence projected into a single confidence class with a machine-readable
   `rule_trace`, rather than leaving the user to reconcile a `.tracking` file, an NDR
   score, and a QC table by hand.

In short: PanIsoGuard's pitch should lead with the **rescue + firewall + unified
adjudication**, and present multi-caller consensus as a **correctly-operationalized,
matcher-robust standard input** — not as a novel merge.
