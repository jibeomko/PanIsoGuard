# Validation

PanIsoGuard is an adjudication layer, so "validation" means two distinct things,
kept separate on purpose:

- **What the `benchmark` / `ablate` subcommands report** — *non-redundancy* vs the
  SQANTI3 filter and *per-axis contribution*. These run on real data with **no truth
  labels**, so by construction they show disagreement and which axis drives which
  call — **not correctness**.
- **Ground-truth validation** — whether a reclassification is *right*. This needs
  labelled truth and is done by the dedicated protocols under
  [`../benchmark/`](../benchmark/), summarized below.

## Unit suite

`ctest` (htslib-only build): intron-chain fingerprint, interval index, all readers
on controlled fixtures (incl. SJ.tab strand discrimination), cross-reader coordinate
convergence, the full 2-axis rule projection (mapping override + variant/pangenome
circularity firewall), BAM read-level features, the splice-motif core, the pangenome
rescue (all-or-nothing + provenance firewall), and the non-redundancy/ablation math.

## Ground-truth validation (completed)

Each protocol generates a labelled truth set and scores PanIsoGuard's verdicts
against it. They validate the **adjudication logic** — given caller + SQANTI3 output
— at chr22 / SIRV / synthetic scale, on simulated reads and real human variants.

| Protocol | Setup | Result |
|----------|-------|--------|
| [controlled_truth](../benchmark/controlled_truth) | GENCODE chr22, incomplete-reference (hidden = genuine, shifted-exon = false); short-read axis | on decisive calls **precision(genuine)=1.000, specificity(false)=1.000**; genuine NIC-like isoforms honestly **abstained** (no novel junction to corroborate) |
| [synthetic_axes](../benchmark/synthetic_axes) | synthetic contig, 4 labelled categories × 40, run through `adjudicate` in 4 configs | full config classifies **all per truth (0 errors)**; per-config deltas isolate each axis (BAM → mapping artifact, haplotype → reference-bias rescue) |
| [bam_axis](../benchmark/bam_axis) | chr22 SQANTI-SIM with `--bam`; wire the previously-measured-but-unused **indel-near / soft-clip** read fractions into the mapping mechanism | the **indel-near** signal cleanly separates novel junctions (genuine max **0.190** vs false mean **0.350**); gated at 0.5 it catches **14/124 false at 100 % precision** (0 genuine flagged; 41/124 before the 0.0.6 rule that a junction feeds the gate only with ≥ 3 spanning reads). Wired ON vs OFF: **0 genuine lost**; 9 false calls relabelled `UNCONFIRMED → ARTIFACT` with an explicit mapping mechanism. Its earlier specificity lift (0.946 → 0.966) is now delivered by the `PARTIAL` fix without the BAM, so ON == OFF at the operating point. Soft-clip is **default-OFF** (terminal, not junction-proximal → risks demoting genuine novels on noisy reads; opt-in only) |
| [gm12878_realdata](../benchmark/gm12878_realdata) | the BAM axis + reference-bias rescue on a **real published GM12878 ONT direct-RNA** dataset (ENCFF440ZML) | **real-data robustness.** Soft-clip default-OFF **empirically vindicated** — on real ONT it would fire on **187/191 (98 %)** of junctions (adapter/poly-A). The HiFi-calibrated `indel_near` 0.5 gate **over-flags on ONT** (175/191 = 92 %; barely separates single 0.925 vs ≥2 0.732) while low-MAPQ/supplementary barely fire (6 %/2 %) → the BAM thresholds are **chemistry-dependent** (opt-in, re-calibrate per platform). The **consensus** axis separates the same data chemistry-independently (254 vs 31). Reference-bias rescue **0/1364** (GM12878 = NA12878 is reference-grade; value reserved for divergent genomes) The 0.0.7 **whole-chain read support works on real ONT with exact coordinates**: 94.5 % of splice sites among near-matching reads are exact, 155/285 chains have ≥ 3 exact full-chain reads (± 5 bp would add 20; no tolerance shipped); under `--platform ont` they become `NOVEL_PROBABLE` (148) instead of `UNCONFIRMED` under the HiFi gate (129; `ARTIFACT` before ruleset 0.0.5) |
| [end2end](../benchmark/end2end) | PBSIM3 → minimap2 → FLAIR3 → SQANTI3 → `adjudicate`, GENCODE chr22 (seeded, re-run with 0.0.7) | on the SQANTI-novel set (20 genuine / 249 false): **precision 1.000, specificity 1.000**, genuine recall among decisive calls 1.000 (19 `NOVEL_PROBABLE` through whole-chain reads, 1 abstain); 242 of 249 false novels rejected, 7 abstain — correct on *real* caller + SQANTI3 output |
| [hg002](../benchmark/hg002) | real GIAB HG002 v4.2.1 SNVs → personalized haplotype, `--haplotype-provenance wgs` | headline `REFERENCE_BIAS` validated on **real, independent (WGS-derived) variants** — non-circular; **0 false rescues** |
| [hg002_wholegenome](../benchmark/hg002) (yield) | the *same* scan extended to **whole-genome** HG002 (chr1–22, 597,781 GENCODE introns) | **rescue YIELD quantified**: reference bias at splice junctions is rare but real — **33** genuine reference-bias junctions genome-wide, **33/33 rescued** (100 % sensitivity), **0 false**, **33/33 held** by the firewall under circular-risk provenance. The rescue is a **high-specificity/sensitivity guardrail with a low base rate** (~1 in 18k introns), not a high-yield discovery engine |
| [giab_cohort_rescue](../benchmark/giab_cohort_rescue) | the rescue across a **four-individual GIAB cohort** (HG001/HG002/HG003/HG004) | **cohort impact**: 31–37 reference-bias junctions per personalized genome, **137 total → 137/137 rescued** (100 % sensitivity), **0 false**, **137/137 held** by the firewall. A consistent, perfectly-specific safeguard across individuals — the impact case for PanIsoGuard's differentiator |
| [refbias_cohort](../benchmark/refbias_cohort) | the rescue on **real PacBio Iso-Seq across TWO divergent West African individuals** (HG03516 ESN + HG02717 GWD; HPRC R2 — same-individual RNA + HiFi assembly each) | **POSITIVE real-data application, REPLICATED** (the differentiator firing on real RNA calls, perfectly specific across both genomes). **86 reference-bias junctions** (40 + 46; **77 novel vs GENCODE** — discoveries a reference-only pipeline would report), each non-canonical on GRCh38 but canonical on that individual's own assembly haplotype via a real SNV; all read-supported, **86/86 rescued (0 false)**, **86/86 firewall-held**. vs reference-grade GM12878 **0**. *Directional contrast, not a controlled rate.* See [hg03516_refbias](../benchmark/hg03516_refbias) for the single-individual deep-dive + 5-angle adversarial audit (read-support filter; assembly-SNV / no-indel / placeholder-QC caveats) |
| [refbias_phasing](../benchmark/refbias_phasing) | **non-circular** check of the rescue on the same two individuals: for heterozygous reference-bias junctions, are the long reads that use the junction from the haplotype on which it is canonical (phase from the trio-phased assembly)? | **27 of 28 testable junctions consistent** (572/575 phased junction reads from the motif haplotype; both haplotypes expressed at the loci). 1 contradicted (HG03516 chr3:184709997, junction reads from the non-canonical haplotype, all with a 2-bp insertion next to the donor — an alignment artifact the rescue does not screen for). 54 heterozygous, 34 homozygous (untestable) |
| [hg03516_refbias](../benchmark/hg03516_refbias) | the single-individual **deep-dive + adversarial audit** for HG03516 (one arm of the cohort above) | the audited HG03516 result (**40** reference-bias junctions, 35 novel, 40/40 rescued, 0 false, 40/40 held) with the full 5-angle audit that set the cohort's caveats: dropped 1 displaced-coordinate hit, exact-coordinate read support, and the directional (not controlled-rate) framing of the GM12878 contrast |
| [chain_support](../benchmark/chain_support) | the 0.0.7 **whole-chain long-read support** axis ON vs OFF on the five-caller chr22 set (same alignment, SQANTI-SIM truth), with 1 M short-read pairs and without | with short reads: FLAIR recall **0.545 → 0.745** at precision 0.993 (146 genuine / **0** false lifted), IsoQuant 111 / 0, ESPRESSO 132 / 1, Bambu 83 / 15, TALON 148 / 54; without short reads PanIsoGuard goes from deciding nothing to 63–92 % of isoforms at 0.7–1.3 % risk for FLAIR / IsoQuant / ESPRESSO. The signal removes caller-level false chains, not alignment-level ones; it is the callers' own reads, so `NOVEL_PROBABLE` only. `panisoguard decoys` (recombined annotated transcripts, a wrong combination of known junctions): in simulation wrong combinations reach ≥ 3 reads 0–0.1 % of the time and the deleted genuine transcripts 100 % (31/31); on real PacBio / Kinnex data 33–46 % of recombinants of two expressed isoforms do, with as many reads as their parents (unannotated combinatorial isoforms, not chimeras), so on real data the decoy rate is an upper bound and a per-dataset knob |
| [lrgasp_wtc11](../benchmark/lrgasp_wtc11) | the first check against **experimentally validated isoforms on real data**: LRGASP's 89 WTC11 RT-PCR + amplicon-sequencing targets, our own IsoQuant / SQANTI3 / PanIsoGuard run on the public WTC11 PacBio cDNA, ONT cDNA and Illumina reads | of the supported NIC/NNC targets IsoQuant reconstructed and that are still novel in GENCODE v49 (17 PacBio / 16 ONT): **0 `ARTIFACT`** on either platform, with or without short reads; 16 / 14 confident with short reads (the `UNCONFIRMED` ones have no short-read support at the junction in this one Illumina replicate); **without short reads all 17 / 16 are `NOVEL_PROBABLE` through whole-chain reads**; the 4 validated isoforms with SQANTI3's RT-switching sign carry it in the reason (`;rt_switch_sign`, ruleset 0.0.5) instead of being held. Recall-side only (79 positives / 7 negatives, none of the negatives reconstructed); 29 / 34 supported targets were not called by IsoQuant at all (caller recall, rarely detected targets). The GENCODE-reviewed sheet adds 26 negatives, but the caller reconstructed none of the 7 `rejected` models and the rest are annotated transcripts (`KNOWN` by pass-through); the one unsupported novel chain it did produce was held `UNRESOLVED`. SQANTI3's RTS sign is on 12–18 % of the validated isoforms (a measured false-sign rate) |
| [refbias_3way](../benchmark/refbias_3way) | the rescue against an **independent standard**: the same reads (HG03516, HG02717) through GRCh38 (A), an SNV-personalized diploid (B1), the HPRC phased assembly (B2) and PanIsoGuard post hoc on A (C); junction truth = where the **same reads** splice on the personal genome | rule level: **31/32** and **37/39** SNV-class junctions on B1 (31/35, 38/43 on B2), precision 0.88–0.91 (the misses are reads from the non-canonical haplotype); verdict level **12/32** and **21/39**, precision 12/13 and 21/22 (the gap is the all-junctions rule and non-adjudicated categories, not the criterion); shifted-known junctions (indel/SV) never caught, by construction; firewall 0 rescues under `unknown`. Cost: personalized pipelines 7–12 h / 38–78 GB vs `adjudicate` 3–6 min per run / 0.5 GB |
| [sirv_multicaller](../benchmark/multicaller) | the multi-caller consensus on a **second reference** — Lexogen SIRV-Set4 (dense, multi-contig) | **generalization**: single-caller novel precision **0.037** → ≥ 2-caller **1.000** vs the SIRV hidden-chain truth, mirroring the chr22 result (0.012 → ≥3 0.980). The consensus signal is a property of multi-caller data, not of one reference/simulator |
| [sirv](../benchmark/sirv) | Lexogen SIRV-Set4 spike-in control, dense overlapping isoforms | specificity(false)=0.944, recall=1.000; the 3 FPs are FLAIR mis-collapses of *individually real* junctions — a documented short-read-axis limitation |
| [calibration](../benchmark/calibration) | class → empirical P(genuine), train/test split | **well-calibrated**: controlled Brier 0.0000 / ECE 0.0013; end-to-end Brier 0.0121 / ECE 0.0114 |
| [multicaller](../benchmark/multicaller) | **five** callers (FLAIR + IsoQuant + Bambu + ESPRESSO + TALON) on one shared alignment, integrated by `combine`, scored vs SQANTI-SIM truth | single-caller novel precision **0.012** (per-caller precision 0.40–0.98); PR curve over "≥ k callers" peaks at **≥3 (P 0.980, R 0.910, F1 0.944)** — the optimal consensus threshold scales with caller count. Engine reproduces it: `adjudicate --caller-support` (long-read only) takes FLAIR from **0 confident novels (all UNRESOLVED) → 0.975 precision** by promoting cross-caller-agreed novels |
| [wholegenome_multicaller](../benchmark/wholegenome_multicaller) | **real** public GM12878 whole-genome (IsoQuant + Bambu + ESPRESSO), no truth → orthogonal **canonical-motif** validation + **decision-impact** framing | *Honest, small-N.* Disagreement holds at genome scale (**254 single-caller vs 31 consensus** of 285 novel chains). Consensus novels **100 %** canonical vs single-caller **97.6 %** — direction correct but small lift (+0.024, ceiling effect: stringent production callers are already clean). **Decision impact:** the reported novel count swings **31 → 285 (9.2×)** purely on integration policy — **89 % (254/285)** of the naive-union novels are caller-*private* (no second caller corroborates). On real data the value is the reproducible high-confidence core + making the 9× policy choice explicit/auditable, not bulk artifact removal |
| [merge_comparison](../benchmark/merge_comparison) | head-to-head vs the incumbent merge tools on the **same** 5-caller chr22 set: `combine` vs **gffcompare -i**, **TAMA**, and a 0–20 bp wobble sweep | **`combine` ≡ gffcompare -i exactly** (same n_callers distribution + PR curve) — a clean re-implementation, not a novel merge. The consensus ≥3 precision is **matcher-robust** (0.975–0.980 across exact / gffcompare / TAMA / wobble), and **0** genuine novels cross the single↔multi boundary under any fuzzy matcher → the consensus gate is a property of the **data**, not of exact matching. Honest: consensus is established practice; the **reference-bias rescue** is the differentiator |
| [sqanti_sim](../benchmark/sqanti_sim) | canonical SQANTI-SIM (GENCODE chr22): delete 769 transcripts → PBSIM3 HiFi → FLAIR → SQANTI3 → `adjudicate` | genuine-novel detection **precision 0.989, specificity 0.966, AUPRC 0.971** (0.831 = positive rate, i.e. a random ranking — not a competing method; see [sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h) for the real baselines); NNC recall 1.000; known→KNOWN 1.000; **0 genuine→ARTIFACT**. Moderate overall recall = honest abstention on NIC-combinatorial / ISM-partial, not misclassification |
| [pangenome](../benchmark/pangenome) (GATE-1) | real **HPRC v1.1** MC GRCh38 chr22 graph → `vg deconstruct` → population-deletion junctions; sample **out-of-graph** (non-circular) | **0 false rescues** on real FLAIR novel junctions; rescue **fires** when a novel intron is a real population deletion (population provenance); circularity firewall **holds** all rescues `UNRESOLVED` under circular-risk provenance |

**Scope of the reference-bias numbers (read before citing them).** In `hg002`,
`hg002_wholegenome`, `giab_cohort_rescue`, `hg03516_refbias` and `refbias_cohort`, a
junction is labelled `CREATED` (reference bias) when it is non-canonical on GRCh38 and
canonical on the individual's haplotype — the same criterion the variant rescue applies.
So "N/N rescued, 0 false, N/N firewall-held" shows the rule is implemented as specified
and that the provenance gate works; it is **not** independent evidence that these
junctions are reference-bias artifacts. The motif-creating base is always intronic, so
reads that use the junction never contain it. The independent check is
[refbias_phasing](../benchmark/refbias_phasing): for heterozygous junctions, do the long
reads that use the junction carry, at the gene's other heterozygous SNVs, the alleles of
the haplotype that makes it canonical? In HG03516 + HG02717, 27 of 28 testable junctions
pass (572 of 575 phased junction reads from that haplotype, while both haplotypes are
expressed at the loci); the one failure looks like an alignment artifact. The 34
homozygous junctions cannot be phased and remain unchecked. Yield is low: ~30–45
junctions per divergent genome.

**Accuracy against real baselines.** On the SQANTI-SIM truth set,
[sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h) compares PanIsoGuard with the
SQANTI3 rules / ML filter and a one-line short-read rule, with oracle and simulated
Illumina junctions. At ranking PanIsoGuard ties the one-line rule, with or without the BAM;
at the operating point it is more precise and, since the 0.0.7 whole-chain read support,
recalls 0.75 of the genuine isoforms with the BAM (0.55 without). With no short reads and the
BAM it ranks above SQANTI3's rules filter (0.953 vs 0.909); without either it confirms
nothing. See there, including the risk–coverage points.

Together these cover the validated short-read, BAM mapping, and variant axes,
the full pipeline (real FLAIR/SQANTI3), the headline reference-bias rescue on
real human variation, calibration, and multi-caller integration. The file-based
pangenome axis is now validated on the **real HPRC v1.1 chr22 graph** (GATE-1, below):
0 false rescues on real FLAIR novel junctions, correct rescue on real population
deletions, and the circularity firewall holding under circular-risk provenance.

## Tracked result artifacts

These numbers are no longer prose-only. The metric for every protocol is committed
under [`../benchmark/results/`](../benchmark/results/) as a schema-checked
`metrics.json`, indexed by [`MANIFEST.tsv`](../benchmark/results/MANIFEST.tsv):

- **`status=tracked`, self-contained** — `controlled_truth` and `synthetic_axes` are
  regenerated and drift-checked by [`collect.py`](../benchmark/results/collect.py), and
  the same drivers run in CI as the `integration_controlled_truth` /
  `integration_synthetic_axes` CTests, so a metric and its pass/fail assertion share one
  code path.
- **`status=tracked`, heavy pipeline** — `sqanti_sim`, `end2end`, `hg002`, and
  `pangenome` carry the result of a real run with the current binary, produced by their
  documented `command`. The SQANTI-SIM **threshold sweep**
  ([`sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv), 30 configs via
  [`sweep.py`](../benchmark/sqanti_sim/sweep.py)) closed the calibration gate above;
  `end2end` gives specificity 1.000 / precision 1.000 on a seeded pbsim run; `hg002` reproduces the
  real-variant rescue with 0 false rescues; `pangenome` is the GATE-1 result on the real
  HPRC v1.1 chr22 graph (its `gate1_check` sensitivity/firewall half is the
  `integration_pangenome_gate1` CTest).
- **`status=transcribed_pending_tracked_run` / `pending`** — reserved for protocols
  whose number is still only transcribed from prose or not yet produced (none of the
  core axes are in this state now).

## What remains

| Item | Status | Why |
|------|--------|-----|
| **Pangenome real-graph rescue (GATE-1)** | **done (chr22)** | Validated on the real HPRC v1.1 Minigraph-Cactus chr22 graph (`vg deconstruct` → population-deletion junctions → `--pangenome-junctions`); see [pangenome](../benchmark/pangenome). Remaining: genome-wide scale, insertion/inversion-based junctions, and the in-process GBZ traversal (`-DWITH_PANGENOME_LIB`). |
| **LRGASP real data** ([lrgasp/](../benchmark/lrgasp)) | blocked | The LRGASP pre-run caller GTFs are **Synapse-gated**. The multi-caller capability is already validated on genuine **five-caller** output (FLAIR + IsoQuant + Bambu + ESPRESSO + TALON), truth-scored (see `multicaller`). |
| **HG002 long-read RNA** | blocked | No clean public HG002/GM24385 long-read RNA-seq dataset was found (ENCODE/ENA empty). The variant axis is instead validated on real HG002 *variants* + the full pipeline on `end2end`/`sirv`. |

**Calibration caveat.** `config/rules.default.toml` thresholds are conservative
defaults. They are well-calibrated on the run truth sets above (low ECE/Brier) and
the SQANTI-SIM AUPRC (0.971) confirms strong genuine-vs-false ranking — on an oracle
short-read set, and no better than a one-line short-read rule (see
[sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h)). A threshold
**sweep** on the SQANTI-SIM v49 chr22 truth set
([`benchmark/results/sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv),
produced by [`benchmark/sqanti_sim/sweep.py`](../benchmark/sqanti_sim/sweep.py)) shows
AUPRC is **robust (0.970–0.971)** across the `sj_min_uniq_reads × require_canonical_motif
× perc_A_degradation` grid and the shipped `default-0.0.1` config is within 1e-4 of the
grid-best — the only degradation is the expected support cliff once `sj_min_uniq_reads`
exceeds the data's coverage (AUPRC → 0.919). So the conservative defaults are
**near-optimal there**; they are not yet swept on additional datasets. Treat the
confidence classes as labels, not a tuned probability and not an ordered scale: on SQANTI-SIM
without the BAM, `UNRESOLVED` is more precise than `UNCONFIRMED` at every short-read depth
([benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h), `per_class.py`).
