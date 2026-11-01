# Validation

PanIsoGuard is an adjudication layer, so two different things get called validation here
and are kept apart:

- What a run on real data shows: how PanIsoGuard's verdicts differ from the SQANTI3 filter's
  and which evidence axis drives which call. Without truth labels that is disagreement, not
  correctness.
- Ground-truth validation: whether a reclassification is right. This needs labelled truth
  and is done by the protocols under [`../benchmark/`](../benchmark/), summarised below.

## Unit suite

`ctest` (htslib-only build) covers the intron-chain fingerprint, every reader on controlled
fixtures (including SJ.tab strand discrimination), cross-reader coordinate convergence, the
two-axis rule projection (mapping override, variant and pangenome circularity firewall),
BAM read-level features, the splice-motif core, the pangenome rescue (all-or-nothing plus the
provenance firewall), the decoy generator, and the input and output guards (write failures,
an output that is an input, malformed GTF/BED12, config values out of range).

## Ground-truth validation

Each protocol builds a labelled truth set and scores PanIsoGuard's verdicts against it. They
test the adjudication logic, given caller and SQANTI3 output, at chr22, SIRV and synthetic
scale, on simulated reads and on real human variants.

| Protocol | Setup | Result |
|----------|-------|--------|
| [controlled_truth](../benchmark/controlled_truth) | GENCODE chr22 with an incomplete reference (hidden transcripts = genuine, shifted-exon = false); short-read axis | on decisive calls precision(genuine) 1.000, specificity(false) 1.000; genuine NIC-like isoforms are held (no novel junction to corroborate) |
| [synthetic_axes](../benchmark/synthetic_axes) | synthetic contig, 4 labelled categories × 40, run through `adjudicate` in 4 configurations | the full configuration classifies all per truth (0 errors); the per-configuration differences isolate each axis (BAM → mapping artifact, haplotype → reference-bias rescue) |
| [bam_axis](../benchmark/bam_axis) | chr22 SQANTI-SIM with `--bam`; the indel-near read fraction wired into the mapping mechanism | indel-near separates novel junctions (genuine max 0.191 vs false mean 0.129 over junctions with ≥ 3 spanning reads); gated at 0.5 it catches 14 of 124 false junctions with no genuine one flagged (before the 0.0.6 spanning-read gate, which counted every junction with ≥ 1 read: false mean 0.350, 41 of 124 caught) |
| [gm12878_realdata](../benchmark/gm12878_realdata) | the BAM axis and the reference-bias rescue on the published GM12878 ONT direct-RNA data (ENCFF440ZML) | the HiFi-calibrated `indel_near` gate fires on 175 of 191 ONT junctions, so `--platform ont` turns it off; the reference-bias scan finds 0 of 1,364 (a reference-grade genome); whole-chain reads work on ONT with exact coordinates (155 of 285 chains at ≥ 3 reads) |
| [end2end](../benchmark/end2end) | PBSIM3 → minimap2 → FLAIR3 → SQANTI3 → `adjudicate`, GENCODE chr22 (seeded, re-run with 0.0.7) | on the SQANTI-novel set (20 genuine, 249 false): precision 1.000, specificity 1.000; 19 genuine isoforms `NOVEL_PROBABLE` through whole-chain reads, 1 held; 242 of 249 false novels rejected, 7 held |
| [hg002](../benchmark/hg002) | GIAB HG002 v4.2.1 SNVs → personalized haplotype, `--haplotype-provenance wgs` | the `REFERENCE_BIAS` rescue on real, WGS-derived (non-circular) variants: the one genuine case rescued, 0 false rescues |
| [hg002_wholegenome](../benchmark/hg002) | the same scan over the whole genome (chr1–22, 597,781 GENCODE introns) | 33 reference-bias junctions genome-wide; 33 of 33 rescued, 0 false, 33 of 33 held by the firewall under circular-risk provenance |
| [giab_cohort_rescue](../benchmark/giab_cohort_rescue) | the same scan on the four GIAB genomes (HG001–HG004) | 31–37 reference-bias junctions per genome, 137 in total; 137 rescued, 0 false, 137 held |
| [refbias_cohort](../benchmark/refbias_cohort) | the rescue on real PacBio Iso-Seq from two West African individuals (HG03516 ESN, HG02717 GWD; HPRC Release 2, same-individual RNA and HiFi assembly) | 86 reference-bias junctions (40 + 46; 77 novel to GENCODE); 86 rescued, 0 false, all 86 held under circular-risk provenance; the aggregator was corrected on 2026-10-31 and the recount gives the same numbers |
| [refbias_phasing](../benchmark/refbias_phasing) | a non-circular check on the same two individuals: for heterozygous reference-bias junctions, do the reads that use the junction come from the haplotype on which it is canonical (phase from the trio-phased assembly)? | 27 of 28 testable junctions consistent (572 of 575 phased junction reads from the motif haplotype, with both haplotypes expressed at the loci); the one failure looks like an alignment artifact |
| [hg03516_refbias](../benchmark/hg03516_refbias) | the single-individual analysis and review behind the cohort protocol | 40 reference-bias junctions (35 novel), 40 rescued, 0 false, 40 held, after dropping one displaced-coordinate hit; the review set the cohort's caveats |
| [chain_support](../benchmark/chain_support) | the whole-chain long-read support axis on and off on the five-caller chr22 set, with 1 M short-read pairs and without | With short reads: FLAIR recall 0.545 → 0.745 at precision 0.993 (146 genuine, 0 false lifted), IsoQuant 111 / 0, ESPRESSO 132 / 1, Bambu 83 / 15, TALON 148 / 54; without short reads the lift separates caller-level false chains (FLAIR, IsoQuant), not alignment-level ones |
| [caller_filters](../benchmark/caller_filters) | three post-caller filters on the same output (none, SQANTI3 rules, PanIsoGuard + BAM), five callers, no short reads; and each caller's own knob against PanIsoGuard | SQANTI3's rules filter adds no precision to 4 of 5 callers without short reads and costs 10–13 % of their genuine isoforms; PanIsoGuard raises precision on every caller (FLAIR 0.839 → 0.993 at recall 0.75, TALON 0.396 → 0.716 at 0.71) |
| [tusco_negatives](../benchmark/tusco_negatives) | the available real-data negatives: NIC/NNC isoforms IsoQuant reports in the 46 TUSCO single-isoform genes on LRGASP WTC11 and two HPRC Kinnex samples | 13 presumed-false novel isoforms: PanIsoGuard confirms 3 (HG02717, whole-chain reads at 3–5, no short reads), SQANTI3's rules filter keeps 8; with short reads 0 of 4 confirmed and 3 `ARTIFACT` |
| [sqanti_sim_heldout](../benchmark/sqanti_sim_heldout) | the chr22 protocol repeated unchanged on chr20 and chr19 with the chr22 thresholds | Without short reads PanIsoGuard + BAM AUPRC 0.962 / 0.957 vs SQANTI3 rules 0.946 / 0.928; with 1 M short-read pairs a tie with the one-line rule, as on chr22 |
| [lrgasp_targets](../benchmark/lrgasp_targets) | Every amplicon-tested LRGASP WTC11 target chain adjudicated as a would-be isoform, no caller in the loop | 26 unsupported chains: PanIsoGuard confirms 0 on PacBio and 3 on ONT (8–60 whole-chain reads), SQANTI3's rules filter keeps 18; 64 supported chains: 27 lifted on PacBio, 36 on ONT |
| [verdict_distribution](../benchmark/verdict_distribution) | descriptive: the verdict classes of the novel isoforms in the lrgasp_wtc11 and refbias_3way runs | Without short reads 37–66 % of a library's novel isoforms are lifted to probable and the rest held; with short reads the held ones become confirmed or artifact calls, except novel combinations of known junctions |
| [rts_sensitivity](../benchmark/rts_sensitivity) | what the ruleset-0.0.5 RT-switch policy (a whole-chain or consensus lift stands with `;rt_switch_sign`) moves in every tracked run | SQANTI-SIM 26 / 0 and 57 / 0 genuine; end2end 4 / 0; LRGASP keeps the 4 validated isoforms the hold would have lost |
| [lrgasp_wtc11](../benchmark/lrgasp_wtc11) | LRGASP's 89 WTC11 RT-PCR and amplicon-sequencing targets; our own IsoQuant / SQANTI3 / PanIsoGuard run on the public WTC11 PacBio cDNA, ONT cDNA and Illumina reads | of the supported NIC/NNC targets IsoQuant reconstructed and that are still novel in GENCODE v49 (17 PacBio, 16 ONT), PanIsoGuard discards none; SQANTI3's rules filter discards 2 and 4, all for `RTS_stage` |
| [refbias_3way](../benchmark/refbias_3way) | the rescue against an independent standard: the same reads through GRCh38 (A), an SNV-personalized diploid (B1), the HPRC phased assembly (B2), and PanIsoGuard post hoc on A (C) | at the rule level 31 of 32 and 37 of 39 SNV-class junctions on B1 agree with the personal-genome alignment (the misses are reads from the other haplotype); at the verdict level 12 of 32 and 21 of 39, because the `REFERENCE_BIAS` verdict needs every novel junction of the isoform explained; shifted-known (indel / SV) junctions are out of scope by construction |
| [sirv_multicaller](../benchmark/multicaller) | the multi-caller consensus on a second reference, Lexogen SIRV-Set4 | single-caller novel precision 0.037 → ≥ 2 callers 1.000 against the SIRV hidden-chain truth (chr22: 0.012 → ≥ 3 callers 0.980) |
| [sirv](../benchmark/sirv) | Lexogen SIRV-Set4 spike-in control | Specificity(false) 0.944, recall 1.000; the 3 false positives are FLAIR mis-collapses of individually real junctions, a documented limit of the short-read axis |
| [calibration](../benchmark/calibration) | class → empirical P(genuine), train/test split | controlled Brier 0.0000 / ECE 0.0013; end-to-end Brier 0.0121 / ECE 0.0114 |
| [multicaller](../benchmark/multicaller) | five callers (FLAIR, IsoQuant, Bambu, ESPRESSO, TALON) on one alignment, integrated by `combine`, scored against SQANTI-SIM truth | single-caller novel precision 0.012 (per-caller 0.40–0.98); the "≥ k callers" curve peaks at ≥ 3 (P 0.980, R 0.910, F1 0.944) |
| [wholegenome_multicaller](../benchmark/wholegenome_multicaller) | public GM12878 whole genome (IsoQuant, Bambu, ESPRESSO), no truth; canonical-motif check and decision impact | Small N: 254 single-caller vs 31 consensus of 285 novel chains; consensus novels 100 % canonical vs 97.6 % for single-caller ones |
| [merge_comparison](../benchmark/merge_comparison) | `combine` against gffcompare -i, TAMA and a 0–20 bp wobble sweep on the same five-caller chr22 set | `combine` groups like gffcompare -i (same n_callers distribution and PR curve); ≥ 3-caller precision similar under every matcher on this set (0.975–0.980; at ≥ 2 callers 0.836–0.933); 4 of 769 genuine novel chains (TAMA) or 1 (wobble ≥ 5 bp) cross the single/multi-caller boundary (recounted 2026-11-01) |
| [sqanti_sim](../benchmark/sqanti_sim) | SQANTI-SIM on GENCODE chr22: 769 transcripts deleted → PBSIM3 HiFi → FLAIR → SQANTI3 → `adjudicate` | precision 0.989, specificity 0.966, AUPRC 0.971 on the oracle junction set; the real baselines are in [sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h) |
| [pangenome](../benchmark/pangenome) | HPRC v1.1 Minigraph-Cactus chr22 graph → `vg deconstruct` → population-deletion junctions; the sample is not in the graph | 0 false rescues on real FLAIR novel junctions; the rescue fires when a novel intron is a population deletion (population provenance) and is held `UNRESOLVED` under circular-risk provenance |

**What the reference-bias numbers mean.** In `hg002`, `hg002_wholegenome`,
`giab_cohort_rescue`, `hg03516_refbias` and `refbias_cohort` a junction is labelled
`CREATED` when it is non-canonical on GRCh38 and canonical on the individual's haplotype,
which is the criterion the variant rescue applies. "N of N rescued, 0 false, N of N held"
therefore shows that the rule is implemented as specified and that the provenance gate
works; it is not independent evidence that these junctions are reference-bias artifacts.
The motif-creating base is intronic, so reads that use the junction never contain it. The
independent check is [refbias_phasing](../benchmark/refbias_phasing): for heterozygous
junctions, do the reads that use the junction carry, at the gene's other heterozygous SNVs,
the alleles of the haplotype that makes it canonical? In HG03516 and HG02717, 27 of 28
testable junctions pass; the 34 homozygous junctions cannot be phased. The yield is
30–45 junctions per divergent genome.

**Accuracy against real baselines.** On the SQANTI-SIM truth set,
[sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h) compares PanIsoGuard with the SQANTI3
rules and ML filters and with a one-line short-read rule, on oracle and on simulated
Illumina junctions. At ranking PanIsoGuard ties the one-line rule, with or without the BAM;
at the operating point it is more precise and, with the whole-chain read support of 0.0.7,
recalls 0.75 of the genuine isoforms with the BAM (0.55 without). With no short reads and
the BAM it ranks above SQANTI3's rules filter (AUPRC 0.953 vs 0.909); without either it
confirms nothing.

## Tracked result files

The number behind every protocol is committed under
[`../benchmark/results/`](../benchmark/results/) as a schema-checked `metrics.json`, indexed
by [`MANIFEST.tsv`](../benchmark/results/MANIFEST.tsv):

- `status=tracked`, self-contained: `controlled_truth` and `synthetic_axes` are regenerated
  and drift-checked by [`collect.py`](../benchmark/results/collect.py); the same drivers run
  in CI as the `integration_controlled_truth` and `integration_synthetic_axes` CTests.
- `status=tracked`, heavy pipeline: the other protocols carry the result of a run with the
  binary named in the envelope, produced by their documented `command`. The SQANTI-SIM
  threshold sweep is [`sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv)
  (30 configurations via [`sweep.py`](../benchmark/sqanti_sim/sweep.py)); the pangenome
  protocol's sensitivity and firewall half is the `integration_pangenome_gate1` CTest.
- `status=transcribed_pending_tracked_run` / `pending`: reserved for numbers still only
  transcribed from prose or not yet produced (none at present).

## What remains

| Item | Status | Notes |
|------|--------|-------|
| Pangenome rescue on a real graph | done for chr22 | HPRC v1.1 Minigraph-Cactus chr22 (`vg deconstruct` → population-deletion junctions → `--pangenome-junctions`); see [pangenome](../benchmark/pangenome). Not done: genome-wide scale, insertion- and inversion-based junctions, an in-process graph traversal. |
| LRGASP caller output ([lrgasp/](../benchmark/lrgasp)) | blocked | The pre-run LRGASP caller GTFs are Synapse-gated. The multi-caller axis is scored on five-caller output in `multicaller` instead, and the validated LRGASP targets are used in `lrgasp_wtc11` and `lrgasp_targets`. |
| HG002 long-read RNA | blocked | No public HG002/GM24385 long-read RNA-seq set was found (ENCODE, ENA). The variant axis is checked on HG002 variants and the full pipeline on `end2end` and `sirv`. |

**Thresholds.** `config/rules.default.toml` holds conservative defaults. A sweep on the
SQANTI-SIM chr22 truth set
([`benchmark/results/sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv),
[`sweep.py`](../benchmark/sqanti_sim/sweep.py)) finds AUPRC flat (0.970–0.971) across the
`sj_min_uniq_reads × require_canonical_motif × perc_A_degradation` grid, with the shipped
defaults within 1e-4 of the best cell; the only drop is the support cliff once
`sj_min_uniq_reads` exceeds the data's coverage (AUPRC 0.919). The defaults have not been
swept on other data sets. The confidence classes are labels, not a probability and not an
ordered scale: on SQANTI-SIM without the BAM, `UNRESOLVED` is more precise than
`UNCONFIRMED` at every short-read depth
([sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h), `per_class.py`).
