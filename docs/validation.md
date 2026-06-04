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
| [end2end](../benchmark/end2end) | PBSIM3 → minimap2 → FLAIR3 → SQANTI3 → `adjudicate`, GENCODE chr22 | on the SQANTI-novel set: **specificity(false)=0.996**, genuine-recall(decisive)=1.000 — correct on *real* caller + SQANTI3 output |
| [hg002](../benchmark/hg002) | real GIAB HG002 v4.2.1 SNVs → personalized haplotype, `--haplotype-provenance wgs` | headline `PAN_REF_RESCUED_FALSE_NOVEL` validated on **real, independent (WGS-derived) variants** — non-circular; **0 false rescues** |
| [hg002_wholegenome](../benchmark/hg002) (yield) | the *same* scan extended to **whole-genome** HG002 (chr1–22, 597,781 GENCODE introns) | **rescue YIELD quantified**: reference bias at splice junctions is rare but real — **33** genuine reference-bias junctions genome-wide, **33/33 rescued** (100 % sensitivity), **0 false**, **33/33 held** by the firewall under circular-risk provenance. The rescue is a **high-specificity/sensitivity guardrail with a low base rate** (~1 in 18k introns), not a high-yield discovery engine |
| [sirv](../benchmark/sirv) | Lexogen SIRV-Set4 spike-in control, dense overlapping isoforms | specificity(false)=0.944, recall=1.000; the 3 FPs are FLAIR mis-collapses of *individually real* junctions — a documented short-read-axis limitation |
| [calibration](../benchmark/calibration) | class → empirical P(genuine), train/test split | **well-calibrated**: controlled Brier 0.0000 / ECE 0.0013; end-to-end Brier 0.0121 / ECE 0.0114 |
| [multicaller](../benchmark/multicaller) | **five** callers (FLAIR + IsoQuant + Bambu + ESPRESSO + TALON) on one shared alignment, integrated by `combine`, scored vs SQANTI-SIM truth | single-caller novel precision **0.012** (per-caller precision 0.40–0.98); PR curve over "≥ k callers" peaks at **≥3 (P 0.980, R 0.910, F1 0.944)** — the optimal consensus threshold scales with caller count. Engine reproduces it: `adjudicate --caller-support` (long-read only) takes FLAIR from **0 confident novels (all AMBIGUOUS) → 0.975 precision** by promoting cross-caller-agreed novels |
| [wholegenome_multicaller](../benchmark/wholegenome_multicaller) | **real** public GM12878 whole-genome (IsoQuant + Bambu + ESPRESSO), no truth → orthogonal **canonical-motif** validation | *Honest, small-N.* Disagreement holds at genome scale (**254 single-caller vs 31 consensus** of 285 novel chains). Consensus novels **100 %** canonical vs single-caller **97.6 %** — direction correct but small lift (+0.024, ceiling effect: stringent production callers are already clean). On real data the value is the reproducible high-confidence core, not bulk artifact removal |
| [merge_comparison](../benchmark/merge_comparison) | head-to-head vs the incumbent merge tools on the **same** 5-caller chr22 set: `combine` vs **gffcompare -i**, **TAMA**, and a 0–20 bp wobble sweep | **`combine` ≡ gffcompare -i exactly** (same n_callers distribution + PR curve) — a clean re-implementation, not a novel merge. The consensus ≥3 precision is **matcher-robust** (0.975–0.980 across exact / gffcompare / TAMA / wobble), and **0** genuine novels cross the single↔multi boundary under any fuzzy matcher → the consensus gate is a property of the **data**, not of exact matching. Honest: consensus is established practice; the **reference-bias rescue** is the differentiator |
| [sqanti_sim](../benchmark/sqanti_sim) | canonical SQANTI-SIM (GENCODE chr22): delete 769 transcripts → PBSIM3 HiFi → FLAIR → SQANTI3 → `adjudicate` | genuine-novel detection **precision 0.982, specificity 0.946, AUPRC 0.970** (baseline 0.831); NNC recall 1.000; known→HIGH_CONF_KNOWN 1.000; **0 genuine→ARTIFACT**. Moderate overall recall = honest abstention on NIC-combinatorial / ISM-partial, not misclassification |
| [pangenome](../benchmark/pangenome) (GATE-1) | real **HPRC v1.1** MC GRCh38 chr22 graph → `vg deconstruct` → population-deletion junctions; sample **out-of-graph** (non-circular) | **0 false rescues** on real FLAIR novel junctions; rescue **fires** when a novel intron is a real population deletion (population provenance); circularity firewall **holds** all rescues `AMBIGUOUS` under circular-risk provenance |

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
  `end2end` reproduces spec 0.996 on a seeded pbsim run; `hg002` reproduces the
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
the SQANTI-SIM AUPRC (0.970) confirms strong genuine-vs-false ranking. A threshold
**sweep** on the SQANTI-SIM v49 chr22 truth set
([`benchmark/results/sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv),
produced by [`benchmark/sqanti_sim/sweep.py`](../benchmark/sqanti_sim/sweep.py)) shows
AUPRC is **robust (0.969–0.970)** across the `sj_min_uniq_reads × require_canonical_motif
× perc_A_degradation` grid and the shipped `default-0.0.1` config is within 1e-4 of the
grid-best — the only degradation is the expected support cliff once `sj_min_uniq_reads`
exceeds the data's coverage (AUPRC → 0.919). So the conservative defaults are
**near-optimal there**; they are not yet swept on additional datasets. Treat the
confidence classes as calibrated *ordinal* evidence integration, not a tuned probability.
