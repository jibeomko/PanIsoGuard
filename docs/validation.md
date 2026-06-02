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
| [sirv](../benchmark/sirv) | Lexogen SIRV-Set4 spike-in control, dense overlapping isoforms | specificity(false)=0.944, recall=1.000; the 3 FPs are FLAIR mis-collapses of *individually real* junctions — a documented short-read-axis limitation |
| [calibration](../benchmark/calibration) | class → empirical P(genuine), train/test split | **well-calibrated**: controlled Brier 0.0000 / ECE 0.0013; end-to-end Brier 0.0121 / ECE 0.0114 |
| [multicaller](../benchmark/multicaller) | real FLAIR + IsoQuant output integrated by `combine` | 397 + 157 → **426 unique intron chains**, 128 agreed by both callers (validates running-ID integration on genuine multi-caller output) |
| [sqanti_sim](../benchmark/sqanti_sim) | canonical SQANTI-SIM (GENCODE chr22): delete 769 transcripts → PBSIM3 HiFi → FLAIR → SQANTI3 → `adjudicate` | genuine-novel detection **precision 0.982, specificity 0.946, AUPRC 0.970** (baseline 0.831); NNC recall 1.000; known→HIGH_CONF_KNOWN 1.000; **0 genuine→ARTIFACT**. Moderate overall recall = honest abstention on NIC-combinatorial / ISM-partial, not misclassification |

Together these cover the validated short-read, BAM mapping, and variant axes,
the full pipeline (real FLAIR/SQANTI3), the headline reference-bias rescue on
real human variation, calibration, and multi-caller integration. The file-based
pangenome axis has unit-level coverage and remains marked experimental until GATE-1.

## Tracked result artifacts

These numbers are no longer prose-only. The metric for every protocol is committed
under [`../benchmark/results/`](../benchmark/results/) as a schema-checked
`metrics.json`, indexed by [`MANIFEST.tsv`](../benchmark/results/MANIFEST.tsv):

- **`status=tracked`** — the two self-contained protocols (`controlled_truth`,
  `synthetic_axes`) are regenerated and drift-checked by
  [`collect.py`](../benchmark/results/collect.py), and the same drivers run in CI as
  the `integration_controlled_truth` / `integration_synthetic_axes` CTests, so a
  metric and its pass/fail assertion share one code path.
- **`status=transcribed_pending_tracked_run`** — the heavy protocols (e.g.
  `sqanti_sim`) carry the published number plus a `command` to reproduce it, pending a
  committed tracked run. The SQANTI-SIM **threshold sweep** that turns the conservative
  defaults into a tuned operating point is templated in
  [`sqanti_sim/sweep.tsv`](../benchmark/results/sqanti_sim/sweep.tsv) (the release
  blocker below).

## What remains

| Item | Status | Why |
|------|--------|-----|
| **Pangenome real-graph rescue (GATE-1)** | in progress | The variant axis is validated on real HG002 variants; the **file-based pangenome tier** is being validated against the HPRC v1.1 Minigraph-Cactus graph (`vg deconstruct` of the GRCh38 chr22 path → graph-supported junctions → `--pangenome-junctions`). Heavy external tooling; the axis is marked *experimental*. |
| **LRGASP real data** ([lrgasp/](../benchmark/lrgasp)) | blocked | The LRGASP pre-run caller GTFs are **Synapse-gated**. The multi-caller capability is already validated on genuine FLAIR+IsoQuant output (see `multicaller`). |
| **HG002 long-read RNA** | blocked | No clean public HG002/GM24385 long-read RNA-seq dataset was found (ENCODE/ENA empty). The variant axis is instead validated on real HG002 *variants* + the full pipeline on `end2end`/`sirv`. |

**Calibration caveat.** `config/rules.default.toml` thresholds are conservative
defaults. They are well-calibrated on the run truth sets above (low ECE/Brier) and
the SQANTI-SIM AUPRC (0.970) confirms strong genuine-vs-false ranking, but the
thresholds have **not** been swept/tuned per dataset; treat the confidence
classes as calibrated *ordinal* evidence integration, not a tuned probability.
