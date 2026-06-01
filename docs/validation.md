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
| [hg002](../benchmark/hg002) | real GIAB HG002 v4.2.1 SNVs → personalized haplotype, `--haplotype-provenance wgs` | headline `PAN_REF_RESCUED` validated on **real, independent (WGS-derived) variants** — non-circular; **0 false rescues** |
| [sirv](../benchmark/sirv) | Lexogen SIRV-Set4 spike-in control, dense overlapping isoforms | specificity(false)=0.944, recall=1.000; the 3 FPs are FLAIR mis-collapses of *individually real* junctions — a documented short-read-axis limitation |
| [calibration](../benchmark/calibration) | class → empirical P(genuine), train/test split | **well-calibrated**: controlled Brier 0.0000 / ECE 0.0013; end-to-end Brier 0.0121 / ECE 0.0114 |
| [multicaller](../benchmark/multicaller) | real FLAIR + IsoQuant output integrated by `combine` | 397 + 157 → **426 unique intron chains**, 128 agreed by both callers (validates running-ID integration on genuine multi-caller output) |

Together these cover all four implemented axes (short-read, BAM mapping, variant) and
the full pipeline (real FLAIR/SQANTI3), plus the headline reference-bias rescue on
real human variation, calibration, and multi-caller integration.

## What remains

| Item | Status | Why |
|------|--------|-----|
| **SQANTI-SIM AUPRC sweep** ([sqanti_sim/](../benchmark/sqanti_sim)) | not run | The canonical simulator would add a per-class precision/recall + **AUPRC / ΔAUPRC** threshold sweep (letting `benchmark`/`ablate` report AUPRC, not only contingency counts). Equivalent truth-based P/R is already shown by `controlled_truth` + `end2end`; this adds the standard-tool curve and a calibration anchor for the default thresholds. |
| **Pangenome real-graph rescue (GATE-1)** | not run | The variant axis is validated on real HG002 variants; the **file-based pangenome tier** still needs HPRC v1.1 + `vg`/`rpvg` to extract graph-supported junctions for an out-of-graph individual. Heavy external tooling; the axis is marked *experimental*. |
| **LRGASP real data** ([lrgasp/](../benchmark/lrgasp)) | blocked | The LRGASP pre-run caller GTFs are **Synapse-gated**. The multi-caller capability is already validated on genuine FLAIR+IsoQuant output (see `multicaller`). |
| **HG002 long-read RNA** | blocked | No clean public HG002/GM24385 long-read RNA-seq dataset was found (ENCODE/ENA empty). The variant axis is instead validated on real HG002 *variants* + the full pipeline on `end2end`/`sirv`. |

**Calibration caveat.** `config/rules.default.toml` thresholds are conservative
defaults. They are well-calibrated on the run truth sets above (low ECE/Brier), but
have **not** been tuned against a full SQANTI-SIM AUPRC sweep; treat the confidence
classes as calibrated *ordinal* evidence integration, not a tuned probability.
