# Post-caller filters: PanIsoGuard next to what the caller and SQANTI3 already offer

Two questions a reader of the head-to-head asks next. If I already have a caller's output, what do
I gain by running PanIsoGuard on it instead of SQANTI3's rules filter, or nothing? And if my caller
has its own precision knob, is tightening it not the same thing? Both on the chr22 SQANTI-SIM
truth set, five callers on one alignment ([../multicaller](../multicaller)), no short reads unless
said (the common case for a long-read study), labels as in
[../sqanti3_filter_h2h](../sqanti3_filter_h2h). `run.sh` + `score_filters.py`; envelope
[../results/caller_filters/metrics.json](../results/caller_filters/metrics.json).

## Three post-caller filters on the same output

Recall is the share of the caller's own genuine novel isoforms that survive the filter; the
caller's recall against the truth is outside this comparison. SQANTI3's rules filter has no
coverage column here, so it is the perc_A + RTS rules.

| caller (novel isoforms, genuine) | none, P | SQANTI3 rules, P / R | PanIsoGuard + BAM, P / R (coverage / risk) | + 1 M short-read pairs, P / R |
|---|---:|---:|---:|---:|
| FLAIR (870, 730) | 0.839 | 0.918 / 0.886 | **0.993 / 0.748** (0.63 / 0.007) | 0.993 / 0.745 |
| IsoQuant (618, 607) | 0.982 | 0.982 / 0.886 | 0.991 / 0.708 (0.70 / 0.009) | 0.991 / 0.708 |
| Bambu NDR 0.5 (552, 409) | 0.741 | 0.734 / 0.875 | 0.764 / 0.951 (0.92 / 0.236) | 0.748 / 0.951 |
| ESPRESSO (573, 555) | 0.969 | 0.966 / 0.868 | 0.987 / 0.939 (0.92 / 0.013) | 0.987 / 0.935 |
| TALON, filtered (1,760, 697) | 0.396 | 0.406 / 0.897 | **0.716 / 0.709** (0.39 / 0.284) | 0.709 / 0.697 |

- **Without short reads SQANTI3's rules filter does not add precision** to four of the five
  callers (± 0.01) and costs 10–13 % of their genuine isoforms; FLAIR is the exception (+0.08).
- **PanIsoGuard raises precision on every caller**: by 0.01–0.02 where the caller is already
  strict (IsoQuant, ESPRESSO, Bambu), by 0.15 on FLAIR and 0.32 on TALON, keeping 71–95 % of the
  genuine isoforms. The selective risk says where the gain is cheap: 0.7–1.3 % of its decisive
  calls are wrong on FLAIR, IsoQuant and ESPRESSO, 24–28 % on Bambu and TALON, whose false
  chains are carried by reads ([../chain_support](../chain_support)).
- **For a strict caller the gain is reasons, not precision.** On IsoQuant PanIsoGuard moves
  precision from 0.982 to 0.991 and holds 29 % of the genuine isoforms for review; what it adds
  there is the per-isoform evidence, the reference-bias check and the verdict trail.

## The caller's own knob

Recall here is relative to the genuine novel chains of the permissive setting (Bambu NDR 0.5:
409; TALON before its filter: 761), so a stricter setting and PanIsoGuard on the permissive
output are on one axis.

| setting | novel isoforms | as is, P / R | SQANTI3 rules | PanIsoGuard + BAM on this output |
|---|---:|---:|---:|---:|
| Bambu NDR 0.5 (tracked) | 552 | 0.741 / 1.000 | 0.734 / 0.875 | 0.764 / 0.951 |
| Bambu NDR 0.3 | 288 | 0.819 / 0.577 | | |
| Bambu NDR 0.2 | 167 | 0.904 / 0.369 | | |
| Bambu NDR 0.1 | 71 | 0.986 / 0.171 | | |
| TALON before `talon_filter_transcripts` | 6,904 | 0.110 / 1.000 | 0.117 / 0.886 | 0.353 / 0.712 |
| TALON after it (tracked) | 1,760 | 0.396 / 0.916 | 0.406 / 0.821 | **0.716 / 0.649** |

- **Bambu's NDR buys precision faster than PanIsoGuard does, at a steep recall price.** NDR 0.3
  reaches 0.82 by discarding 42 % of the genuine chains; PanIsoGuard on the NDR 0.5 output keeps
  95 % of them at 0.76. They are different trades, not the same filter: the NDR ranks by
  Bambu's own model, PanIsoGuard by read-level evidence and leaves the discarded ones labelled.
- **PanIsoGuard is not a substitute for a caller's read-count filter.** On TALON's unfiltered
  output it reaches 0.353 at recall 0.71, below TALON's own filter (0.396 at 0.92); on top of
  that filter it reaches 0.716. The two stack: the caller's filter removes one-read noise, the
  adjudicator removes what survived it with a reason attached.

## A supervised ceiling: SQANTI3's ML filter trained on the truth

The SQ3ml rows of the head-to-head use SQANTI3's default labels (FSM = true, NNC = false), which
are wrong for a simulation where most NNC are genuine. `../sqanti3_filter_h2h/ml_truth_cv.py`
trains the same random forest on the simulated truth itself (genuine = TP, false = TN) in
two-fold cross-validation and scores the held-out fold. The filter refuses fewer than 250
isoforms per class, which rules out every caller but TALON (697 genuine / 1,063 false):

| TALON, no short reads | precision | recall | AUPRC |
|---|---:|---:|---:|
| SQANTI3 ML, default labels | 0.40 | 0.90 | – |
| SQANTI3 ML, trained on the truth (2-fold CV) | 0.740 | 0.518 | 0.737 |
| PanIsoGuard + BAM, no training | 0.716 | 0.709 | 0.629 |

Given ground-truth labels, SQANTI3's own features rank TALON's isoforms better than
PanIsoGuard's rules do (AUPRC 0.737 vs 0.629) and reach a similar precision at lower recall.
That is the honest ceiling of a supervised filter on a permissive caller; it needs labelled
truth, which no real data set provides, and it says nothing about why an isoform was kept.

## Reproduce

```bash
WORK=/path/to/sqanti_sim/work MC=/path/to/multicaller CHAIN=/path/to/chain_support_out \
H2H=/path/to/h2h_work GENOME=/path/to/chr22.fa OUT=/tmp/caller_filters bash run.sh
```

Needs the SQANTI3 env (`sqanti3_qc.py`, `sqanti3_filter.py`), the `bambu` R env for the NDR
sweep, and `scikit-learn` for the scorer's imports. Public/simulated data only.
