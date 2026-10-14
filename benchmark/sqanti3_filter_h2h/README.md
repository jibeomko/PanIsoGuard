# Head-to-head vs the SQANTI3 filter (and a one-line short-read rule)

**Question:** on the same truth set and the same inputs, is PanIsoGuard more accurate than
what users already run — the SQANTI3 rules / ML filter — or than the simplest possible
short-read rule? **Answer: no.** Its ranking ties a one-line rule; at its default
operating point it trades a lot of recall for a little precision; with no short reads it
is worse than the SQANTI3 rules filter. The only measurable edge is small and comes from
the BAM mapping axis.

## Setup

Truth set = [../sqanti_sim](../sqanti_sim) (GENCODE v49 chr22, SQANTI-SIM, FLAIR, SQANTI3
6.0.1), scored with the same exact-intron-chain labels. Task: genuine vs false among the
**870 non-FSM multi-exon isoforms** (730 genuine / 140 false).

Short-read settings (every method gets the same junctions):

| setting | junctions |
|---|---|
| `none` | no short reads |
| `oracle` | the truth-derived SJ.tab of `../sqanti_sim` (every simulated junction, 50 reads) |
| `sr_N` | **realistic**: N wgsim 2×100 pairs from the simulated transcripts, weighted by their long-read counts → STAR 2-pass → `SJ.out.tab` |

Methods:

| name | what |
|---|---|
| `SQ3rules` | `sqanti3_filter.py rules`, default JSON |
| `SQ3strict` | `sqanti3_filter.py rules` with `min_cov ≥ 3` required |
| `SQ3ml` | `sqanti3_filter.py ml`, defaults |
| `allSJ3` | one line: keep an isoform iff **every junction has ≥ 3 unique short reads** (SQANTI3 `min_cov`) |
| `PIG` / `PIG+bam` | `panisoguard adjudicate` (± `--bam`); positive = `HIGH/NOVEL_PROBABLE`, ranked by class |

## Result (non-FSM, n = 870)

| short reads | method | precision | recall | F1 | AUPRC |
|---|---|---:|---:|---:|---:|
| none | SQ3rules | 0.918 | 0.886 | **0.902** | **0.909** |
| none | PIG | — (0 calls) | 0.000 | — | 0.824 |
| none | PIG+bam | — (0 calls) | 0.000 | — | 0.844 |
| sr_1M | SQ3rules | 0.918 | 0.890 | 0.904 | 0.909 |
| sr_1M | SQ3strict | 0.978 | 0.840 | 0.903 | 0.955 |
| sr_1M | allSJ3 | 0.976 | **0.945** | **0.960** | 0.968 |
| sr_1M | PIG | **0.990** | 0.545 | 0.703 | 0.969 |
| sr_1M | PIG+bam | **0.990** | 0.545 | 0.703 | **0.973** |
| oracle | allSJ3 | 0.976 | **1.000** | **0.988** | 0.976 |
| oracle | PIG | 0.989 | 0.589 | 0.738 | 0.971 |
| oracle | PIG+bam | 0.989 | 0.589 | 0.738 | 0.975 |

250k and 4M pairs give the same picture (chr22 junction discovery saturates by ~1M); full
tables, per SQANTI category, in
[../results/sqanti3_filter_h2h/metrics.json](../results/sqanti3_filter_h2h/metrics.json).

AUPRC difference vs `allSJ3`, bootstrap 95 % CI (2000 resamples):

| | oracle | sr_250k | sr_1M | sr_4M |
|---|---|---|---|---|
| PIG | [−0.015, +0.005] | [−0.008, +0.013] | [−0.010, +0.010] | [−0.010, +0.010] |
| PIG+bam | [−0.011, +0.009] | [−0.004, +0.017] | [−0.006, +0.014] | [−0.006, +0.014] |

With no short reads, vs `SQ3rules`: PIG [−0.105, −0.066], PIG+bam [−0.089, −0.042] (worse).
Before 0.0.6 the BAM gain was significant at 250k pairs ([+0.002, +0.021]); the 0.0.6 rule
that a junction feeds the BAM gate only with ≥ 3 spanning reads removed that.

## Reading

- **Ranking:** without the BAM, PanIsoGuard ties `allSJ3` in every setting. With the BAM,
  the indel-near mapping signal adds ~+0.004 AUPRC, within the bootstrap interval at
  every depth.
- **Operating point:** PanIsoGuard's precision is the highest (0.990 vs 0.976), but
  it holds nearly all NIC isoforms `UNRESOLVED` (no individually-novel junction to
  corroborate) and all ISM isoforms `UNCONFIRMED`. So it recalls ~0.55 against
  `allSJ3`'s ~0.94. On NNC it has the same TP/FP counts as `allSJ3`.
- **No short reads:** PanIsoGuard makes no positive calls at all (everything is held
  `UNRESOLVED` unless `--caller-support` is given), and its ranking is below SQANTI3's
  default rules.
- `SQ3ml` is not used for conclusions: its default training labels (FSM = TP, NNC = TN)
  are wrong for this simulation, where most NNC are genuine.

So PanIsoGuard is **not a more accurate filter**. What it adds is a per-call verdict
with the rule trace and the mechanism behind it, and explicit abstention.

## Fixed along the way

The first run of this comparison exposed a projection bug: `PARTIAL` short-read support
(some but not all novel junctions corroborated) was promoted to `NOVEL_PROBABLE`.
Those were exactly PanIsoGuard's extra NNC false positives (oracle run: 6 vs 3 for the
short-read rule; 3 vs 3 after the fix). Since 0.0.4, `PARTIAL` → `UNCONFIRMED`.

## Caveats

- Single simulated truth set (chr22, PBSIM3 HiFi, FLAIR). Simulated short reads are clean
  (no intron retention, pre-mRNA, or mapping-hard regions), so `allSJ3` may look better
  here than on real data.
- The BAM axis thresholds are HiFi-calibrated and over-fire on ONT
  ([../gm12878_realdata](../gm12878_realdata)).

## Reproduce

```bash
bash ../sqanti_sim/run.sh     # the truth set (once)
bash run.sh                   # edit the paths at the top first
```

## Per-class precision (`per_class.py`)

The same truth, split by PanIsoGuard class (1 M simulated pairs; `pig` = no BAM, `pigbam` = with
the BAM). The classes are not an ordered scale: `UNRESOLVED` (mostly NIC isoforms with no novel
junction to test) is more precise than `UNCONFIRMED` at every depth (250 k, 1 M, 4 M pairs and
the oracle give the same order).

| class | genuine | false | precision (no BAM) | genuine | false | precision (BAM) |
|---|---:|---:|---:|---:|---:|---:|
| NOVEL_CONFIRMED | 357 | 3 | 0.992 | 351 | 3 | 0.992 |
| NOVEL_PROBABLE | 41 | 1 | 0.976 | 47 | 1 | 0.979 |
| UNRESOLVED | 152 | 7 | 0.956 | 151 | 6 | 0.962 |
| UNCONFIRMED | 176 | 61 | 0.743 | 172 | 46 | 0.789 |
| ARTIFACT | 4 | 68 | 0.056 | 9 | 84 | 0.097 |

`UNCONFIRMED` is two populations: ISM pass-through (154 genuine / 13 false, 0.92; the simulator
does not model intra-priming, so this is an upper bound) and novel isoforms with weak support
(no BAM: 22 / 48, 0.31; with BAM: 18 / 33, 0.35), below the 0.83 base rate. Numbers from the
0.0.6 binary (ruleset `builtin-0.0.3`); `per_class_v006.txt` of `run.sh`'s work dir has every
depth.
