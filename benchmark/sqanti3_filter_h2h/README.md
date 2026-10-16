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
| none | PIG+bam | **0.992** | 0.670 | 0.800 | **0.947** |
| sr_1M | SQ3rules | 0.918 | 0.890 | 0.904 | 0.909 |
| sr_1M | SQ3strict | 0.978 | 0.840 | 0.903 | 0.955 |
| sr_1M | allSJ3 | 0.976 | **0.945** | **0.960** | 0.968 |
| sr_1M | PIG | **0.990** | 0.545 | 0.703 | 0.969 |
| sr_1M | PIG+bam | **0.992** | 0.710 | 0.827 | **0.972** |
| oracle | allSJ3 | 0.976 | **1.000** | **0.988** | 0.976 |
| oracle | PIG | 0.989 | 0.589 | 0.738 | 0.971 |
| oracle | PIG+bam | 0.991 | 0.747 | 0.852 | 0.974 |

250k and 4M pairs give the same picture (chr22 junction discovery saturates by ~1M); full
tables, per SQANTI category, in
[../results/sqanti3_filter_h2h/metrics.json](../results/sqanti3_filter_h2h/metrics.json).

AUPRC difference vs `allSJ3`, bootstrap 95 % CI (2000 resamples):

| | oracle | sr_250k | sr_1M | sr_4M |
|---|---|---|---|---|
| PIG | [−0.015, +0.005] | [−0.008, +0.013] | [−0.010, +0.010] | [−0.010, +0.010] |
| PIG+bam | [−0.011, +0.008] | [−0.003, +0.017] | [−0.006, +0.013] | [−0.006, +0.013] |

With no short reads, vs `SQ3rules`: PIG [−0.105, −0.066] (worse), PIG+bam **[+0.023, +0.055]**
(better, since the 0.0.7 whole-chain read support). Before 0.0.6 the BAM gain with short
reads was significant at 250k pairs ([+0.002, +0.021]); the 0.0.6 rule that a junction feeds
the BAM gate only with ≥ 3 spanning reads removed that.

### Risk–coverage (1 M pairs)

PanIsoGuard abstains (`UNRESOLVED`, `UNCONFIRMED`); the other methods decide every isoform.
Coverage = share of the 870 isoforms given a decisive class (`NOVEL_CONFIRMED`,
`NOVEL_PROBABLE`, `ARTIFACT`); selective risk = wrong decisive calls among them.

| method | coverage | selective risk |
|---|---:|---:|
| SQ3rules | 1.000 | 0.159 |
| SQ3strict | 1.000 | 0.151 |
| allSJ3 | 1.000 | 0.066 |
| PIG | 0.545 | 0.017 |
| PIG+bam | 0.707 | 0.021 |
| PIG+bam, no short reads | 0.602 | 0.029 |

Every setting's point is in the envelope (`coverage`, `selective_risk` per method).

## Reading

- **Ranking:** PanIsoGuard ties `allSJ3` in every setting, with or without the BAM
  (every AUPRC difference is within the bootstrap interval).
- **Operating point:** PanIsoGuard's precision is the highest (0.990–0.992 vs 0.976).
  Without the BAM it holds nearly all NIC isoforms `UNRESOLVED` (no individually-novel
  junction to corroborate) and all ISM isoforms `UNCONFIRMED`, so it recalls ~0.55 against
  `allSJ3`'s ~0.94. With the BAM (0.0.7), NIC isoforms whose whole intron chain is carried
  by ≥ 3 long reads become `NOVEL_PROBABLE`: recall 0.71, 120 genuine gained, 0 false
  ([../chain_support](../chain_support)). ISM stays held. On NNC it has the same TP/FP
  counts as `allSJ3`.
- **No short reads:** without the BAM PanIsoGuard makes no positive calls (everything is
  held `UNRESOLVED` unless `--caller-support` is given) and ranks below SQANTI3's default
  rules. With the BAM it keeps 489 genuine / 4 false (recall 0.67 at precision 0.992) and
  ranks above them; those calls rest on the caller's own reads, hence `NOVEL_PROBABLE`.
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
the BAM). The classes are not an ordered scale: without the BAM, `UNRESOLVED` (mostly NIC
isoforms with no novel junction to test) is more precise than `UNCONFIRMED` at every depth
(250 k, 1 M, 4 M pairs and the oracle give the same order).

| class | genuine | false | precision (no BAM) | genuine | false | precision (BAM) |
|---|---:|---:|---:|---:|---:|---:|
| NOVEL_CONFIRMED | 357 | 3 | 0.992 | 351 | 3 | 0.992 |
| NOVEL_PROBABLE | 41 | 1 | 0.976 | 167 | 1 | 0.994 |
| UNRESOLVED | 152 | 7 | 0.956 | 0 | 6 | 0.000 |
| UNCONFIRMED | 176 | 61 | 0.743 | 203 | 46 | 0.815 |
| ARTIFACT | 4 | 68 | 0.056 | 9 | 84 | 0.097 |

With the BAM (0.0.7), the genuine `UNRESOLVED` NIC isoforms move to `NOVEL_PROBABLE` through
whole-chain read support (120 of the 167), and `UNRESOLVED` keeps the 6 false ones.
`UNCONFIRMED` is two populations: ISM pass-through (154 genuine / 13 false, 0.92; the simulator
does not model intra-priming, so this is an upper bound) and novel isoforms with weak support
(no BAM: 22 / 48, 0.31; with BAM: 49 / 33, 0.60), below the 0.83 base rate. Numbers from the
0.0.7 binary (ruleset `builtin-0.0.4`); `per_class_v006.txt` of `run.sh`'s work dir has every
depth.
