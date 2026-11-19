# Head-to-head vs the SQANTI3 filter (and a one-line short-read rule)

Question: on the same truth set and the same inputs, is PanIsoGuard more accurate than
what users already run (the SQANTI3 rules / ML filter) or than the simplest possible
rule? Answer: no. With short reads its ranking ties a one-line short-read rule; without them
it ranks above SQANTI3's default rules only when it is given the BAM, and then level with the
caller's own read count. At its default operating point it trades recall for precision.

## Setup

Truth set = [../sqanti_sim](../sqanti_sim) (GENCODE v49 chr22, SQANTI-SIM, FLAIR, SQANTI3
6.0.1), scored with the same exact-intron-chain labels. Task: genuine vs false among the
870 non-FSM multi-exon isoforms (730 genuine / 140 false).

Short-read settings (every method gets the same junctions):

| setting | junctions |
|---|---|
| `none` | No short reads |
| `oracle` | the truth-derived SJ.tab of `../sqanti_sim` (every simulated junction, 50 reads) |
| `sr_N` | realistic: N wgsim 2×100 pairs from the simulated transcripts, weighted by their long-read counts, STAR 2-pass, `SJ.out.tab` |

Methods:

| name | what |
|---|---|
| `SQ3rules` | `sqanti3_filter.py rules`, default JSON |
| `SQ3strict` | `sqanti3_filter.py rules` with `min_cov ≥ 3` required |
| `SQ3ml` | `sqanti3_filter.py ml`, defaults |
| `allSJ3` | One line: keep an isoform iff every junction has ≥ 3 unique short reads (SQANTI3 `min_cov`) |
| caller's read count, `SQ3rules` with `FL ≥ k` | no short reads only: [below](#is-this-a-read-count-filter-no-short-reads) |
| `PIG` / `PIG+bam` | `panisoguard adjudicate` (± `--bam`); positive = `NOVEL_CONFIRMED` / `NOVEL_PROBABLE`, ranked by class |

## Result (non-FSM, n = 870)

| short reads | method | precision | recall | F1 | AUPRC |
|---|---|---:|---:|---:|---:|
| none | SQ3rules | 0.918 | 0.886 | 0.902 | 0.909 |
| none | PIG | n/a (0 calls) | 0.000 | n/a | 0.824 |
| none | PIG+bam | 0.993 | 0.748 | 0.853 | 0.953 |
| sr_1M | SQ3rules | 0.918 | 0.890 | 0.904 | 0.909 |
| sr_1M | SQ3strict | 0.978 | 0.840 | 0.903 | 0.955 |
| sr_1M | allSJ3 | 0.976 | 0.945 | 0.960 | 0.968 |
| sr_1M | PIG | 0.990 | 0.545 | 0.703 | 0.970 |
| sr_1M | PIG+bam | 0.993 | 0.752 | 0.856 | 0.975 |
| oracle | allSJ3 | 0.976 | 1.000 | 0.988 | 0.976 |
| oracle | PIG | 0.989 | 0.589 | 0.738 | 0.971 |
| oracle | PIG+bam | 0.991 | 0.782 | 0.874 | 0.977 |

250k and 4M pairs give the same picture (chr22 junction discovery saturates by ~1M); full
tables, per SQANTI category, in
[../results/sqanti3_filter_h2h/metrics.json](../results/sqanti3_filter_h2h/metrics.json).

AUPRC difference vs `allSJ3`, bootstrap 95 % CI (2000 resamples):

| | oracle | sr_250k | sr_1M | sr_4M |
|---|---|---|---|---|
| PIG | [-0.015, +0.005] | [-0.008, +0.014] | [-0.010, +0.011] | [-0.010, +0.011] |
| PIG+bam | [-0.009, +0.011] | [-0.001, +0.020] | [-0.003, +0.017] | [-0.004, +0.017] |

With no short reads, vs `SQ3rules`: PIG [-0.105, -0.066] (worse), PIG+bam [+0.029, +0.061]
(better, since the 0.0.7 whole-chain read support). Before 0.0.6 the BAM gain with short
reads was significant at 250k pairs ([+0.002, +0.021]); the 0.0.6 rule that a junction feeds
the BAM gate only with ≥ 3 spanning reads removed that.

### Risk-coverage (1 M pairs)

PanIsoGuard abstains (`UNRESOLVED`, `UNCONFIRMED`); the other methods decide every isoform.
Coverage = share of the 870 isoforms given a decisive class (`NOVEL_CONFIRMED`,
`NOVEL_PROBABLE`, `ARTIFACT`); selective risk = wrong decisive calls among them.

| method | coverage | selective risk |
|---|---:|---:|
| SQ3rules | 1.000 | 0.159 |
| SQ3strict | 1.000 | 0.151 |
| allSJ3 | 1.000 | 0.066 |
| PIG | 0.544 | 0.015 |
| PIG+bam | 0.737 | 0.014 |
| PIG+bam, no short reads | 0.632 | 0.007 |

Every setting's point is in the envelope (`coverage`, `selective_risk` per method).

## Is this a read-count filter? (no short reads)

Without short reads PanIsoGuard's positive calls rest on reads that carry the whole intron
chain. A caller already reports a read count per isoform, and SQANTI3's rules filter can
require one (the `FL` descriptor), so the question is what the whole-chain evidence adds to a
plain count. `count_baseline.py` scores, on the same truth set: the caller's own count as a
filter at several thresholds and as a ranking; SQANTI3's rules filter run with its default
rules plus `FL ≥ k` in every rule (the `FL` column filled from the caller's counts, as
`sqanti3_qc.py --fl_count` fills it); the whole-chain read count with no other rule; and
PanIsoGuard with the BAM. The AUPRC of a count row ranks by the count itself; a rules-filter
row has one operating point. Envelope block `read_count_baseline`.

| no short reads | FLAIR: genuine / false kept | precision | recall | AUPRC | IsoQuant: genuine / false kept | precision | recall | AUPRC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SQANTI3 rules, default | 647 / 58 | 0.918 | 0.886 | 0.909 | 538 / 10 | 0.982 | 0.886 | 0.982 |
| caller's read count ≥ 3 | 730 / 140 | 0.839 | 1.000 | 0.958 | 607 / 10 | 0.984 | 1.000 | 0.992 |
| caller's read count ≥ 5 | 730 / 108 | 0.871 | 1.000 |  | 607 / 7 | 0.989 | 1.000 |  |
| caller's read count ≥ 10 | 729 / 68 | 0.915 | 0.999 |  | 605 / 7 | 0.989 | 0.997 |  |
| caller's read count ≥ 20 | 726 / 43 | 0.944 | 0.995 |  | 604 / 6 | 0.990 | 0.995 |  |
| SQANTI3 rules with FL ≥ 3 | 647 / 58 | 0.918 | 0.886 | 0.909 | 538 / 9 | 0.984 | 0.886 | 0.983 |
| SQANTI3 rules with FL ≥ 10 | 646 / 34 | 0.950 | 0.885 | 0.937 | 536 / 7 | 0.987 | 0.883 | 0.987 |
| whole-chain reads ≥ 3, nothing else | 730 / 41 | 0.947 | 1.000 | 0.964 | 607 / 8 | 0.987 | 1.000 | 0.988 |
| PanIsoGuard + BAM | 546 / 4 | 0.993 | 0.748 | 0.953 | 430 / 4 | 0.991 | 0.708 | 0.988 |

- **As a ranking it is level with the count.** The caller's read count alone ranks FLAIR's
  isoforms above SQANTI3's default rules (AUPRC 0.958 against 0.909, difference [+0.028, +0.070]),
  and PanIsoGuard does not rank them better than the count does (0.953, difference
  [-0.023, +0.015]; IsoQuant 0.988 against 0.992, [-0.013, +0.007]). "Above SQANTI3's
  rules filter" in the table above is a statement about the default rules, not about what a
  count can do.
- **At the same threshold the two counts are different evidence.** FLAIR reports no isoform
  with fewer than three supporting reads, so a count of three removes nothing: all 140 false
  isoforms pass. 41 of them have three reads that carry the whole chain; every one of the
  730 genuine isoforms does. A caller assigns a read to the isoform it is compatible with; a
  whole-chain read has every intron of the chain at exact coordinates in one alignment. The
  count needs a threshold near 20 to come close (precision 0.944 against 0.947), and that
  number moves with depth.
- **The verdict is a stricter operating point than any count threshold.** PanIsoGuard keeps
  546 genuine and 4 false isoforms (precision 0.993 at recall 0.748). Raising the caller's count
  until recall falls to that level keeps 518 and 19 (count ≥ 64, precision 0.965); the
  whole-chain count alone, 537 and 14 (precision 0.975). The difference comes from the rules around the
  count: an artifact sign holds the call, and an incomplete splice match is a known fragment,
  not a novel call. On IsoQuant, whose output has 11 false isoforms, every method is within
  about 0.01 of the others.
- **SQANTI3's rules with a count are a sound filter.** With `FL ≥ 10` they reach precision
  0.950 at recall 0.885 on FLAIR. A user who wants one keep-or-discard list without short reads
  and tunes that threshold on their own data gets most of the precision. What PanIsoGuard adds
  to that is not accuracy: the evidence per junction, the held class, and the genome.

## Reading

- Ranking: PanIsoGuard ties `allSJ3` in every setting, with or without the BAM
  (every AUPRC difference is within the bootstrap interval).
- Operating point: PanIsoGuard's precision is the highest (0.990-0.992 vs 0.976).
  Without the BAM it holds nearly all NIC isoforms `UNRESOLVED` (no individually-novel
  junction to corroborate) and all ISM isoforms `UNCONFIRMED`, so it recalls ~0.55 against
  `allSJ3`'s ~0.94. With the BAM (0.0.7), NIC isoforms whose whole intron chain is carried
  by ≥ 3 long reads become `NOVEL_PROBABLE`: recall 0.75, 151 genuine gained, 0 false
  ([../chain_support](../chain_support)). ISM carried end to end becomes `KNOWN_FRAGMENT`,
  which is not a novel call. On NNC it has the same TP/FP
  counts as `allSJ3`.
- No short reads: without the BAM PanIsoGuard makes no positive calls (everything is
  held `UNRESOLVED` unless `--caller-support` is given) and ranks below SQANTI3's default
  rules. With the BAM it keeps 546 genuine / 4 false (recall 0.75 at precision 0.993) and
  ranks above them, level with the caller's read count (the section above); those calls rest
  on the caller's own reads, hence `NOVEL_PROBABLE`.
- `SQ3ml` is not used for conclusions: its default training labels (FSM = TP, NNC = TN)
  are wrong for this simulation, where most NNC are genuine.

So PanIsoGuard is not a more accurate filter. SQANTI3's rules filter also says why it
discards an isoform (`*_filtering_reasons.txt`: the descriptor that failed). What PanIsoGuard
adds is evidence that filter does not read (the long-read alignments, the sample's genome,
other callers), a record for every isoform and junction including the ones it keeps, and
explicit abstention.

## Fixed along the way

The first run of this comparison exposed a projection bug: `PARTIAL` short-read support
(some but not all novel junctions corroborated) was promoted to `NOVEL_PROBABLE`.
Those were exactly PanIsoGuard's extra NNC false positives (oracle run: 6 vs 3 for the
short-read rule; 3 vs 3 after the fix). Since 0.0.4, `PARTIAL` to `UNCONFIRMED`.

## Caveats

- Single simulated truth set (chr22, PBSIM3 HiFi, FLAIR). Simulated short reads are clean
  (no intron retention, pre-mRNA, or mapping-hard regions), so `allSJ3` may look better
  here than on real data.
- The BAM axis thresholds are HiFi-calibrated and over-fire on ONT
  ([../gm12878_realdata](../gm12878_realdata)).

## Reproduce

```bash
bash ../sqanti_sim/run.sh     # the truth set (once)
bash run.sh                   # edit the paths at the top first; MC= and CHAIN= add IsoQuant to the count baselines
```

## Per-class precision (`per_class.py`)

The same truth, split by PanIsoGuard class (1 M simulated pairs; `pig` = no BAM, `pigbam` = with
the BAM). The classes are not an ordered scale: without the BAM, `UNRESOLVED` (mostly NIC
isoforms with no novel junction to test) is more precise than `UNCONFIRMED` at every depth
(250 k, 1 M, 4 M pairs and the oracle give the same order).

| class | genuine | false | precision (no BAM) | genuine | false | precision (BAM) |
|---|---:|---:|---:|---:|---:|---:|
| NOVEL_CONFIRMED | 357 | 3 | 0.992 | 351 | 3 | 0.992 |
| NOVEL_PROBABLE | 41 | 1 | 0.976 | 198 | 1 | 0.995 |
| KNOWN_FRAGMENT | n/a | n/a | n/a | 153 | 11 | 0.933 |
| UNRESOLVED | 160 | 7 | 0.958 | 0 | 6 | 0.000 |
| UNCONFIRMED | 169 | 61 | 0.735 | 23 | 36 | 0.390 |
| ARTIFACT | 3 | 68 | 0.042 | 5 | 83 | 0.057 |

With the BAM (0.0.7), the genuine `UNRESOLVED` NIC isoforms move to `NOVEL_PROBABLE` through
whole-chain read support (151 of the 198), the ISM isoforms carried end to end move to
`KNOWN_FRAGMENT` (153 of the 154 genuine and 11 of the 13 false: the simulator does not model
intra-priming, so the 3' sign cannot separate them here), and `UNRESOLVED` keeps the 6 false
ones. Without the BAM `UNCONFIRMED` is two populations: ISM pass-through (154 / 13, 0.92) and
novel isoforms with weak support (15 / 48, 0.24); with the BAM what remains is the weak-support
novel set (23 / 36, 0.39), below the 0.83 base rate. Numbers from the 0.1.1 binary (ruleset
`builtin-0.0.8`); `per_class_v011.txt` of `run.sh`'s work dir has every depth.
