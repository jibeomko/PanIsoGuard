# Whole-chain long-read support: what it lifts, for which caller

Since 0.0.7 `adjudicate --bam` counts, per isoform, the primary reads whose alignment carries
the **entire intron chain** (`evidence.chain_n_full_reads`). When short-read support is
`UNKNOWN` (no `SJ.tab`, or a NIC isoform whose junctions are all individually known, so
there is no novel junction to test) and no artifact sign is present, ≥ 3 such reads lift the
call from `UNRESOLVED` to `NOVEL_PROBABLE` (`chain_support`). The reads are the caller's own,
so the lift stops at `PROBABLE` and never overrides a short-read verdict.

This protocol asks what that rule does on a truth set, per caller: five callers on one chr22
SQANTI-SIM alignment (the [multicaller](../multicaller) set), each adjudicated with the axis
ON (default) and OFF (`min_full_reads` raised beyond any read count), with 1 M simulated
short-read pairs (`sr1M`) and with none (`nosr`). Positive = `NOVEL_CONFIRMED` /
`NOVEL_PROBABLE`; labels as in [sqanti3_filter_h2h](../sqanti3_filter_h2h).

## Result

| caller | short reads | n (genuine) | precision OFF → ON | recall OFF → ON | lifted to confident: genuine / false | selective risk @ coverage, OFF → ON |
|---|---|---:|---|---|---:|---|
| FLAIR | 1 M | 870 (730) | 0.990 → 0.993 | 0.545 → **0.745** | **146 / 0** | 0.024 @ 0.57 → 0.019 @ 0.73 |
| FLAIR | none | 870 (730) | — → 0.993 | 0.000 → **0.748** | 546 / 4 | — @ 0.00 → 0.007 @ 0.63 |
| IsoQuant | 1 M | 618 (607) | 0.988 → 0.991 | 0.526 → 0.708 | 111 / 0 | 0.018 @ 0.53 → 0.014 @ 0.71 |
| IsoQuant | none | 618 (607) | — → 0.991 | 0.000 → 0.708 | 430 / 4 | — @ 0.00 → 0.009 @ 0.70 |
| ESPRESSO | 1 M | 573 (555) | 0.985 → 0.987 | 0.697 → 0.935 | 132 / 1 | 0.025 @ 0.70 → 0.021 @ 0.93 |
| ESPRESSO | none | 573 (555) | — → 0.987 | 0.000 → 0.939 | 521 / 7 | — @ 0.00 → 0.013 @ 0.92 |
| Bambu | 1 M | 552 (409) | 0.725 → 0.748 | 0.748 → 0.951 | 83 / 15 | 0.281 @ 0.78 → 0.257 @ 0.96 |
| Bambu | none | 552 (409) | — → 0.764 | 0.000 → 0.951 | 389 / 120 | — @ 0.00 → 0.236 @ 0.92 |
| TALON | 1 M | 1,760 (697) | 0.700 → 0.709 | 0.485 → 0.697 | 148 / 54 | 0.299 @ 0.29 → 0.290 @ 0.40 |
| TALON | none | 1,760 (697) | — → 0.716 | 0.000 → 0.709 | 494 / 196 | — @ 0.00 → 0.284 @ 0.39 |

Coverage = fraction of non-FSM isoforms that received a decisive class (`NOVEL_CONFIRMED`,
`NOVEL_PROBABLE`, `ARTIFACT`, `REFERENCE_BIAS`); selective risk = wrong decisive calls among
them. Full numbers, including every OFF→ON transition by truth label, in
[../results/chain_support/metrics.json](../results/chain_support/metrics.json).

## Reading

- **Where the recall was.** With 1 M short-read pairs, FLAIR's recall goes from 0.545 to
  0.745 at precision 0.993: 146 genuine isoforms lifted, 0 false. They are the NIC isoforms
  whose junctions are all individually known, which the short-read axis cannot test (26 of
  them carry SQANTI3's RT-switch sign, which ruleset 0.0.5 no longer lets demote the lift;
  the 4 with an intra-priming sign stay `UNCONFIRMED`). The same for IsoQuant (111 / 0) and
  ESPRESSO (132 / 1).
- **Without short reads** PanIsoGuard confirms nothing (and, since ruleset 0.0.5, rejects
  nothing either: the BAM-only `ARTIFACT` calls it used to make were 36 % wrong). With chain
  support it decides 63–92 % of the isoforms at 0.7–1.3 % risk for FLAIR, IsoQuant and
  ESPRESSO.
- **Where it costs precision.** Bambu and TALON's false novel chains are carried by reads
  too (they pass alignment artifacts through), so the lift adds false positives there: 15
  and 54 with short reads, 120 and 196 without. Their precision is low to begin with
  (0.70–0.75) and barely moves, but a user of a permissive caller should know that read
  support cannot separate what the caller did not. Caller-level false chains (a caller
  stitching junctions no single read carries) are what the signal removes; alignment-level
  ones (every read carries the wrong chain) it cannot.
- **Not independent evidence.** The same BAM made the isoforms. That is why the rule stops
  at `NOVEL_PROBABLE`, why a short-read verdict (`UNSUPPORTED`, `PARTIAL`) is never
  overridden, and why a mapping artifact still preempts it. An independent library or a
  molecule-level hold-out would turn this into experimental corroboration; neither exists in
  the tracked data.

## A second truth set: SIRV-Set4

The same read-level signal on the Lexogen SIRV-Set4 spike-in reference (PBSIM3 reads, the
[multicaller](../multicaller) SIRV run; truth = the 13 SIRV transcripts hidden from the
reduced annotation; `score_sirv_chain.py`):

| caller | novel chains (genuine / false) | genuine with ≥ 3 full-chain reads | false with ≥ 3 |
|---|---:|---:|---:|
| FLAIR | 190 (12 / 178) | 8 / 12 | 13 / 178 |
| IsoQuant | 10 (6 / 4) | 5 / 6 | 0 / 4 |

The separation replicates (medians 11.5 vs 1.0 reads for FLAIR). The engine-level ON/OFF
comparison is **not evaluable** on SIRV: SQANTI3 puts 152 of FLAIR's 228 SIRV models, 7 of
the 12 genuine ones among them, into `antisense` / `intergenic` (the SIRV annotation's strand
against FLAIR's), which PanIsoGuard does not adjudicate, and the adjudicated NNC carry
mapping or non-canonical signs from the simulated reads. So on SIRV the rule lifts nothing
either way (recorded as `sirv_set4` in the envelope); the evidence it rests on is the
read-level table above.

## Reproduce

```bash
# SQANTI3 classification per caller (FLAIR's comes from ../sqanti_sim; the others as in ../multicaller):
#   sqanti3_qc.py --isoforms <caller.gtf> --refGTF chr22_modified.gtf --refFasta chr22.fa -o <caller> -d sq_<caller> --report skip
WORK=/path/to/sqanti_sim/work MC=/path/to/multicaller/data SJ=/path/to/h2h_sr/sr_1000000/SJ.out.tab \
  GENOME=/path/to/chr22.fa OUT=/tmp/chain_support bash run.sh
```

```bash
# SIRV-Set4 replication (needs pysam):
score_sirv_chain.py --work /path/to/pig_sirv/work --mc /path/to/pig_sirv_mc --out /tmp/sirv_chain \
  --pig ../../build/panisoguard --cfg ../../config/rules.default.toml --emit-metrics ../results/chain_support/metrics.json
```

`score_chain.py` imports the labelling and metrics of
[../sqanti3_filter_h2h/score_h2h.py](../sqanti3_filter_h2h/score_h2h.py). Public/simulated
data only. The axis itself: `src/evidence/bam_features.cpp` (`BamReader::chain_support`),
`src/core/rules.cpp` (the `UNKNOWN` branch), `config/rules.default.toml` (`[axis_chain]`).
