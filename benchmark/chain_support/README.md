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

## Decoy chains: how often does the lift fire on a chain the annotation does not know?

`panisoguard decoys` recombines each gene's annotated transcripts: the 5′ introns of one joined
to the 3′ introns of another at an intron they share, keeping only chains that are neither
annotated nor a contiguous piece of an annotated chain (at most 10 per gene, seed 1). Every
junction is known, so these are exactly the NIC chains that only the chain axis can lift, and
their whole-chain reads are counted with the same code the axis uses. `score_decoys.py`
summarizes the runs and, on SQANTI-SIM, splits the decoys by the truth.

| dataset (annotation) | decoys | ≥ 3 reads | both parents ≥ 3 → ≥ 3 reads | annotated parents ≥ 3 |
|---|---:|---:|---:|---:|
| SQANTI-SIM chr22, reduced annotation (the truth set) | 3,045 | 35 (1.1 %) | 3 / 29 | 272 / 3,397 (8 %) |
| same, every decoy (no cap) | 40,471 | 89 (0.2 %) | 6 / 214 | 581 / 6,278 (9 %) |
| SQANTI-SIM chr22, full annotation | 3,198 | **0** | 0 / 45 | 479 / 3,572 (13 %) |
| LRGASP WTC11 PacBio cDNA (GENCODE v49) | 139,003 | 3,415 (2.5 %) | 1,872 / 4,323 (**43 %**) | 21,638 / 153,588 (14 %) |
| LRGASP WTC11 ONT cDNA | 139,003 | 848 (0.6 %) | 411 / 1,259 (33 %) | 8,205 (5 %) |
| HG03516 Kinnex Iso-Seq | 139,003 | 1,996 (1.4 %) | 976 / 2,106 (46 %) | 12,825 (8 %) |
| HG02717 Kinnex Iso-Seq | 139,003 | 4,043 (2.9 %) | 2,308 / 5,516 (42 %) | 22,729 (15 %) |
| GM12878 ONT direct RNA | 139,003 | 7 (0.01 %) | 0 / 5 | 884 (0.6 %) |

On the reduced SQANTI-SIM annotation the truth tells the two kinds of decoy apart: 31 of the 35
lifted decoys are deleted genuine transcripts (31 / 31 of the genuine decoys are lifted), and
the wrong ones reach 3 reads 4 / 3,014 times (0.13 %; 19 / 40,401 = 0.05 % without the cap;
0 / 26 and 2 / 210 among decoys whose parents are both expressed). With the full annotation,
where nothing expressed is unannotated, 0 / 3,198.

- **In simulation the lift almost never fires on a wrong combination of known junctions**
  (0–0.1 %) and always on the right ones. The decoy check recovers the truth set's own result
  without the truth.
- **On real data a third to a half of the recombinants of two expressed isoforms are carried by
  ≥ 3 reads.** They do not look like RT or PCR chimeras: the lifted decoys carry about as many
  reads as their parents (median 7 against the weaker parent's 5; read-count ratio to that
  parent median 0.92–1.00, under 0.05 in fewer than 1 %), whereas a chimera is a small
  fraction of the parent molecules. They look like combinatorial isoforms GENCODE does not
  list, in weakly expressed genes where 3 reads is most of the evidence there is. On real data
  the decoy rate is therefore an upper bound that mostly measures annotation incompleteness;
  the simulation, which has no unannotated isoforms, gives the lower bound.
- **They reproduce in an independent library like the annotated transcripts do.** The same
  annotation and seed give the same 139,003 decoys on every dataset, so a decoy lifted in one
  library can be looked up in another (`score_decoys.py --cross`). Of the decoys lifted in
  LRGASP PacBio, 18 % also have ≥ 3 reads in the ONT cDNA library of the same cells, against
  30 % of the annotated parents lifted in PacBio (the ONT library is shallower); from ONT to
  PacBio 73 % against 79 %; between the two Kinnex individuals 86 % / 92 % and 42 % / 52 %.
  Restricted to decoys whose two parents are expressed in both libraries, 75–93 % of the
  lifted ones are lifted in the other library too. A chimera made during one library's
  reverse transcription would not recur in another platform or another person at the rate
  of the real transcripts.
- **Scope.** A decoy is a wrong *combination* of right junctions. Chains with a wrong junction
  (the shifted-coordinate, alignment-level false chains that cost Bambu and TALON precision
  above) are outside its reach by construction.
- **What it says about the rule.** ≥ 3 whole-chain reads is rarely an accident of combination
  (and a 4-bp shifted junction is a different problem), but it cannot tell a biologically
  distinct isoform from a rare recombinant of two expressed ones: that is why the lift stops at
  `NOVEL_PROBABLE`. The conditional column is a per-dataset knob: run `decoys` on your own
  BAM to see how combinatorial your expressed genes are at your depth before trusting the
  lifted NIC calls as a set. 1–12 min genome-wide, 0.4 GB.

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

```bash
# decoy chains (any dataset; the truth split needs the SQANTI-SIM truth table):
panisoguard decoys --ref-gtf /path/to/chr22_modified.gtf --bam /path/to/aln.bam --out /tmp/decoys/sqsim
score_decoys.py --run sqsim=/tmp/decoys/sqsim --truth sqsim=/path/to/truth.truth.tsv \
  --run pb=/tmp/decoys/lrgasp_pb --run ont=/tmp/decoys/lrgasp_ont --cross pb,ont \
  --emit-metrics ../results/chain_support/metrics.json
```

`score_chain.py` imports the labelling and metrics of
[../sqanti3_filter_h2h/score_h2h.py](../sqanti3_filter_h2h/score_h2h.py). Public/simulated
data only. The axis itself: `src/evidence/bam_features.cpp` (`BamReader::chain_support`),
`src/core/rules.cpp` (the `UNKNOWN` branch), `config/rules.default.toml` (`[axis_chain]`).
