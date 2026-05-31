# End-to-end read-simulation E1

The realistic, messy complement to [`controlled_truth`](../controlled_truth) and
[`synthetic_axes`](../synthetic_axes): instead of hand-built inputs, it runs the full
upstream pipeline and checks that PanIsoGuard behaves correctly on **real FLAIR +
SQANTI3 output** derived from simulated reads.

```
reference (1 chromosome) ─ hide ~20% transcripts (= genuine-novel truth)
   │
   ├─ PBSIM3  ─ simulate long reads from all chosen transcripts
   ├─ minimap2 ─ splice-aware alignment
   ├─ FLAIR    ─ collapse reads into isoforms
   ├─ SQANTI3  ─ structural QC vs the REDUCED annotation
   └─ PanIsoGuard adjudicate
            │
            └─ score: an isoform's intron chain matching a HIDDEN transcript = genuine;
               a novel chain matching nothing real = false (FLAIR/alignment artifact)
```

## Run

```bash
# edit the tool/env paths at the top of run_e2e.sh, then:
bash run_e2e.sh e2e_work
```

Tools (each typically its own conda env): `pbsim3`, `minimap2`, `samtools`,
`FLAIR 3.x`, `SQANTI3 5.x`, and the built `panisoguard`.

> **FLAIR caveat:** FLAIR 3.0.0's `bam2Bed12` entry point is broken on bioconda
> (`ModuleNotFoundError: No module named 'flair.bam2Bed12'`); `bam2bed12.py` here is a
> drop-in BAM→BED12 converter that bypasses it and feeds `flair collapse`.
>
> **Truth caveat:** rigorous read-level truth bookkeeping is SQANTI-SIM's job; this
> harness labels truth by intron-chain match to the hidden transcripts instead, which
> is sufficient to score genuine-vs-artifact novelty.

## Result (GENCODE v49 chr22; 150 transcripts, 30 hidden; 6000 PBSIM3 reads)

FLAIR produced 397 isoforms (5019/6000 reads aligned). On the SQANTI-novel
(NIC/NNC) set, scored by intron-chain truth:

```
            pred_genuine  pred_false  abstain
genuine            4           0         16
false              1         241          7

specificity(false) = 0.996   genuine-recall(decisive) = 1.000   precision(genuine) = 0.800
```

- **PanIsoGuard rejected 241/242 (99.6%) of FLAIR's false-novel artifacts** — on real,
  messy caller output it removes essentially all spurious novelty.
- **No decisive genuine novel was wrongly rejected** (recall 1.0; 0 false negatives).
- **16/20 genuine novels were honestly abstained** (`AMBIGUOUS`): they are NIC-like,
  with every junction already in the reduced annotation, so there is no novel junction
  to corroborate — held, not guessed (the documented behaviour, also seen in
  `controlled_truth`).
- One false positive (a single artifact called genuine).

This validates PanIsoGuard end-to-end on real FLAIR/SQANTI3 output; the clean
per-axis correctness is established in the sibling controlled benchmarks.
