# Multi-caller consensus — caller-agnostic integration (5 callers)

The core value of PanIsoGuard is that it sits **above** any single isoform caller.
This protocol makes that concrete: it runs **five** long-read callers on the **same**
alignment, integrates them by intron-chain fingerprint with `panisoguard combine`, and
shows that **multi-caller agreement is a strong, caller-agnostic confidence signal** —
one a single-caller + SQANTI3 pipeline cannot produce.

## Why this matters

A novel isoform call from one caller carries **that caller's** false-discovery rate.
The five callers disagree enormously, each with a different precision/recall profile
(SQANTI-SIM GENCODE v49 chr22, one shared minimap2 alignment, scored vs truth):

| caller | novel calls | genuine | false | precision | recall |
|--------|------------:|--------:|------:|----------:|-------:|
| FLAIR (collapse)        |  862 | 723 |  139 | 0.839 | 0.940 |
| IsoQuant (3.x models)   |  618 | 607 |   11 | 0.982 | 0.789 |
| Bambu (3.x, NDR 0.5)    |  552 | 409 |  143 | 0.741 | 0.532 |
| ESPRESSO (1.4)          |  573 | 555 |   18 | 0.969 | 0.722 |
| TALON (6.0, filtered)   | 1760 | 697 | 1063 | 0.396 | 0.906 |

Precision spans **0.40 (TALON) to 0.98 (IsoQuant)**; recall spans 0.53 to 0.94. There is
no single "right" caller, and you cannot tell *in advance* which of a caller's novel
calls are the artifacts.

## The consensus result

`combine` integrates the five caller sets into **1,982 unique novel multi-exon intron
chains**. Stratifying by **how many callers** recovered each one — vs SQANTI-SIM ground
truth — precision climbs steeply and monotonically with agreement:

| # callers supporting a novel chain | genuine | false | precision |
|-----------------------------------:|--------:|------:|----------:|
| 1 (single-caller)                  |  13 | 1092 | **0.012** |
| 2                                  |  47 |  116 | 0.288 |
| 3                                  | 245 |    9 | 0.965 |
| 4                                  | 126 |    2 | 0.984 |
| 5 (unanimous)                      | 329 |    3 | **0.991** |

A novel chain seen by **only one** of five callers is **almost always an artifact**
(precision 0.012 — 1,092 false vs 13 genuine). Agreement is the signal.

## Precision–recall curve over the consensus threshold

Requiring "≥ k callers" is a tunable gate. Sweeping k traces a clean PR curve:

| threshold | genuine | false | precision | recall | F1 |
|-----------|--------:|------:|----------:|-------:|---:|
| ≥ 1 caller  | 760 | 1222 | 0.383 | 0.988 | 0.553 |
| ≥ 2 callers | 747 |  130 | 0.852 | 0.971 | 0.908 |
| **≥ 3 callers** | 700 | 14 | **0.980** | **0.910** | **0.944** ← max F1 |
| ≥ 4 callers | 455 |    5 | 0.989 | 0.592 | 0.740 |
| ≥ 5 callers | 329 |    3 | 0.991 | 0.428 | 0.598 |

Key finding: **the F1-optimal threshold scales with the number of callers.** With five
callers the best operating point is **≥ 3** (precision 0.98 at recall 0.91) — more
callers means you should *require more agreement*. `consensus_min_callers` is a
configurable rule threshold (TOML); this curve is exactly the calibration a user needs to
pick it. (With three callers, ≥ 2 was already near-optimal.)

## The engine reproduces it (`adjudicate --caller-support`)

The wired **consensus axis** consumes the `combine` matrix. Adjudicating FLAIR's
isoforms **long-read-only** (no short-read `SJ.tab`), with vs without the 5-caller matrix
(default `consensus_min_callers = 2`):

| setting | confident-novel recall | precision |
|---------|-----------------------:|----------:|
| FLAIR alone (no cross-caller info) | 0.000 (all novel held `AMBIGUOUS`) | — |
| **FLAIR + consensus** | 0.685 | **0.975** |

Without orthogonal short-read evidence a single caller's novels are all `AMBIGUOUS` (the
engine refuses to guess). The consensus axis promotes the **593** novels that ≥ 2 callers
agree on (513 → `MEDIUM_CONF_NOVEL`, 80 → `LOW_CONF_PARTIAL`); that promoted set is
**97.5 % precise** against truth.

## Integration device (caller-agnostic merge)

The five callers use completely unrelated native ID schemes, yet `combine` collapses
identical splice chains onto one stable `PIG.NNNNNN` id with every native ID preserved —
agreement is measured on the **splice chain**, not on names:

```
PIG.000233  callers=bambu,espresso,flair,isoquant,talon  novelty=novel
  native_ids: bambu=BambuTx1 ; espresso=ESPRESSO:chr22:2:0 ;
              flair=ENST00000779064.1_..read_1 ; isoquant=transcript2.chr22.nnic ;
              talon=TALONT000010833
```

(Complements the synthetic cross-caller unit test in `tests/unit/consensus_test.cpp`.)

> **`combine` is not a novel merge — and the consensus claim is not a matcher artifact.**
> Multi-caller consensus is established practice (gffcompare, TAMA, LRGASP). The
> [merge_comparison](../merge_comparison) protocol shows `combine` reproduces
> `gffcompare -i` **exactly**, and the ≥3 precision above is **matcher-robust** (0.975–0.980
> across exact / gffcompare / TAMA / a 0–20 bp wobble sweep, with **0** genuine novels
> crossing the single↔multi boundary under any fuzzy matcher). So this protocol measures
> a real property of the *data*; PanIsoGuard's differentiator is the reference-bias rescue,
> not the merge (see [relationship_to_merge_tools.md](../../docs/relationship_to_merge_tools.md)).

## Pipeline

```
shared minimap2 alignment (PBSIM3 reads, chr22)
  ├─ FLAIR collapse             ─┐
  ├─ IsoQuant transcript_models ─┤
  ├─ Bambu (novel BambuTx)      ─┼ panisoguard combine --gtf <caller:gtf>...
  ├─ ESPRESSO (novel_isoform)   ─┤   --ref-gtf chr22_modified.gtf  →  caller-support matrix
  └─ TALON (filtered whitelist) ─┘        │
        adjudicate --caller-support matrix  →  consensus-promoted verdicts
        score_multicaller.py (vs SQANTI-SIM truth)  →  results/multicaller/metrics.json
```

See [run.sh](run.sh) for the exact commands and [score_multicaller.py](score_multicaller.py)
for the scorer (which also emits the `pr_curve`). Result envelope:
[../results/multicaller/metrics.json](../results/multicaller/metrics.json).

## Reproduce

All inputs are public/simulated (SQANTI-SIM GENCODE v49 chr22; no private data). The five
callers run on the shared alignment from the [sqanti_sim](../sqanti_sim) protocol; needs
the `isoquant`, `bambu` (bioconductor-bambu + r-biocmanager), `espresso`, and `talon`
conda envs. Notes baked into [run.sh](run.sh): TALON needs MD-tagged SAM (`samtools calmd`)
and the standard `talon_filter_transcripts` whitelist; Bambu uses a fixed NDR to skip its
online recommendation step; TALON 6.0 on py3.7 needs a `pyranges`/`importlib.metadata`
shim.

## Generalization to a second reference (SIRV-Set4)

To check the consensus signal is **not specific to GENCODE chr22**, the same analysis was
run on the **Lexogen SIRV-Set4** spike-in reference — a denser, multi-contig transcriptome
(7 loci, ~70 dense overlapping isoforms) with its own hidden-chain truth (FLAIR + IsoQuant
+ Bambu + ESPRESSO; [score_sirv.py](score_sirv.py)). The same pattern holds:

| supporting callers | precision (vs SIRV hidden truth) |
|--------------------|:---:|
| 1 (single-caller)  | **0.037** |
| ≥ 2                | **1.000** |

Single-caller novels are again overwhelmingly artifacts (precision 0.04, mirroring the 0.012
on chr22), and ≥ 2-caller agreement is perfectly precise. So "agreement = confidence" is a
property of multi-caller data, not of one reference or simulator. Recorded in
[../results/sirv_multicaller](../results/sirv_multicaller).

## Scope / what remains

Validated for **precision stratification** on one controlled truth set (chr22, five
callers). The consensus axis is deliberately a *weaker, methodological* signal than
orthogonal short-read corroboration (it promotes to `MEDIUM`, never `HIGH`); when a
short-read `SJ.tab` is supplied it takes precedence. A whole-genome real-data extension
(no truth → orthogonal canonical-motif validation) is the companion
[wholegenome_multicaller](../wholegenome_multicaller) protocol.
