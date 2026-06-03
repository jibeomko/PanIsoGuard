# Multi-caller consensus — caller-agnostic integration (FLAIR + IsoQuant + Bambu)

The core value of PanIsoGuard is that it sits **above** any single isoform caller.
This protocol makes that concrete: it runs **three** long-read callers on the **same**
alignment, integrates them by intron-chain fingerprint with `panisoguard combine`, and
shows that **multi-caller agreement is a strong, caller-agnostic confidence signal** —
one a single-caller + SQANTI3 pipeline cannot produce.

## Why this matters

A novel isoform call from one caller carries **that caller's** false-discovery rate.
The three callers disagree substantially, and each emits a different mix of genuine
discoveries and artifacts (SQANTI-SIM GENCODE v49 chr22, one shared minimap2 alignment):

| caller | novel calls | genuine | false | precision | recall |
|--------|------------:|--------:|------:|----------:|-------:|
| FLAIR (collapse)      | 862 | 723 | 139 | 0.839 | 0.940 |
| IsoQuant (3.x models) | 618 | 607 |  11 | 0.982 | 0.789 |
| Bambu (3.x, NDR 0.5)  | 552 | 409 | 143 | 0.741 | 0.532 |

There is no single "right" caller: FLAIR has the best recall, IsoQuant the best
precision, Bambu trails both — and you cannot tell *in advance* which of a caller's
novel calls are the artifacts.

## The consensus result

`combine` integrates the three caller sets into **1,023 unique novel multi-exon intron
chains**. Stratifying those chains by **how many callers** recovered each one — against
SQANTI-SIM ground truth — precision rises monotonically and steeply with agreement:

| # callers supporting a novel chain | genuine | false | precision |
|-----------------------------------:|--------:|------:|----------:|
| 1 (single-caller)                  | 141 | 263 | **0.349** |
| 2                                  | 220 |   9 | **0.961** |
| 3 (unanimous)                      | 386 |   4 | **0.990** |

A novel chain seen by **only one** caller is *more likely an artifact than real*
(precision 0.35). Requiring **≥ 2** callers is the decisive gate:

| novel-call set | genuine | false | precision | recall |
|----------------|--------:|------:|----------:|-------:|
| union (≥ 1 caller — "trust any caller") | 747 | 276 | 0.730 | 0.971 |
| **consensus (≥ 2 callers)**             | 606 |  13 | **0.979** | 0.788 |

The consensus gate raises novel-call precision **+0.25 (0.730 → 0.979)** for a recall
cost of ~0.18 — exactly the trade a confidence layer should expose.

## The engine reproduces it (`adjudicate --caller-support`)

The wired **consensus axis** consumes the `combine` matrix. Adjudicating FLAIR's
isoforms **long-read-only** (no short-read `SJ.tab`), with vs without the matrix:

| setting | confident-novel recall | precision |
|---------|-----------------------:|----------:|
| FLAIR alone (no cross-caller info) | 0.000 (all novel held `AMBIGUOUS`) | — |
| **FLAIR + consensus** | 0.560 | **0.978** |

Without orthogonal short-read evidence a single caller's novels are all `AMBIGUOUS`
(the engine refuses to guess). The consensus axis promotes the **483** novels that ≥ 2
callers agree on (418 → `MEDIUM_CONF_NOVEL`, 65 → `LOW_CONF_PARTIAL`), and that promoted
set is **97.8 % precise** against truth.

## Integration device (caller-agnostic merge)

The three callers use completely unrelated native ID schemes, yet `combine` collapses
identical splice chains onto one stable `PIG.NNNNNN` id with every native ID preserved:

```
pig_id      n_callers  callers                novelty  native_ids
PIG.000033  3          bambu,flair,isoquant   novel    bambu=BambuTx1; flair=ENST00000779064.1_..read_1; isoquant=transcript2.chr22.nnic
PIG.000038  3          bambu,flair,isoquant   novel    bambu=BambuTx3; flair=ENST00000754829.1_..read_2; isoquant=transcript18.chr22.nic
```

This is what makes the consensus count meaningful: agreement is measured on the **splice
chain**, not on names. (Complements the synthetic cross-caller unit test in
`tests/unit/consensus_test.cpp`.)

## Pipeline

```
shared minimap2 alignment (PBSIM3 reads, chr22)
  ├─ FLAIR collapse             ─┐
  ├─ IsoQuant transcript_models ─┤  panisoguard combine --gtf flair:.. --gtf isoquant:.. --gtf bambu:..
  └─ Bambu (novel BambuTx)      ─┘     --ref-gtf chr22_modified.gtf  →  caller-support matrix
                                            │
        adjudicate --caller-support matrix  →  consensus-promoted verdicts
        score_multicaller.py (vs SQANTI-SIM truth)  →  results/multicaller/metrics.json
```

See [run.sh](run.sh) for the exact commands and [score_multicaller.py](score_multicaller.py)
for the scorer. Result envelope: [../results/multicaller/metrics.json](../results/multicaller/metrics.json).

## Reproduce

All inputs are public/simulated (SQANTI-SIM GENCODE v49 chr22; no private data). Fill in
the placeholder paths at the top of [run.sh](run.sh) and run it; it needs the `isoquant`
and `bambu` (bioconductor-bambu + r-biocmanager) conda envs and the shared alignment from
the [sqanti_sim](../sqanti_sim) protocol.

## Scope / what remains

Validated for **specificity stratification** on one controlled truth set (chr22, three
callers). The consensus axis is deliberately a *weaker, methodological* signal than
orthogonal short-read corroboration (it promotes to `MEDIUM`, never `HIGH`); when a
short-read `SJ.tab` is supplied it takes precedence. Extending to additional callers
(ESPRESSO, TALON) and to whole-genome real data is future work.
