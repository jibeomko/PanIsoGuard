# Whole-genome multi-caller consensus on real data (GM12878)

Extends the chr22 [multicaller](../multicaller) protocol from a controlled simulation to
**real, whole-genome, fully public** data: ENCODE GM12878 ONT direct-RNA, three callers
(IsoQuant, Bambu, ESPRESSO) run on one shared genome alignment and integrated by
`panisoguard combine`. **This is an honest, small-N result — read it with its caveats.**

## No truth on real data → orthogonal validation

Real data has no SQANTI-SIM ground truth, so the chr22 precision/recall scoring does not
apply. Instead we use a caller-**independent** quality signal: **splice-motif
canonicality**. Genuine spliceosomal introns are overwhelmingly canonical (GT-AG, or the
minor GC-AG / AT-AC); alignment/caller artifacts are far more often non-canonical. If
multi-caller agreement is a real quality signal, the fraction of novel chains whose
*every* junction is canonical should rise with the number of supporting callers.

## What we found

`combine` integrated the three callers into **285 unique novel multi-exon chains**
(whole genome, GENCODE v49). Two things stand out — one strong, one modest:

**1. Disagreement holds at genome scale (strong, mirrors chr22).**

| supporting callers | novel chains |
|--------------------|-------------:|
| 1 (single-caller)  | 254 |
| 2                  | 30 |
| 3 (unanimous)      | 1 |

Of 285 novel chains, **254 (89 %) are single-caller** and only **31 are recovered by
≥ 2** callers. Even on real data the callers agree on a small minority of novel calls —
the same motivation for integration as the simulation. (Note also the *absolute* count is
small: GM12878 is a well-characterized cell line, so most transcripts are already known
and genuine novel discovery is modest — IsoQuant 113, Bambu 195, ESPRESSO 9 novel.)

**2. Quality direction is correct, but the headroom is small (ceiling effect).**

| novel-call set | all-canonical fraction |
|----------------|-----------------------:|
| single-caller (exactly 1) | 248/254 = **0.976** |
| consensus (≥ 2 callers)   | 31/31 = **1.000** |

Consensus reaches 100 % canonical; the 6 non-canonical chains are all single-caller
(exactly the kind of artifact consensus excludes). But single-caller chains are *already*
97.6 % canonical, so the absolute lift is only **+0.024**.

## Honest reading

The strong artifact-filtering on chr22 (single-caller precision **0.012**) was driven
largely by the **liberal** end of the caller spectrum — unfiltered TALON's 1,063 false
novels. Here all three callers run at **recommended settings** and are internally
stringent, so they emit few non-canonical artifacts and there is little for consensus to
filter — a ceiling effect, on a small novel set.

So on **real data with production callers**, the practical value of multi-caller consensus
is **less "remove artifacts" and more "identify the reproducible high-confidence core"**:
the 31 cross-caller-agreed novels are the ones you would report with highest confidence,
while the 254 single-caller novels — though mostly canonical — are unreproducible across
methods and warrant the `AMBIGUOUS`/lower-confidence treatment the engine gives them. The
disagreement itself (89 % single-caller) is the headline finding that generalizes from
simulation to real, whole-genome data. (Consistent with the project's honest-results ethos,
cf. [pangenome_public](../results/pangenome_public).)

## Decision impact — what you would *report* changes 9× with the policy

The point above, stated as the decision a user actually faces. On this **real** GM12878
data, the number of novel isoforms you would put in a paper depends entirely on how you
integrate the three callers:

| integration policy | novel isoforms reported |
|---------------------|------------------------:|
| Bambu alone | 195 |
| IsoQuant alone | 113 |
| ESPRESSO alone | 9 |
| **naive union** (trust *any* caller) | **285** |
| **PanIsoGuard consensus** (≥ 2 callers) | **31** |

The reported count swings **31 → 285 (9.2×)** purely on policy. The gap is **caller-private
novels** — calls made by exactly one of the three callers and corroborated by **none** of
the others: **254 of the 285 union chains (89 %)** are private (Bambu 164, IsoQuant 84,
ESPRESSO 6). Cross-caller agreement is sparse and asymmetric — the 31-chain consensus core
is co-supported mostly by Bambu+IsoQuant (29), with Bambu+ESPRESSO (3) and ESPRESSO+IsoQuant
(1) almost never agreeing, and a single unanimous chain.

This is the **conclusions-change case**: a naive multi-caller merge would report 285 novel
isoforms, 89 % of which no second method reproduces; PanIsoGuard reports the **31-chain
reproducible core** (also the only 100 %-canonical subset) and makes the policy **explicit
and auditable** via `n_callers` in the verdict, instead of burying a 9× discretion behind
an undocumented union. The choice is no longer hidden — it is a logged rule threshold
(`consensus_min_callers`).

## Reproduce

All inputs are fully public (ENCODE GM12878 ONT dRNA `ENCSR368UNC`; GRCh38/GENCODE v49);
no private data. See [run.sh](run.sh). The BAM is first restricted to the main chromosomes
(chr1-22/X/Y/M) so ESPRESSO does not choke on alt/random contigs; callers run on that
shared alignment; `combine` integrates them; [score_wholegenome.py](score_wholegenome.py)
reads donor/acceptor dinucleotides from the genome (strand-aware) and stratifies the
all-canonical fraction by caller agreement. Result envelope:
[../results/wholegenome_multicaller/metrics.json](../results/wholegenome_multicaller/metrics.json).

## Scope / what remains

One real public sample (GM12878); K562 lacks a retained alignment in this workspace, and
re-alignment + a 5-caller whole-genome sweep (incl. TALON) is future work. The orthogonal
canonical-motif signal has a low ceiling on stringent-caller data; a higher-dynamic-range
orthogonal signal (e.g. CAGE/polyA proximity, short-read junction support) would sharpen
the real-data quality claim.
