# Merge-tool head-to-head — is the consensus claim a matcher artifact?

The [multicaller](../multicaller) result rests on **how caller agreement is counted**.
`panisoguard combine` groups isoforms by **exact intron-chain fingerprint**. Two fair
objections a reviewer will raise:

1. **"`combine` just re-derives `gffcompare -i`."** Is it a novel algorithm or a
   re-implementation of the established N-way comparison?
2. **"Exact matching splits true cross-caller agreements that a fuzzy merge (TAMA)
   would keep — so the ≥k precision is an artifact of your matcher, not the data."**

This protocol settles both by putting every matcher on **one footing**: each produces
caller-support groups over the **same 5 chr22 caller GTFs** (FLAIR + IsoQuant + Bambu +
ESPRESSO + TALON), and we score each with the **identical** "≥ k callers" PR curve vs
SQANTI-SIM truth. Matchers compared: `combine` (exact), **gffcompare -i** (the incumbent),
**TAMA merge** (fuzzy, junction wobble 10 bp), and a controlled **exact→fuzzy wobble
sweep** (0–20 bp).

## Result 1 — `combine` ≡ `gffcompare -i` (exactly)

| matcher | novel chains | ≥2 P/R | ≥3 P/R | n_callers dist (1·2·3·4·5) |
|---------|-------------:|:------:|:------:|:--------------------------:|
| **PanIsoGuard `combine` (exact)** | 1982 | 0.852 / 0.971 | **0.980 / 0.910** | 1105·163·254·128·332 |
| **gffcompare -i** | 1982 | 0.852 / 0.971 | **0.980 / 0.910** | 1105·163·254·128·332 |
| wobble 0 bp (control) | 1982 | 0.852 / 0.971 | 0.980 / 0.910 | 1105·163·254·128·332 |

`combine` reproduces gffcompare's grouping **byte-for-byte** (same novel count, same
n_callers distribution, same PR curve). So `combine` is **not a novel merge algorithm** —
it is a clean, **htslib-native re-implementation** of the established intron-chain N-way
comparison that feeds the adjudicator directly (no Perl, no GTF round-trip, no extra
dependency). We state this honestly rather than imply novelty.

## Result 2 — the consensus precision is matcher-robust (not an exact-matching artifact)

| matcher | ≥2 precision | ≥3 precision |
|---------|:-----------:|:-----------:|
| combine / gffcompare (exact) | 0.852 | **0.980** |
| TAMA (fuzzy, 10 bp) | 0.933 | **0.976** |
| wobble 5 / 10 / 20 bp | 0.839 / 0.837 / 0.836 | **0.975 / 0.975 / 0.975** |

≥3-caller novel precision is **0.975–0.980 across every matcher** — exact, the incumbent,
the fuzzy standard, and the whole wobble range. The headline ("requiring more callers
raises precision") is a property of the **data** (the callers genuinely disagree on
artifacts), not of PanIsoGuard's exact matcher.

## Result 3 — exact matching does not lose agreements at the boundary that matters

The decisive metric: how many **genuine** novels does exact matching call *single-caller*
that a fuzzy matcher upgrades to *≥ 2 callers*?

| fuzzy matcher | genuine novels crossing single → ≥2 | genuine upgraded (any amount) |
|---------------|:----:|:----:|
| gffcompare | **0** | 0 |
| TAMA (10 bp) | **0** | 55 |
| wobble 5 / 10 / 20 bp | **0** | 38 / 52 / 67 |

**Zero** at every wobble level and for TAMA. The chains a fuzzy merge upgrades (55 for
TAMA) all move *within* the multi-caller pile (e.g. 3→4 callers); **none** crosses the
single↔multi boundary the consensus gate depends on. The single-caller artifacts are
genuinely single-caller — no other caller reports anything within 20 bp — so fuzzy
merging does not rescue them and does not inflate the ≥k precision. (TAMA's more
aggressive merge collapses the call set to 974 novel chains and trades some recall at ≥3,
0.831 vs 0.910 — a different operating point, not a more accurate one.)

## Honest takeaway

- Multi-caller consensus is **established field practice** (LRGASP; *Nat Commun* 2024),
  not a PanIsoGuard invention. `combine` **operationalizes it correctly** — provably
  equal to gffcompare, robust to the exact-vs-fuzzy choice — and feeds it into the
  adjudicator as one evidence axis.
- The consensus axis is therefore a **sound, validated** component, but **not the
  differentiator**. PanIsoGuard's genuine moat is the **reference-bias rescue + circularity
  firewall** (no surveyed merge/QC tool does this) and the unified auditable verdict. See
  [../../docs/relationship_to_merge_tools.md](../../docs/relationship_to_merge_tools.md).

## Reproduce

```bash
conda create -n mergecmp -c conda-forge -c bioconda gffcompare ucsc-gtftogenepred ucsc-genepredtobed
gffcompare -i gtflist.txt -o gffcmp                          # incumbent N-way
# TAMA (github.com/GenomeRIK/tama, python2): GTF->BED12 (gene_id;transcript_id in col4)
python2 tama_merge.py -f filelist.txt -p tama -a 50 -m 10 -z 50 -d merge_dup
python3 score_merge.py --callers callers.tsv --truth truth.truth.tsv --ref chr22_modified.gtf \
  --paniso matrix5.tsv --gffcompare gffcmp.tracking --tama tama_trans_report.txt \
  --wobble 0,5,10,20 --emit-metrics ../results/merge_comparison/metrics.json
```

All inputs are the **same public/simulated** SQANTI-SIM chr22 5-caller GTFs as the
[multicaller](../multicaller) protocol; no private data. Result envelope:
[../results/merge_comparison/metrics.json](../results/merge_comparison/metrics.json).
