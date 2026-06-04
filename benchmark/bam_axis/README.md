# BAM mapping axis — wiring the indel-near / soft-clip read signals

The BAM read-level axis (`adjudicate --bam`) inspects the spanning reads at every novel
junction and measures four fractions: low-MAPQ, supplementary/secondary, **terminal
soft-clip**, and **indel-adjacent-to-the-junction**. Before v0.0.4 only the first two fed
the verdict; the soft-clip and indel-near fractions were computed and emitted to
`attribution.jsonl` for inspection but **never acted on**. This protocol asks whether they
*should* be — and wires the answer in.

## The mechanism (and why it is the "pangenome insertion" case)

A genuine spliceosomal intron aligns as a clean `N` (ref-skip) in the read CIGAR. A common
**false** novel junction is instead an **alignment-ambiguous indel** — a small insertion or
deletion that the aligner could equally render as a short intron, so the caller reports a
"novel junction" that is really an indel. That is exactly what `bam_frac_indel_near`
measures: the fraction of exact-spanning reads that carry an indel within the junction
window. (A sample **insertion** relative to the reference does not create a reference-
coordinate intron the way a deletion does — it surfaces as an indel/soft-clip near the
junction. So the BAM indel-near signal is how insertion-driven apparent-novelty gets
caught, complementing the deletion-driven reference-bias rescue of the variant/pangenome
axis.)

## Does it discriminate? (chr22 SQANTI-SIM truth)

Adjudicating the chr22 SQANTI-SIM FLAIR output **with `--bam`** and stratifying
`bam_frac_indel_near` by truth (novel junctions, reads actually spanning):

| novel-junction set | n | mean indel_near | max / p90 |
|--------------------|--:|----------------:|----------:|
| **genuine** | 430 | 0.032 | max **0.190** |
| **false**   | 124 | **0.350** | p90 **1.000** |

Genuine junctions top out at **0.190**; false ones pile up at **1.000** (≈half the false
set has an indel next to the junction on *every* spanning read). An 11× mean separation on
a signal the engine was throwing away.

**Soft-clip** is ≈0 on this clean PBSIM3 HiFi data (genuine max 0.042, false max 0.020) —
it does not discriminate here. It is wired as the same mapping mechanism (a heavily
soft-clipped junction = reads that could not align through) but its yield is **data-
dependent** and only expected to matter on noisier (e.g. ONT) reads; no yield is claimed
on this set.

## The operating point (zero genuine cost)

Because genuine `indel_near` never exceeds 0.190, **any** gate in [0.20, 0.80] flags zero
genuine junctions. At the shipped default **0.5** (matching the low-MAPQ / supplementary
gates):

| gate (`max_indel_near_frac`) | false caught | genuine flagged | precision |
|:---:|:---:|:---:|:---:|
| 0.5 | **41 / 124** | **0** | **1.000** |

## Impact — wired ON vs OFF (same binary, same inputs)

`adjudicate --bam` with the default config (ON) vs a config that raises the indel/softclip
gates above 1.0 (OFF, never fires):

| | OFF (low-MAPQ/supp only) | ON (+ indel-near) |
|---|:---:|:---:|
| false-novel **specificity** | 0.946 | **0.966** (+0.020) |
| genuine **sensitivity** | 0.589 | **0.589** (unchanged) |
| genuine isoforms lost from confident-novel | — | **0** |

False-novel transitions OFF→ON: **24** `LOW_CONF_PARTIAL → ARTIFACT` and **3**
`MEDIUM_CONF_NOVEL → LOW_CONF_PARTIAL` (three false novels that previously escaped as
*confident* novel). Of the 41 caught, **29 had no artifact mechanism flagged at all** under
the old rule — the indel-near signal is the only axis that explains them. **Zero** genuine
isoforms changed class.

So the wiring is a strictly-conservative correctness gain: it removes false novels the
existing mapping signals miss entirely, at no cost to genuine recall — the same
high-specificity ethos as the reference-bias rescue.

## Reproduce

```bash
WORK=/path/to/sqanti_sim/work ./run.sh   # needs the sqanti_sim work dir (FLAIR + SQANTI3 + aln.bam)
```

Inputs are public/simulated (SQANTI-SIM GENCODE v49 chr22, the same shared alignment as the
[sqanti_sim](../sqanti_sim) / [multicaller](../multicaller) protocols). The unit test
`rules: indel-near and soft-clip are mapping-artifact triggers` pins the rule logic.
Result envelope: [../results/bam_axis/metrics.json](../results/bam_axis/metrics.json).
