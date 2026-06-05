# PanIsoGuard on a real published ONT dataset (GM12878) — robustness

Applies PanIsoGuard's axes to a **real, public, whole-genome ONT direct-RNA** dataset —
ENCODE GM12878 (`ENCFF440ZML`, `minimap2 -ax splice -uf`) — and reports, honestly, what
each axis does on real data. The headline is a **chemistry-calibration** finding that
**empirically vindicates a conservative design choice**, plus the honest scope of the
reference-bias rescue on a reference-grade sample.

The inputs are the **285 caller-novel chains** from the
[wholegenome_multicaller](../wholegenome_multicaller) combine (IsoQuant + Bambu + ESPRESSO),
adjudicated with the real ONT BAM, plus a reference-bias scan of every novel junction
against NA12878's own genome (GM12878 = NA12878 = GIAB **HG001**).

## 1. The BAM mapping thresholds are chemistry-dependent (and soft-clip default-OFF is right)

Running the mapping axis (`--bam`) on the real ONT reads, fraction of the **191
BAM-evaluable** novel junctions each signal flags at the shipped 0.5 gate:

| BAM signal | fires on real ONT (frac > 0.5) | status |
|------------|:------------------------------:|--------|
| **soft-clip** | **187 / 191 (98 %)** | shipped **default-OFF** |
| **indel-near** | 175 / 191 (92 %) | default-on (HiFi-calibrated) |
| low-MAPQ | 11 / 191 (6 %) | default-on |
| supplementary | 4 / 191 (2 %) | default-on |

Two things this proves:

- **Soft-clip default-OFF is empirically vindicated.** ONT direct-RNA reads carry terminal
  soft-clips (adapter / poly-A) on **98 %** of junction-spanning reads. Had the soft-clip
  gate shipped **on** (it did, briefly, before the code review), it would have demoted
  nearly every novel call on real ONT data. This is exactly the false-positive the review
  flagged: the soft-clip feature is terminal-anywhere, not junction-proximal. It is now
  off by default (`max_softclip_frac = 1.01`); this real-data run is the receipt.
- **The HiFi-calibrated `indel_near` gate over-flags on ONT.** ONT's intrinsic indel error
  rate makes "indel adjacent to junction" pervasive, so the 0.5 gate (genuine max **0.19**
  on clean HiFi — see [bam_axis](../bam_axis)) fires on **92 %** of ONT junctions and barely
  separates single-caller (mean **0.925**) from ≥ 2-caller (mean **0.732**) novels
  (artifact rate **62 %** vs **58 %**). Meanwhile the *older* low-MAPQ / supplementary
  signals barely fire (6 % / 2 %). **Takeaway:** the BAM mapping axis is **opt-in (`--bam`)**
  and its thresholds **must be re-calibrated per chemistry**; the shipped defaults are
  HiFi-tuned. (Making `indel_near` rate-normalized for noisy chemistries is future work.)

## 2. The consensus axis is the chemistry-independent signal

The same 285 real novel chains separate cleanly at the **sequence level**, no read-chemistry
calibration required: **254 single-caller** vs **31 reproducible (≥ 2-caller)** — the 9.2×
decision impact quantified in [wholegenome_multicaller](../wholegenome_multicaller). This is
why PanIsoGuard leads with consensus + reference-bias (sequence-level) and treats the BAM
read-level axis as an optional, chemistry-sensitive corroborator.

## 3. Reference-bias rescue: 0 / 1364 (correctly) on a reference-grade sample

Scanning **all 1,364** caller-reported novel junctions against NA12878's own SNV-consensus
haplotype (HG001 GIAB v4.2.1) finds **0** reference-bias (`CREATED`) junctions. This is the
**expected, honest** result: GM12878 = NA12878 is the GIAB **reference-grade** sample (about
as concordant with GRCh38 as a human genome gets), so the rescue's base rate is ~0 — matching
the `pangenome_public` null. The rescue's value **rises with divergence** from the reference;
on personalized genomes it fires consistently and perfectly specifically (the
[giab_cohort_rescue](../giab_cohort_rescue): **137 across 4 individuals, 0 false**). A
divergent / non-reference long-read RNA sample with a matched personal genome is the data
that would show the rescue firing on real RNA calls; sourcing one is the open data-access item.

## Reproduce

```bash
WORK=/path/to/gm12878_mc ./run.sh   # needs the wholegenome_multicaller combine + the ONT BAM
                                     # + the HG001 SNV-consensus haplotypes (see giab_cohort_rescue)
```

All inputs are public (ENCODE GM12878 ONT dRNA `ENCFF440ZML`; GIAB HG001 v4.2.1; GENCODE v49);
no private data. Result envelope:
[../results/gm12878_realdata/metrics.json](../results/gm12878_realdata/metrics.json).
