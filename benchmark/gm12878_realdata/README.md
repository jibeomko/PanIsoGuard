# PanIsoGuard on a real published ONT dataset (GM12878) — robustness

Applies PanIsoGuard's axes to a **real, public, whole-genome ONT direct-RNA** dataset —
ENCODE GM12878 (`ENCFF440ZML`, `minimap2 -ax splice -uf`) — and reports, honestly, what
each axis does on real data. The headline is a **chemistry-calibration** finding that
**empirically vindicates a conservative design choice**, plus the honest scope of the
reference-bias rescue on a reference-grade sample.

Section 4 (0.0.7) adds the whole-chain read support on the same reads. The inputs are the
**285 caller-novel chains** from the
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
  (the mapping sign fires on 17 / 31 = 55 % of ≥ 2-caller vs 112 / 254 = 44 % of single-caller chains). Meanwhile the *older* low-MAPQ / supplementary
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
on personalized genomes it fires at a similar rate per genome with no false rescue in those scans (the
[giab_cohort_rescue](../giab_cohort_rescue): **137 across 4 individuals, 0 false**). A
divergent / non-reference long-read RNA sample with a matched personal genome is the data
that would show the rescue firing on real RNA calls; sourcing one is the open data-access item.

## 4. Whole-chain read support works on real ONT with exact coordinates (0.0.7)

The 0.0.7 chain axis counts reads whose alignment carries the isoform's entire intron chain
at **exact** coordinates. ONT reads are expected to wobble at splice sites, so this is the
test of whether the axis is usable on ONT at all. On the 285 real novel chains
(`chain_wobble.py`, reads at each chain's first donor site, median 17):

| per-site tolerance | chains with ≥ 1 full-chain read | chains with ≥ 3 |
|---|---:|---:|
| exact (what PanIsoGuard counts) | 268 / 285 | **155 / 285** |
| ± 1 bp | 269 | 155 |
| ± 2 bp | 271 | 155 |
| ± 5 bp | 281 | 175 |
| ± 10 bp | 282 | 188 |

Among reads that match a chain within 10 bp, **94.5 %** of the splice sites are exact; 1–2 bp
offsets are essentially absent (0.1 %), the rest sit at 3–7 bp. These are caller chains built
from the same splice-aware alignment (`minimap2 -ax splice -uf`), so the aligner has already
placed the sites consistently; a tolerance of ± 1–2 bp would change nothing, ± 5 bp would add
20 chains at the cost of accepting neighbouring alternative sites. **No tolerance is shipped.**
Exact full-chain support is more frequent among the ≥ 2-caller chains (22 / 31, 71 %) than
the single-caller ones (133 / 254, 52 %), the same direction as the consensus signal.

What the two presets do with it (no short reads exist for this sample, so every call is
`UNKNOWN` support):

| preset | UNCONFIRMED (`mapping_without_short_reads`) | NOVEL_PROBABLE | UNRESOLVED |
|---|---:|---:|---:|
| `--platform hifi` (default; indel-near gate on) | 129 | 53 | 103 |
| `--platform ont` (indel-near gate off) | 9 | **148** | 128 |

Under the HiFi defaults the over-firing indel-near gate (section 1) holds most chain-supported
isoforms `UNCONFIRMED` (before ruleset 0.0.5 they were `ARTIFACT`; this dataset is one reason
a mapping sign alone no longer rejects a call); under the ONT preset they become
`NOVEL_PROBABLE`. There is no truth
for these chains, so this section is about whether the axis functions on ONT (it does), not
about accuracy.

## Reproduce

```bash
WORK=/path/to/gm12878_mc ./run.sh   # needs the wholegenome_multicaller combine + the ONT BAM
# the data dir of wholegenome_multicaller/run.sh: gm.main.bam + matrix.tsv; override with BAM=... MATRIX=... when named otherwise
                                     # + the HG001 SNV-consensus haplotypes (see giab_cohort_rescue)
                                     # PY=<python with pysam> for chain_wobble.py
```

All inputs are public (ENCODE GM12878 ONT dRNA `ENCFF440ZML`; GIAB HG001 v4.2.1; GENCODE v49);
no private data. Result envelope:
[../results/gm12878_realdata/metrics.json](../results/gm12878_realdata/metrics.json).
