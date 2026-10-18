# LRGASP WTC11: experimentally validated isoforms, real PacBio and ONT data

The first check of PanIsoGuard against **experimentally validated** isoforms on **real**
long-read data. LRGASP (Pardo-Palacios et al., *Nature Methods* 2024) amplified 89 WTC11
isoform targets by RT-PCR and sequenced the amplicons on ONT and PacBio; a target counted as
**supported** when at least one gapless (< 2 bp indel), high-identity amplicon read matched
its intron chain (Supplementary Data 15; 79 supported, 3 likely, 7 unsupported; 25 FSM, 36
NIC, 28 NNC by the LRGASP annotation). We join those targets to the WTC11 consolidated models
(their exact intron chains), run the public WTC11 reads through IsoQuant, SQANTI3 and
PanIsoGuard, and ask what class each validated chain received.

| step | PacBio cDNA (ENCSR507JOF, 3 reps) | ONT cDNA (ENCSR539ZXJ, rep 1) |
|---|---:|---:|
| reads aligned (minimap2 `splice:hq` / `splice`, primary) | 8.9 M | 17.1 M |
| IsoQuant 3.13 models (`pacbio_ccs` / `nanopore`) | 78,256 | 39,542 |
| SQANTI3 vs GENCODE v49: FSM / NIC / NNC / ISM | 59,704 / 7,971 / 4,478 / 2,383 | 34,828 / 1,495 / 2,073 / 816 |

Short reads: Illumina rep 1 (ENCSR673UKZ), STAR 2-pass on GENCODE v49 (279,044 junctions).
Adjudicated with `--bam`, `--platform hifi` / `ont`, with (`sr`) and without (`nosr`) the
SJ.tab. Everything public (ENCODE, GitHub LRGASP/supplementary_data, cgl.gi.ucsc.edu).

## Result

Of the 56 **supported NIC/NNC** targets, IsoQuant reconstructed the exact intron chain for
27 (PacBio) and 22 (ONT). Ten (PacBio) / six (ONT) of those are full splice matches in
GENCODE v49 (annotated since LRGASP's v38) and come out `KNOWN`. The rest are novel against
v49:

| supported novel targets, called exactly | PacBio + Illumina | PacBio, no short reads | ONT + Illumina | ONT, no short reads |
|---|---:|---:|---:|---:|
| n | 17 | 17 | 16 | 16 |
| `NOVEL_CONFIRMED` | 10 | 0 | 9 | 0 |
| `NOVEL_PROBABLE` | 6 (4 `chain_support`, 1 `chain_support;rt_switch_sign`, 1 `supported_with_rt_switch`) | **17** (`chain_support`, 2 of them `;rt_switch_sign`) | 5 (1 `chain_support`, 2 `chain_support;rt_switch_sign`, 2 `supported_with_rt_switch`) | **16** (`chain_support`, 4 of them `;rt_switch_sign`) |
| `UNCONFIRMED` | 1 | 0 | 2 | 0 |
| `ARTIFACT` | **0** | **0** | **0** | **0** |

- **No experimentally supported isoform is called an artifact**, on either platform, with
  or without short reads.
- With short reads, 16 of 17 (PacBio) and 14 of 16 (ONT) are confident. The `UNCONFIRMED`
  ones are CALB1 (NNC) and, on ONT, SAC3D1 (NIC): their novel junction has no short-read
  support in this one Illumina replicate (`unsupported_no_sign`), never `ARTIFACT`. Four
  validated isoforms (MTHFS, CENPK, USF1, RND2) carry SQANTI3's RT-switching sign; since
  ruleset 0.0.5 the sign is noted in the reason (`chain_support;rt_switch_sign`,
  `supported_with_rt_switch`) and does not hold them `UNCONFIRMED` (under 0.0.4 MTHFS and,
  on ONT, USF1 were held).
- **Without short reads, whole-chain read support does the work**: all 17 / 16 become
  `NOVEL_PROBABLE` through `chain_support` (2 / 4 with `;rt_switch_sign`), CALB1 and SAC3D1
  included: the long reads carry their chains, and only the short-read axis could not see
  them. That is the real-data counterpart of the SQANTI-SIM result in
  [chain_support](../chain_support).
- The called FSM targets (14 of 23 on PacBio, 18 of 23 on ONT) are all `KNOWN`, as are the
  10 / 6 NIC/NNC targets that GENCODE v49 has since annotated; the 2 / 1 called `likely`
  targets are `KNOWN` or `NOVEL_PROBABLE`.
- **Unsupported targets (7)**: IsoQuant reconstructed none of the six novel ones on either
  platform (so PanIsoGuard had nothing to judge), and the one FSM (ZBTB45, annotated) is
  `KNOWN` on PacBio. The validation set carries too few negatives to say anything about
  precision; this protocol is a recall-side check.

## The GENCODE-reviewed targets: do the negatives help?

Supplementary Data 15's second sheet holds 79 targets that GENCODE annotators reviewed
(`known` 27, `novel` 43, `rejected` 9; chains from LRGASP's manual-annotation BEDs, 76
resolvable), with **26 unsupported** by the amplicons. This is the only experimentally
confirmed negative set available, so it was the natural place to look for real-data
precision. It does not deliver it:

| GENCODE-reviewed, PacBio + Illumina | n | called exactly by IsoQuant | PanIsoGuard |
|---|---:|---:|---|
| supported | 46 | 16 | 14 `KNOWN`, 1 `NOVEL_PROBABLE`, 1 `UNCONFIRMED` |
| likely | 4 | 2 | 1 `KNOWN`, 1 `UNCONFIRMED` |
| unsupported, GENCODE `rejected` (the true negatives) | 7 | **0** | — |
| unsupported, GENCODE `known` / `novel` | 19 | 6 | 5 `KNOWN` (annotated transcripts whose amplicon failed), 1 `UNRESOLVED` (DLX4, NIC) |

ONT: 11 / 1 / 0 / 3 called, all `KNOWN`. The seven `rejected` models, the negatives that
matter, were not reconstructed by IsoQuant on either platform, so there was nothing for
PanIsoGuard to judge; the unsupported `known` ones are annotated transcripts (`KNOWN` by
pass-through, not adjudicated). The one unsupported novel chain the caller did produce
(DLX4+B) was held `UNRESOLVED` on both platforms, with and without short reads (fewer than
3 whole-chain reads). One case, in the right direction, is not a precision estimate.

**SQANTI3's RT-switching sign on validated isoforms.** Of the experimentally supported
targets IsoQuant reconstructed, `RTS_stage = TRUE` on 5 of 41 (PacBio; 3 FSM, 2 NIC) and
7 of 40 (ONT; 3 FSM, 3 NIC, 1 NNC): a 12–18 % false-sign rate on isoforms that exist.
PanIsoGuard consumes the sign as SQANTI3 reports it; this measured rate is why, since
ruleset 0.0.5, the sign only annotates a call that short reads or whole-chain reads support
(`;rt_switch_sign`) instead of holding it `UNCONFIRMED`, and why a sign alone never makes
an `ARTIFACT`.

**What it does not show.** 29 (PacBio) / 34 (ONT) supported novel targets were not
reconstructed by IsoQuant at all, and only one of them sits inside a longer model. They are
the rarely detected targets (median 1 LRGASP PacBio pipeline found them, against 11 for the
ones IsoQuant did call), the same abundance effect LRGASP reports. That is caller recall, not
PanIsoGuard's; the table above is conditional on the chain having been called. A validated
set that is mostly positives cannot estimate precision, and the LRGASP target selection
favoured isoforms detected by many pipelines, which read-support rules favour too.

## Reproduce

```bash
bash raw/dl.sh                         # ENCODE downloads (see the script; ~26 GB)
WORK=/path/to/pig_lrgasp bash run.sh all   # targets, align, STAR, IsoQuant, SQANTI3, adjudicate, score
```

`targets.py` joins Supplementary Data 15 (`supp/wtc11_validation_targets.tsv` and
`supp/gencode_validation_targets.tsv`, parsed from the xlsx) to `models/WTC11_consolidated.bed.gz`
and the two LRGASP manual-annotation BEDs (fetched by the `targets` stage); `score_lrgasp.py` matches our IsoQuant
chains to the targets exactly and writes `targets/per_target_verdicts.tsv` plus the envelope
[../results/lrgasp_wtc11/metrics.json](../results/lrgasp_wtc11/metrics.json). The STAR
binary must be the version that built the genome index.
