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
| `NOVEL_PROBABLE` | 5 (4 `chain_support`, 1 `supported_with_rt_switch`) | **15** (`chain_support`) | 3 | **12** (`chain_support`) |
| `UNCONFIRMED` | 2 | 2 | 4 | 4 |
| `ARTIFACT` | **0** | **0** | **0** | **0** |

- **No experimentally supported isoform is called an artifact**, on either platform, with
  or without short reads.
- With short reads, 15 of 17 (PacBio) and 12 of 16 (ONT) are confident. The `UNCONFIRMED`
  ones are CALB1 (NNC; its novel junction has no short-read coverage in this one Illumina
  replicate) and, on both platforms, MTHFS and (ONT) USF1 and SAC3D1, where SQANTI3 flags
  an RT-switching sign or short reads are absent at the junction. The validated status of
  those four says the RTS sign is a false sign at least sometimes; PanIsoGuard keeps them as
  `UNCONFIRMED` (a sign plus no corroboration), never `ARTIFACT`.
- **Without short reads, whole-chain read support does the work**: the same 15 / 12 become
  `NOVEL_PROBABLE` through `chain_support`, the 2 / 4 with an RTS sign `UNCONFIRMED`
  (`chain_support_with_rt_switch`). That is the real-data counterpart of the SQANTI-SIM
  result in [chain_support](../chain_support).
- All 24 (PacBio) / 24 (ONT) called FSM targets are `KNOWN`; the 2 / 1 called `likely`
  targets are `KNOWN` or `NOVEL_PROBABLE`.
- **Unsupported targets (7)**: IsoQuant reconstructed none of the six novel ones on either
  platform (so PanIsoGuard had nothing to judge), and the one FSM (ZBTB45, annotated) is
  `KNOWN` on PacBio. The validation set carries too few negatives to say anything about
  precision; this protocol is a recall-side check.

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

`targets.py` joins Supplementary Data 15 (`supp/wtc11_validation_targets.tsv`, parsed from
the xlsx) to `models/WTC11_consolidated.bed.gz`; `score_lrgasp.py` matches our IsoQuant
chains to the targets exactly and writes `targets/per_target_verdicts.tsv` plus the envelope
[../results/lrgasp_wtc11/metrics.json](../results/lrgasp_wtc11/metrics.json). The STAR
binary must be the version that built the genome index.
