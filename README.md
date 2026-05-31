# PanIsoGuard

**Caller-agnostic adjudication of long-read RNA-seq novel isoforms.**

[![CI](https://github.com/jibeomko/PanIsoGuard/actions/workflows/ci.yml/badge.svg)](https://github.com/jibeomko/PanIsoGuard/actions/workflows/ci.yml)
[![License: MIT AND BSL-1.0](https://img.shields.io/badge/license-MIT%20AND%20BSL--1.0-blue.svg)](LICENSE)

PanIsoGuard is a post-processing / decision layer that ingests the *novel* isoform
calls produced by long-read isoform callers (FLAIR, IsoQuant, Bambu, ESPRESSO,
TALON, …) and/or [SQANTI3](https://github.com/ConesaLab/SQANTI3) output, together
with a BAM and **optional** evidence inputs (short-read `SJ.tab`, personalized
haplotype FASTA), and re-classifies each novel call into a confidence class with a
**machine-readable mechanistic attribution** and a **provenance / circularity flag**.

It does not recompute SQANTI3's QC features and does not claim to beat the SQANTI3
random-forest filter or FLAIR2. Its contribution is the **adjudication logic** —
how orthogonal evidence axes are integrated into a transparent, auditable verdict
(thresholds live in a runtime [`config/rules.default.toml`](config/rules.default.toml),
and each verdict carries a `rule_trace`; see [docs/decision_engine.md](docs/decision_engine.md)).

## At a glance

![PanIsoGuard overview: long-read novel isoform calls (real or artifact?) are checked against four kinds of evidence — SQANTI QC priors, short-read junctions, long-read mapping, and variants/reference bias — and sorted into plain-language verdicts: real novel, reference-bias rescued, uncertain (held), or artifact.](docs/figures/overview.png)

<sub>Consumes caller + SQANTI3 output (does not replace them); adds an orthogonal,
auditable verdict layer, each call carrying a machine-readable `rule_trace`. The four
output classes shown are a simplified grouping — the full set is listed
[below](#confidence-classes). Vector source:
[`docs/figures/overview.svg`](docs/figures/overview.svg).</sub>

> **Status: alpha.** Four evidence axes (SQANTI priors, short-read junctions, BAM
> read-level mapping, variant/reference-bias) and the `adjudicate` / `benchmark` /
> `ablate` / `combine` subcommands are implemented and tested. Rule thresholds are
> conservative defaults that are **not yet calibrated** against simulated truth
> (SQANTI-SIM); treat the confidence classes as orthogonal evidence integration,
> not a calibrated probability.

## Subcommands

| Command | Purpose |
|---------|---------|
| `adjudicate` | classify novel isoforms → confidence class + mechanistic attribution + provenance (3 output files) |
| `benchmark`  | non-redundancy vs SQANTI3 on novel isoforms (2×2, Jaccard, McNemar) |
| `ablate`     | per-evidence-axis class-change (which axis drives which calls) |
| `combine`    | integrate multiple callers' isoforms by intron-chain fingerprint → caller-support matrix |
| `version`    | version, linked htslib, compiled-in capabilities |

`panisoguard <command> --help` for options. Input formats: [docs/input_formats.md](docs/input_formats.md).

## Usage

`adjudicate` is the main entry point. The **only** hard requirements are a SQANTI3
classification and an output prefix — every evidence input below is optional and
simply switches on another axis (see [Evidence tiers](#evidence-tiers)).

**Baseline** — SQANTI priors only (no short/long-read evidence; novel calls are
flagged or held, never positively confirmed):

```bash
panisoguard adjudicate \
  --classification sample_classification.txt \
  --isoforms-bed   flair.isoforms.bed \
  --ref-gtf        gencode.v49.annotation.gtf \
  --out-prefix     out/sample
```

**Recommended** — add short-read junctions (`--sj-tab`) and the long-read BAM
(`--bam`), the two axes that let a novel isoform be confirmed or rejected on evidence:

```bash
panisoguard adjudicate \
  --classification sample_classification.txt \
  --isoforms-bed   flair.isoforms.bed \
  --ref-gtf        gencode.v49.annotation.gtf \
  --sj-tab         star.SJ.out.tab \
  --bam            aligned.sorted.bam \
  --reference      GRCh38.fa \
  --out-prefix     out/sample
```

**Reference-bias rescue** — add a personalized haplotype FASTA. Provenance gates the
circularity firewall: `wgs`/`external` may promote to a rescue verdict, while
`rna_derived`/`unknown` are held as `AMBIGUOUS`:

```bash
panisoguard adjudicate \
  --classification sample_classification.txt \
  --isoforms-bed   flair.isoforms.bed \
  --ref-gtf        gencode.v49.annotation.gtf \
  --reference      GRCh38.fa \
  --reference-haplotype sample.hap1.fa \
  --reference-haplotype sample.hap2.fa \
  --haplotype-provenance wgs \
  --out-prefix     out/sample
```

**Combine several callers first** (optional) — merge isoforms by intron-chain
fingerprint into a caller-support matrix, then feed the union to `adjudicate`:

```bash
panisoguard combine \
  --gtf flair:flair.gtf \
  --gtf isoquant:isoquant.gtf \
  --gtf bambu:bambu.gtf \
  --ref-gtf gencode.v49.annotation.gtf \
  --out caller_support_matrix.tsv
```

### Outputs

`adjudicate` writes three files at `<out-prefix>`:

| File | Contents |
|------|----------|
| `<prefix>.adjudicated.tsv`   | one row per isoform — confidence class, primary mechanism, novel-junction support counts |
| `<prefix>.attribution.jsonl` | per-isoform `rule_trace`: every rule that fired, in order (fully auditable) |
| `<prefix>.provenance.log`    | which axes were active + circularity status of the run |

### Reviewer-facing analyses

```bash
# non-redundancy vs the SQANTI3 filter on novel isoforms (2x2, Jaccard, McNemar)
panisoguard benchmark [adjudicate options] --out bench/

# per-axis contribution: which calls change when an axis is removed
panisoguard ablate    [adjudicate options] --axes short_read,mapping,variant --out abl/
```

## Confidence classes

`HIGH_CONF_KNOWN` · `HIGH_CONF_NOVEL` · `MEDIUM_CONF_NOVEL` · `LOW_CONF_PARTIAL` ·
`PAN_REF_RESCUED_FALSE_NOVEL` · `AMBIGUOUS` · `ARTIFACT` — emitted as a
deterministic projection of a 2-axis evidence grid (novelty-support ×
artifact-mechanism).

## Evidence tiers

| Tier | Input | Required? | Mechanism |
|------|-------|-----------|-----------|
| 0 | SQANTI3 classification (priors) + STAR `SJ.tab` | recommended | short-read junction corroboration |
| 1 | BAM (HiFi/ONT) | recommended | read-level mapping / chimera / soft-clip features |
| 2 | personalized haplotype FASTA (`--reference-haplotype`) | optional | variant-created/destroyed splice-site motif |

## Build

Requires a C++17 compiler, CMake ≥ 3.20, and **htslib ≥ 1.18**.

```bash
git clone <repo> && cd PanIsoGuard
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # unit suite
./build/panisoguard --version
```

htslib is discovered from `$CONDA_PREFIX`; override with `-DCMAKE_PREFIX_PATH=/prefix`.

## Install (conda)

A bioconda recipe is provided under [`recipes/bioconda/`](recipes/bioconda/). Once
released:

```bash
conda install -c bioconda -c conda-forge panisoguard
```

## Runtime & memory

Single-threaded, on one cohort's merged set (188,912 isoforms; GRCh38 + GENCODE v49):

| Step | Wall | Peak RAM |
|------|------|----------|
| reference catalog (GENCODE v49) | ~3.2 s | ~0.34 GB |
| `adjudicate` (SQANTI priors + short-read SJ) | ~5 s | ~0.55 GB |
| + variant axis (faidx motif) | ~40 s | ~0.6 GB |
| + BAM mapping axis (4.5 GB BAM) | ~1.6 min | ~0.6 GB |

The upstream callers + alignment dominate end-to-end time; PanIsoGuard's own
adjudication is the fast tail.

## Orchestration (optional)

[`workflow/`](workflow/) provides a Snakemake pipeline that runs the upstream
callers in parallel over a single shared alignment and pipes into PanIsoGuard
(`combine` + `adjudicate`). It is a thin convenience wrapper, not the core tool.

## Validation

See [docs/validation.md](docs/validation.md) for what is verified and the
truth-based validation plan (SQANTI-SIM, HG002/HPRC, LRGASP).

## License

PanIsoGuard is MIT-licensed — see [LICENSE](LICENSE).

It bundles three third-party single-header components under `thirdparty/`, with
their license texts included: **cgranges** (`IITree.h`, MIT), **toml++** (MIT), and
**Catch2** (Boost Software License 1.0, test-only). See [THIRDPARTY.txt](THIRDPARTY.txt)
for attribution. The effective combined license of the redistributed source is
**MIT AND BSL-1.0**. htslib is a dynamically-linked dependency, not vendored.
