# PanIsoGuard

**Evaluates novel isoforms from long-read RNA-seq and explains each verdict.**

[![CI](https://github.com/jibeomko/PanIsoGuard/actions/workflows/ci.yml/badge.svg)](https://github.com/jibeomko/PanIsoGuard/actions/workflows/ci.yml)
[![License: MIT AND BSL-1.0](https://img.shields.io/badge/license-MIT%20AND%20BSL--1.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23131084.svg)](https://doi.org/10.5281/zenodo.23131084)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/bioconda/panisoguard)
[![bioconda version](https://img.shields.io/conda/vn/bioconda/panisoguard.svg)](https://anaconda.org/bioconda/panisoguard)
[![platforms](https://img.shields.io/conda/pn/bioconda/panisoguard.svg)](https://anaconda.org/bioconda/panisoguard/files)

Long-read RNA-seq callers report novel isoforms that may reflect genuine transcripts,
sequencing or alignment artifacts, or differences from the reference genome. PanIsoGuard
evaluates these candidates using available read evidence and QC information after isoform
calling.

It integrates short-read junction support, long-read alignments and SQANTI3 metrics, with
optional caller agreement, personalized genomes and pangenome evidence. Each isoform receives
a verdict, the evidence behind it and a rule trace; candidates with insufficient or conflicting
evidence are left for review. PanIsoGuard uses existing caller outputs and alignments. It also
provides a built-in classification mode for runs without SQANTI3.

> **Status: alpha.** Default thresholds were set on simulated data and checked on real
> PacBio and ONT data. See [benchmarks and limitations](#benchmark--limitations).

## Key features

- **Combine evidence:** assess junction support, whole-chain reads, QC signals and caller
  agreement in one run.
- **Assess possible reference bias:** use personalized genomes or pangenome junctions to
  investigate apparent novelty, with checks on the declared source of the evidence.
- **Explain each verdict:** report the supporting evidence and decision rules, and leave
  candidates with insufficient evidence for review.

<a id="install"></a>

## Quick start

Install from [Bioconda](https://anaconda.org/bioconda/panisoguard), then run a bundled example:

```bash
conda install -c conda-forge -c bioconda panisoguard
panisoguard example tiny
bash tiny/run.sh
```

The second command copies the example into `./tiny`; the third runs it and checks its results
against the bundled expected output. It finishes with `example output matches expected files`.
Results are in `tiny/output/`.

For a separate conda environment, Docker, source builds and other bundled examples, see
[installation](docs/installation.md).

<a id="at-a-glance"></a>

## How it works

![Five steps: read caller output and available evidence, identify novel junctions, gather evidence, apply decision rules, and report a verdict with its reasons. The example isoform has a novel junction supported by short reads and receives a confirmed verdict.](docs/figures/overview.png)

The figure shows the SQANTI3-based workflow. With `--no-sqanti3`, PanIsoGuard computes
structural categories from the isoforms and reference annotation instead. Follow one toy gene
through its inputs and verdicts in the [worked example](examples/worked_locus/README.md).

<a id="in-one-screen"></a>
<a id="usage"></a>
<a id="outputs"></a>

## Inputs & outputs

For your own data, provide caller isoforms in GTF or BED12 format, the reference annotation,
and either a SQANTI3 classification table or `--no-sqanti3`. When using SQANTI3, isoform IDs
must match its classification table and the reference GTF must match the annotation used in
that run. Short-read junctions and long-read alignments are optional evidence inputs.

This example uses SQANTI3 output; replace the file names with your own paths:

```bash
panisoguard adjudicate \
    --classification classification.txt \
    --isoforms-gtf isoforms.gtf \
    --ref-gtf reference.gtf \
    --sj-tab SJ.out.tab \
    --bam aligned.bam \
    --out-prefix out/sample
```

`--sj-tab` and `--bam` are optional. For a run without SQANTI3, replace `--classification
classification.txt` with `--no-sqanti3`; add `--reference genome.fa` to evaluate splice motifs
and downstream A content. SQANTI3's RT-switch flag is not evaluated in that mode.

<a id="confidence-classes"></a>

The main result, `out/sample.adjudicated.tsv`, contains one row per isoform:

| `verdict` | Interpretation |
|---|---|
| `known` | Matches an annotated transcript. |
| `known_fragment` | A fragment of an annotated transcript supported by whole-chain reads. |
| `confirmed` | Meets the evidence rules for a supported novel isoform; the confidence class distinguishes junction support from whole-chain or caller support. |
| `reference_bias` | Genomic evidence supports a possible reference-bias explanation for the novel junctions. |
| `unconfirmed` | Evidence is insufficient, unavailable or conflicting; review the candidate. |
| `artifact` | Unsupported novel junctions plus an artifact signal meet the rejection rules. |
| `not_adjudicated` | The structural category is outside the current adjudication scope. |

`reference_bias` does not independently confirm the full isoform. In particular, a junction
found in a population pangenome does not establish that the sample carries or expresses it.

The run also writes `.attribution.jsonl` (evidence and rule traces), `.junctions.tsv`
(per-junction evidence) and `.provenance.log` (inputs, versions and settings).
Add `--write-kept kept.gtf` to export a filtered annotation: by default it retains `known`,
`known_fragment`, `confirmed` and `reference_bias`. Change this policy with `--keep-verdicts`.
See the [output and confidence-class reference](docs/usage.md#outputs) for details.

<a id="how-it-compares"></a>
<a id="validation"></a>

## Benchmark & limitations

On SQANTI-SIM chr22, 870 isoforms that do not fully match an annotated transcript (730 genuine,
140 false) were evaluated after FLAIR and SQANTI3. The following comparison uses 1 million
simulated Illumina read pairs; PanIsoGuard also receives the long-read BAM:

| Method | Precision | Recall | AUPRC |
|---|---:|---:|---:|
| SQANTI3 rules filter | 0.918 | 0.890 | 0.909 |
| Every junction ≥ 3 short reads | 0.976 | 0.945 | 0.968 |
| PanIsoGuard + BAM | 0.993 | 0.752 | 0.975 |

The default operating point favors precision over recall. AUPRC is similar to the simple
short-read rule within bootstrap uncertainty. In this setting, 641 of 870 candidates
(73.7 %) receive a supported-novel or artifact class, with 1.4 % errors within that group.
Another 164 are known fragments, and 65 are held for review. These are simulation results,
not an estimate of accuracy on every sample.

- **Without short reads:** with the BAM, AUPRC is 0.953 versus 0.909 for SQANTI3's default
  rules filter on chr22, and 0.958 for ranking by the caller's own read count. Separate chr20
  and chr19 evaluations support the same ordering against SQANTI3's rules. Without
  read support, caller agreement or usable genomic evidence, novel candidates remain
  unconfirmed or unresolved.
- **Real-data checks:** in LRGASP WTC11, none of the experimentally supported novel targets
  reconstructed by IsoQuant (17 PacBio, 16 ONT) is called an artifact. This subset checks
  retention of supported transcripts; it does not establish precision on all novel calls.
- **Reference-bias scope:** judged on two individuals' own assemblies, 96 % of the junctions
  attributed to the sample's genome (a variant at the splice site, or a deletion equal to the
  junction) are confirmed, about half of those the assembly accounts for. With short-read
  variant calls of the same people, 93-94 % are confirmed and about a fifth to a quarter are
  found, because short reads miss most of the deletions. 27 of 28 testable heterozygous
  junctions agree with the haplotype carrying the canonical motif. The pangenome rescue
  benchmark is limited to HPRC chr22. A genomic explanation does not prove the full transcript
  is correct.
- **Confidence classes:** these are evidence labels, not probabilities or an ordered score.
  Results depend on the supplied evidence and its coverage.

See [full comparisons, caveats and runtime measurements](docs/benchmarks.md),
[validation protocols](docs/validation.md) and the
[tracked results](benchmark/results/README.md) for denominators, confidence intervals and
reproduction commands.

<a id="documentation"></a>

## Advanced usage & documentation

| Topic | Guide |
|---|---|
| Installation | [Bioconda, Docker and source builds](docs/installation.md) |
| Advanced usage | [Personalized genomes, pangenomes and multiple callers](docs/usage.md) |
| Inputs and options | [Input formats](docs/input_formats.md), [default thresholds](config/rules.default.toml); `panisoguard adjudicate --help` lists CLI options |
| Outputs and evidence | [Output files](docs/usage.md#outputs), [confidence classes](docs/usage.md#confidence-classes), [evidence inputs](docs/usage.md#evidence-it-uses) |
| Decision rules | [Algorithm](docs/algorithm.md), [decision grid and provenance checks](docs/decision_engine.md) |
| Related tools | [SQANTI3](docs/relationship_to_sqanti3.md), [gffcompare, TAMA and other callers](docs/relationship_to_merge_tools.md) |
| Benchmarks | [Comparison and limitations](docs/benchmarks.md), [validation](docs/validation.md) |
| Worked example | [One toy gene from inputs to verdicts](examples/worked_locus/README.md) |
| Development | [Architecture and repository layout](docs/architecture.md), [module interfaces](docs/function_io.md), [design notes](docs/design_notes.md) |
| Releases | [Changelog](CHANGELOG.md), [release process](docs/releasing.md) |

## Citation

Each release is archived on Zenodo. To cite PanIsoGuard, use the DOI for all versions,
[10.5281/zenodo.23131084](https://doi.org/10.5281/zenodo.23131084), which always resolves to
the latest release; each release also has its own DOI on that page. GitHub's "Cite this
repository" button gives the same record from [`CITATION.cff`](CITATION.cff).

## License

PanIsoGuard is MIT-licensed; see [LICENSE](LICENSE).

It includes two small third-party header libraries in `thirdparty/`, with their license
texts: **toml++** (MIT) and **Catch2** (Boost Software License
1.0, used only for tests). See [THIRDPARTY.txt](THIRDPARTY.txt). Taken together, the source as
distributed is under **MIT AND BSL-1.0**. htslib is linked as a separate dependency, not
included.
