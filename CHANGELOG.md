# Changelog

All notable, user-facing changes. Format follows [Keep a Changelog](https://keepachangelog.com);
this project uses [Semantic Versioning](https://semver.org) (pre-1.0: minor/patch only).

## [0.0.5] - 2026-10-04

### Added
- **Study notes** ([notes/](notes/README.md), in Korean): nine notes that follow one toy
  gene from the input files to the verdict (coordinates and novel junctions, short-read
  support, artifact mechanisms, the projection, reference-bias rescue, caller consensus,
  evaluation, and an end-to-end recomputation). The toy inputs are built with real
  SQANTI3 6.0.1 and STAR runs (`notes/data/build_toy.sh`). `notes/check_notes.py`
  recomputes all 81 toy verdicts with the Python standard library and compares them with
  the binary; it runs as the CTest `integration_study_notes`.
- **Copy-paste Quick start.** One block clones, creates a conda env, builds and runs both
  bundled examples; a Docker alternative follows. Both blocks were run verbatim from the README
  in a stock Debian 13 Miniconda container with no system compiler or `ar`. CI now runs the
  Docker route (`container` job) and both examples (CTest `example_tiny`,
  `example_multi_caller`).
- **Non-circular check of the reference-bias rescue**
  ([benchmark/refbias_phasing](benchmark/refbias_phasing)). For heterozygous reference-bias
  junctions in HG03516 and HG02717, the long reads that use the junction are phased with the
  trio-phased HiFi assembly SNVs. 27 of 28 testable junctions pass: 572 of 575 phased junction
  reads come from the haplotype on which the junction is canonical, while both haplotypes are
  expressed at those loci. The one failure (HG03516 chr3:184709997) has junction reads from the
  non-canonical haplotype, all with a 2-bp insertion next to the donor; the rescue is checked
  before the mapping mechanism, so it is rescued anyway. The 34 homozygous junctions cannot be
  phased.
- `adjudicate` warns when SQANTI3 isoform ids do not join: how many are missing from the
  caller isoform file (they fall to `caller_chain_absent`) and how many novel ones are
  missing from the `--caller-support` matrix (no consensus). Both used to be silent.
- `provenance.log` gains a `thresholds` line with every effective rule threshold, so a run
  can be reproduced even if its config file later changes.
- `CITATION.cff` (author, ORCID, license, version): GitHub shows a "Cite this repository"
  button, and Zenodo archives each release with a DOI.

### Changed
- `ruleset_version` is now `builtin-0.0.2` / `default-0.0.2`. The 0.0.4 `PARTIAL` change
  had kept `0.0.1`, so 0.0.3 and 0.0.4 reported the same name for different projections.
  The name is to be bumped with any change to a default threshold or to the projection.
  Verdicts are identical to 0.0.4: no threshold or projection changed in this release.

### Fixed
- **Build failed without system binutils.** conda-forge compilers >= 2.0 install an
  unprefixed `c++` but only a prefixed `ar`, so where no system `/usr/bin/ar` exists (a clean
  container, and this project's own Docker image) CMake found no archiver and linking
  `libpanisoguard_core.a` failed; the Docker build on `main` was broken by it. CI never saw it
  because its runners have system binutils. The core is now an OBJECT library, which needs no
  archiver.
- README / Dockerfile container commands: the command must start with `panisoguard` (a bare
  `adjudicate` is "not found"), and `-u "$(id -u):$(id -g)"` is needed to write outputs into
  the mounted folder.
- `examples/tiny/expected/sample.provenance.log` still said `tool_version 0.0.3`, so the
  README quick example failed its own diff; nothing ran it. The diff now ignores the
  version line, and the example runs as the CTest `example_tiny`.
- `config/rules.default.toml` header still quoted "AUPRC 0.970 vs 0.831 baseline" and an
  unmet GATE-1.
- `docs/decision_engine.md`: the indel-near fraction does drive the mapping mechanism
  (since 0.0.4); `rule_trace` records only the first-matching mechanism, not every flag;
  `ablate --axes consensus` (not `--without`); the `UNKNOWN` reason is
  `short_read/catalog_axis_absent`. Matching stale comment in `src/core/adjudicator.cpp`.
- `docs/validation.md`: the proposed independent check of the reference-bias rescue
  ("do junction reads carry the alt allele?") cannot work, because the motif-creating base
  is intronic and never in those reads; replaced with phasing / carrier comparison.
- README: the bioconda package is on the channel (it said the recipe was still in review).

### CI / tooling
- Bioconda recipe passes `${CMAKE_ARGS}` to CMake (bioconda review suggestion).

## [0.0.4] - 2026-09-27

The flagship change is the **multi-caller consensus axis**, plus a PDF report tool, a
container, a large round of validation, a head-to-head against the SQANTI3 filter, and the
`PARTIAL` projection fix it exposed. PanIsoGuard is now positioned as a traceable-verdict
tool, not a more accurate filter (see README, "How it compares").

### Added
- **BAM mapping axis — indel-near now feeds the verdict.** The read-level
  `bam_frac_indel_near` and `bam_frac_softclip` fractions were measured at every junction
  but never acted on (only low-MAPQ / supplementary did). **`indel_near`** is now a
  `mapping_or_repeat` trigger (`max_indel_near_frac`, default 0.5). Validated on chr22
  SQANTI-SIM: it separates genuine novel junctions (max 0.190) from false (mean 0.350)
  and, gated at 0.5, catches **41/124 false at 100 % precision** with **0 genuine loss**.
  (Its earlier specificity lift 0.946 → 0.966 came from 3 partially-supported false novels
  that the PARTIAL fix below now handles without the BAM; on this set it now only sharpens
  attribution, 24 `LOW_CONF_PARTIAL → ARTIFACT`.)
  **`softclip`** is wired but **disabled by default** (`max_softclip_frac = 1.01`): the
  feature is a *terminal* soft-clip (not junction-proximal), so a default-on gate would
  risk demoting genuine novels on noisy reads — opt-in only. New
  [benchmark/bam_axis](benchmark/bam_axis) + unit test.
- **Multi-caller consensus axis.** `adjudicate --caller-support <combine matrix>` feeds
  cross-caller agreement into the verdict: when short-read support is `UNKNOWN`, a novel
  chain recovered by ≥ `consensus_min_callers` (default 2) callers is promoted
  `AMBIGUOUS → MEDIUM_CONF_NOVEL` (never `HIGH`; a mapping artifact still preempts).
  Ablate with `ablate --axes consensus`; emitted as `consensus_evaluable` / `n_callers`
  in `attribution.jsonl`.
- **`panisoguard-report`** — an optional Python companion that renders a publication-grade,
  SQANTI3-style multi-page **PDF report** from adjudication output (`python/`). Console
  entry point + `pyproject.toml`; matplotlib is its only dependency.
- **Container** — a multi-stage `Dockerfile` that builds the C++ binary from source and
  bundles `panisoguard-report` (works without waiting on the conda release).
- **Quickstarts** — `examples/multi_caller/` (self-contained, offline, shows the consensus
  gradient) and `examples/report/` (an 8-page showcase PDF + generator).
- **Validation benchmarks** (all git-tracked, schema-checked):
  `pangenome` (GATE-1 on the real HPRC v1.1 chr22 graph),
  `pangenome_public` (GM12878 + K562 whole-genome specificity),
  `multicaller` (5 callers, truth-scored; PR curve over the consensus threshold) +
  `sirv_multicaller` (the consensus generalized to a 2nd reference, SIRV-Set4),
  `wholegenome_multicaller` (real GM12878 + the union-vs-consensus decision impact),
  `merge_comparison` (head-to-head vs gffcompare / TAMA), `hg002_wholegenome`
  (variant-rescue yield), `giab_cohort_rescue` (rescue across a 4-individual GIAB
  cohort: 137/137, 0 false), `bam_axis` (wiring the indel-near / soft-clip read
  signals: 41/124 false caught, 0 genuine loss), and
  `gm12878_realdata` (the BAM axis on a **real** ONT dataset: soft-clip default-off
  empirically vindicated — 98 % fire rate on ONT — and the mapping thresholds shown to
  be chemistry-dependent; reference-bias rescue 0/1364 on the reference-grade sample), and
  `refbias_cohort` + `hg03516_refbias` (the **positive real-data application, replicated**:
  the reference-bias rescue fires on real PacBio Iso-Seq from **two** divergent West African
  individuals — HG03516 (ESN) + HG02717 (GWD), HPRC R2 — **86 reference-bias novel junctions
  (77 novel vs GENCODE), 86/86 rescued, 0 false, firewall holds all 86**, perfectly specific
  across both genomes, vs 0 on reference-grade GM12878; independently adversarially audited,
  caveats stated — the differentiator firing on real divergent-genome RNA).

### Fixed
- **`PARTIAL` short-read support no longer promotes to `MEDIUM_CONF_NOVEL`.** An isoform
  with some but not all novel junctions corroborated was called a confident novel. The
  head-to-head against the SQANTI3 filter showed these were exactly PanIsoGuard's extra
  false positives. `PARTIAL` now always projects to `LOW_CONF_PARTIAL`. SQANTI-SIM:
  precision 0.982 → 0.989, specificity 0.946 → 0.966, AUPRC 0.970 → 0.971, recall
  unchanged.
- **`config/rules.default.toml` silently disabled the consensus axis** — `consensus_min_callers`
  was left as a `999` "reserved" sentinel after the axis was wired, so passing
  `--config config/rules.default.toml` turned consensus off. Set to the built-in default (2)
  and guarded by a new `integration_config_equivalence` regression test.
- **Report page-2 projection grid** mis-coloured rescue-mechanism cells
  (`variant_created` / `population_known` bypass the grid); restricted to the artifact
  mechanisms.
- Hardened junction aggregation and provenance.

### Changed
- **Repositioned as a traceable-verdict tool, not a better filter.** New
  [benchmark/sqanti3_filter_h2h](benchmark/sqanti3_filter_h2h) compares PanIsoGuard with
  the SQANTI3 rules / ML filter and a one-line short-read rule on SQANTI-SIM truth,
  using oracle and simulated Illumina junctions. PanIsoGuard ties the one-line rule on
  AUPRC; only the BAM axis adds a small gain. Recall is much lower because of
  abstention. With no short reads it is worse than the SQANTI3 rules filter. The README
  now says this up front. The "AUPRC 0.970 vs 0.831 baseline" claim is dropped: 0.831 was
  the positive rate, not a method. `docs/validation.md` now states that the
  reference-bias benchmarks use the rescue rule's own criterion as truth.
- **Honest positioning** (after head-to-head benchmarking): `combine` is a clean
  re-implementation of `gffcompare -i` (not a novel merge), the consensus precision is
  matcher-robust, and the reference-bias rescue is a **high-specificity, high-sensitivity
  guardrail with a low base rate** (~1 in 18k introns genome-wide) — a correctness
  safeguard, not a high-yield discovery engine. New
  [docs/relationship_to_merge_tools.md](docs/relationship_to_merge_tools.md).
- **SQANTI3 attribution** — independence disclaimer + citation (Pardo-Palacios *et al.*,
  *Nat Methods* 2024) + the "rescue" name-clash disambiguation.
- README: a Table of Contents + a "When to use PanIsoGuard" section.

### CI / tooling
- `scripts/check_version_sync.py` + CI lint: the version must be identical across
  `CMakeLists.txt`, `python/pyproject.toml`, and `recipes/bioconda/meta.yaml`.
- `integration_config_equivalence` CTest (shipped config == built-in defaults; consensus fires).
- `report-tool` CI job: `pip install ./python` + render the showcase fixture to a valid PDF.
- Actions bumped to Node 24 runtimes (checkout v5, setup-micromamba v3).

## [0.0.3] / [0.0.2] / [0.0.1]

Initial pre-release line — the C++ adjudicator (2-axis projection, reference-bias rescue +
circularity firewall), the `adjudicate` / `combine` / `benchmark` / `ablate` subcommands,
and the first tracked benchmark suite. See the git history for details.
