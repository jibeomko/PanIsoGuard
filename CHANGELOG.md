# Changelog

All notable, user-facing changes. Format follows [Keep a Changelog](https://keepachangelog.com);
this project uses [Semantic Versioning](https://semver.org) (pre-1.0: minor/patch only).

## [Unreleased] — targets 0.0.4

The flagship change is the **multi-caller consensus axis**, plus a PDF report tool, a
container, and a large round of honest validation. Release is held until the first
bioconda submission (v0.0.3) merges; see [docs/releasing.md](docs/releasing.md).

### Added
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
  `wholegenome_multicaller` (real GM12878), `merge_comparison` (head-to-head vs
  gffcompare / TAMA), `hg002_wholegenome` (variant-rescue yield), and
  `giab_cohort_rescue` (rescue across a 4-individual GIAB cohort: 137/137, 0 false).

### Fixed
- **`config/rules.default.toml` silently disabled the consensus axis** — `consensus_min_callers`
  was left as a `999` "reserved" sentinel after the axis was wired, so passing
  `--config config/rules.default.toml` turned consensus off. Set to the built-in default (2)
  and guarded by a new `integration_config_equivalence` regression test.
- **Report page-2 projection grid** mis-coloured rescue-mechanism cells
  (`variant_created` / `population_known` bypass the grid); restricted to the artifact
  mechanisms.
- Hardened junction aggregation and provenance.

### Changed
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
