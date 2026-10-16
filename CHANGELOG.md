# Changelog

All notable, user-facing changes. Format follows [Keep a Changelog](https://keepachangelog.com);
this project uses [Semantic Versioning](https://semver.org) (pre-1.0: minor/patch only).

## [Unreleased]

### Added
- **Whole-chain long-read support (ruleset `builtin-0.0.4`).** With `--bam`, the adjudicator
  counts the primary reads whose alignment carries the isoform's entire intron chain
  (`evidence.chain_n_full_reads`, `chain_n_partial_reads`). When short-read support is
  `UNKNOWN` (no SJ.tab, or a NIC isoform whose junctions are all known) and no artifact
  sign is present, ≥ `chain_min_full_reads` (3) lifts the call to `NOVEL_PROBABLE`
  (`chain_support`; `chain_support_with_<mechanism>` → `UNCONFIRMED`). The reads are the
  caller's own, so it never reaches `NOVEL_CONFIRMED` and never overrides a short-read
  verdict. On SQANTI-SIM chr22 with 1 M short-read pairs it lifts FLAIR's recall from
  0.545 to 0.710 at precision 0.992 (120 genuine / 0 false isoforms), IsoQuant 92 / 0,
  ESPRESSO 106 / 1; Bambu (66 / 13) and TALON (122 / 48) pay precision because their false
  chains are read-supported too ([benchmark/chain_support](benchmark/chain_support)).
  `ablate --axes chain`; TOML `[axis_chain] min_full_reads`.
- **`genomic_compatibility`** column (last in `adjudicated.tsv`, also in the JSON): the DNA
  side on its own, independent of the verdict: `haplotype` / `pangenome` (every novel
  junction explained, on one haplotype or path), `partial`, `none`, `NA`. A firewall-held
  isoform reads `UNRESOLVED` + `haplotype`: genomically compatible, RNA unconfirmed.
- **`<prefix>.junctions.tsv`**: one row per junction of the adjudicated isoforms with the
  per-junction evidence the verdict aggregates (isoforms using it, known/novel, short-read
  count, reference and per-haplotype splice motifs with `created` / `disrupted`, pangenome
  support and paths, BAM fractions, `reference_bias`). On HG03516 the rescue criterion
  holds at 35 junctions but only 13 sit in rescued isoforms (the rest share an isoform with
  an unexplained junction); the table keeps them visible.
- `benchmark/pangenome/derive_junctions.py` writes the haplotype/path ids (6th column) from
  the VCF genotypes, so the path-consistency rule runs on the HPRC chr22 set (3,034 of 3,040
  rows carry ids; results unchanged: 0 rescues on the real FLAIR set, gate1 6/6).

### Benchmarks
- The head-to-head ([benchmark/sqanti3_filter_h2h](benchmark/sqanti3_filter_h2h)) re-run
  with the chain axis: at 1 M pairs PanIsoGuard + BAM goes from precision 0.990 / recall
  0.545 to **0.992 / 0.710** (AUPRC 0.972, still tied with the one-line rule); with no
  short reads from 0 calls to 489 genuine / 4 false (recall 0.67, AUPRC 0.947 vs SQANTI3
  rules 0.909, bootstrap CI [+0.023, +0.055]). Per-class at 1 M + BAM: NOVEL_PROBABLE
  167 / 1, UNRESOLVED 0 / 6 (the genuine held NIC moved up), UNCONFIRMED 203 / 46. The
  envelope now carries a risk–coverage point per method (PanIsoGuard + BAM decides 71 %
  of isoforms with 2.1 % wrong; the one-line rule decides all with 6.6 % wrong).
- [benchmark/refbias_3way](benchmark/refbias_3way) arm C re-run with the chain axis: junction-level
  numbers and the 18 / 28 rescues unchanged; 5,784 / 7,011 isoforms that SQANTI3's filter
  keeps and PanIsoGuard used to hold `UNRESOLVED` are now `NOVEL_PROBABLE` (whole-chain
  reads, no short reads in these samples); `adjudicate` takes 3–6 min per genome-wide run
  instead of 2.
- [benchmark/refbias_3way](benchmark/refbias_3way) replicated on HG02717: rule level 37/39
  SNV-class junctions (precision 0.88), verdict level 21/39 (precision 0.95), the same
  miss pattern as HG03516; the envelope now holds both individuals. B2's SQANTI3 step had
  failed silently in both runs (frame-less liftoff CDS lines); fixed in `run.sh`, no metric
  depended on it, the cost row did.

### Fixed
- `panisoguard example` no longer uses `std::filesystem`, which bioconda's osx-64 build
  (macOS 10.9 deployment target, libc++ filesystem needs 10.15) cannot compile; the 0.0.6
  autobump PR failed on osx-64 for this reason. POSIX calls do the same job.

### Changed
- **Variant-disrupted motifs are consumed.** A novel junction canonical on the reference
  but non-canonical on every haplotype of the sample counts as the `noncanonical` sign
  (reason suffix `;haplotype_disrupted`), under the same provenance firewall as the rescue
  and off with `ablate --axes variant`. Before, `HaplotypeProvider` computed it and the
  adjudicator dropped it. `evidence` gains `n_novel_jx_variant_disrupted`. No verdict on the
  tracked sets changes (support is unknown without short reads there).

## [0.0.6] - 2026-10-14

### Changed (breaking: output names and columns)
- **Confidence classes renamed.** Names now say the verdict; the evidence is in the new
  columns. `HIGH_CONF_KNOWN` → `KNOWN`, `HIGH_CONF_NOVEL` → `NOVEL_CONFIRMED`,
  `MEDIUM_CONF_NOVEL` → `NOVEL_PROBABLE`, `LOW_CONF_PARTIAL` → `UNCONFIRMED`,
  `AMBIGUOUS` → `UNRESOLVED`, `PAN_REF_RESCUED_FALSE_NOVEL` → `REFERENCE_BIAS` (the
  rescued junction is real; what was false was the artifact flag, not the novelty: 77 of
  the 86 rescued HG03516/HG02717 junctions are novel on the personal genome too).
  `ARTIFACT` is unchanged. Class membership is unchanged except where the rules below say.
- **Two new columns** in `adjudicated.tsv` / `attribution.jsonl`: `verdict` (known /
  confirmed / reference_bias / unconfirmed / artifact / not_adjudicated) and
  `verdict_reason` (the one rule that decided it, e.g. `partial_support`, `no_coverage`,
  `firewall_held`, `ism_intrapriming`). They split the situations the classes mix:
  `UNCONFIRMED` is ISM (0.92 precision on SQANTI-SIM) plus weak novels (0.36–0.56), and
  `UNRESOLVED` is "nothing to test" (0.96) plus firewall holds plus unadjudicated
  categories. `provenance.log` adds verdict counts.
- **Rules (ruleset `0.0.3`).**
  - The BAM mapping gate needs ≥ `min_spanning` (3) reads at a novel junction; with 1–2
    reads a single indel made the fraction 1.0 and decided the verdict alone (45 % of
    HG03516's mapping-artifact calls rested on < 3 reads).
  - `UNSUPPORTED` needs coverage: if none of the isoform's known junctions has short-read
    support either, support is `UNKNOWN` (`no_coverage`), not evidence against.
  - Intra-priming (`degradation`) needs perc_A ≥ 60 **and** no polyA motif upstream
    (SQANTI3 `polyA_motif_found`, when it ran with `--polyA_motif_list`); 20–30 % of
    perc_A-flagged 3' ends carry a motif and are no longer flagged. The column moved from
    `bio_flags` to `evidence`.
  - A reference-bias rescue carries a BAM mapping sign in its reason
    (`haplotype_variant;mapping_sign`) instead of hiding it; `evidence` gains
    `n_novel_jx_variant_created`, `n_known_junctions`, `n_known_jx_sr_supported`,
    `bam_n_junctions_gated`.
  - `--platform ont` preset: turns the HiFi-calibrated indel-near gate off (it fires on
    ~90 % of ONT junctions).
  - **A reference-bias rescue needs one haplotype (or one pangenome path) carrying every
    novel junction**, not a different one per junction: the per-junction haplotype sets
    are intersected (`evidence.haplotype_consistent`); every-junction-created but no
    common haplotype is held `UNRESOLVED` (`haplotype_inconsistent`). The pangenome
    junction file takes an optional 6th column of path ids; with it, the same rule
    applies (`pangenome_path_consistent`, `path_inconsistent`), without it the counts
    decide as before.
  - **Caller consensus promotes only a chain that was examined**: with the caller chain
    absent or no `--ref-gtf`, `n_callers` is reported but the call stays `UNRESOLVED`.
  - **`--caller-support` is joined by intron-chain fingerprint** (new `fingerprint`
    column in the `combine` matrix), not by native id: two callers reusing one name
    (`transcript_1`) for different chains no longer hand one chain's support to the
    other. Older matrices still join by native id; a name reused across chains is
    warned about.
  - SQANTI3 `RTS_stage` / `all_canonical` / `polyA_motif_found` = `NA` are
    not_evaluable (JSON `null`), no longer `false`.
- `combine` orders monoexonic isoforms by (caller, native id), so their `PIG.NNNNNN`
  ids no longer depend on input order.
- `ablate --axes` rejects an unknown axis name with an error instead of silently ablating
  nothing.
- `provenance.log` records every input (also the reference, haplotype FASTAs and the
  pangenome file) with its size and mtime, and an MD5 per file with `--hash-inputs`.

### Benchmarks (all tracked envelopes re-run with this binary)
- **New: [benchmark/refbias_3way](benchmark/refbias_3way).** The reference-bias rescue
  against an independent standard: HG03516's reads through GRCh38, an SNV-personalized
  diploid, the HPRC phased assembly, and PanIsoGuard post hoc; truth per junction = where
  the same reads splice on the personal genome. Rule level 31/32 SNV-class junctions
  (precision 31/35), verdict level 12/32 (precision 12/13); indel/SV-shifted junctions are
  out of scope by construction; 5 min / 0.4 GB vs 11–12 h / 38–78 GB.
- The `min_spanning` gate changes two earlier numbers: the BAM indel-near trigger now
  catches 14/124 false novel junctions on SQANTI-SIM (was 41/124; still 0 genuine flagged),
  and PanIsoGuard+BAM's AUPRC gain over the one-line short-read rule is within the
  bootstrap interval at every depth (was significant at 250k pairs). Operating points are
  unchanged (precision 0.990, recall 0.545). Per-class precision at 1 M pairs with the
  BAM: NOVEL_CONFIRMED 0.992, NOVEL_PROBABLE 0.979, UNRESOLVED 0.962, UNCONFIRMED 0.789,
  ARTIFACT 0.097.
- The haplotype-consistency rule changed no verdict on the real reference-bias sets
  (HG03516/HG02717 cohort 86/86, GIAB whole-genome cohort 137/137, HG002 chr22, arm C of
  the 3-way benchmark: all unchanged).
- `benchmark/sqanti3_filter_h2h/emit_metrics.py` and `score.py --emit-metrics` of the 3-way
  benchmark now write their envelopes (the head-to-head bootstrap was a scratch script).
- README no longer calls the classes an ordered scale (HIGH > MEDIUM > LOW); on
  SQANTI-SIM at every short-read depth `UNRESOLVED` is more precise than `UNCONFIRMED`
  (`benchmark/sqanti3_filter_h2h/per_class.py`).

### Added
- **Bundled examples.** `cmake --install` now ships `examples/tiny`, `examples/multi_caller`
  and the study notes' toy data under `<prefix>/share/panisoguard`, so the bioconda package
  and the Docker image carry them. `panisoguard example` lists them and
  `panisoguard example tiny` copies one into the current folder, ready to run
  (`tiny/run.sh`). The bioconda recipe test and the CI container job run an example this way.

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
