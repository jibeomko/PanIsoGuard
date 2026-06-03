# The PanIsoGuard decision engine

This document describes the adjudication logic — the methodological contribution.
It is intentionally readable so any verdict can be audited without reading C++:
all thresholds live in [`config/rules.default.toml`](../config/rules.default.toml)
and the engine emits a per-isoform `rule_trace`.

## Framing constraints (non-negotiable)

- PanIsoGuard reports **orthogonal attribution + confidence + provenance**.
- It **never** claims to beat the SQANTI3 random-forest filter or FLAIR2. The
  `benchmark` module reports *non-redundancy* (McNemar / 2×2 / Jaccard), not
  superiority.
- It **does not recompute SQANTI3 QC descriptors** (TSS/TTS, ORF/NMD, polyA/splice
  motif): those are consumed as priors and surfaced verbatim as `bio_flags`, while the
  verdict rests on independent path-level evidence
  ([relationship_to_sqanti3.md](relationship_to_sqanti3.md)).
- Validation uses **public + simulated** data (SQANTI-SIM; HG002/HPRC; LRGASP). A
  clinical cohort is at most a real-world demo / engineering scale test, never the
  source of mechanistic ground truth.

## Two-axis evidence space

Each *novel* isoform (SQANTI `novel_in_catalog` / `novel_not_in_catalog`) is scored
on two axes; the binned pair is projected to one confidence class.

- **Axis A — novelty support** (from short-read `SJ.tab` corroboration of the
  isoform's novel junctions): `SUPPORTED` / `PARTIAL` / `UNSUPPORTED`, or `UNKNOWN`
  when the axis is not evaluable. When short-read corroboration is `UNKNOWN`,
  **multi-caller consensus** (`--caller-support`, from a `panisoguard combine`
  matrix) can stand in as a *weaker, methodological* form of support: a novel chain
  independently recovered by ≥ `consensus_min_callers` callers is corroborated, but
  only enough to reach `MEDIUM_CONF_NOVEL`, never `HIGH` (caller agreement is not
  orthogonal experimental evidence).
- **Axis B — artifact mechanism**: `none` / `mapping_or_repeat` / `noncanonical` /
  `rt_switch` / `degradation` / `variant_created`. An axis whose input is absent is
  `not_evaluable` — a soft-skip, never silently counted for or against.

## Evidence axes and where they come from

| Mechanism | Source | Evaluable when |
|-----------|--------|----------------|
| short-read support | STAR `SJ.tab` exact junction match (strand-aware) vs the reference catalog | `--sj-tab` + `--ref-gtf` |
| mapping_or_repeat | BAM reads spanning the junction — the **low-MAPQ / supplementary** fractions gate this mechanism; soft-clip / indel-near fractions are also computed and reported but do not currently drive the verdict | `--bam` |
| noncanonical / rt_switch / degradation | SQANTI3 priors `all_canonical` / `RTS_stage` / `perc_A_downstream_TTS≥60` | always (from classification) |
| variant_created | reference vs personalized haplotype FASTA splice-motif comparison | `--reference` + `--reference-haplotype` |
| caller consensus | a `panisoguard combine` matrix: how many independent callers recovered the isoform's intron chain (`n_callers`) | `--caller-support` |

## Mechanism priority and its rationale

When several mechanisms fire, one is reported as `primary_mechanism`; the order is
**mapping > noncanonical > rt_switch > degradation** (and the variant rescue, below,
preempts all of them). Rationale: **alignment reliability is the most upstream
concern** — if the spanning reads do not map confidently, the junction may be a
mapping artifact and downstream motif/QC signals are moot; non-canonical motif is the
next strongest structural signal; RT-switching and intra-priming degradation are
softer biochemical signals. Only the first-matching mechanism is the `primary`; all
fired conditions are listed in `rule_trace`.

## Projection to confidence classes (7 reachable)

1. **Reference-bias rescue (highest priority).** A novel-vs-linear-reference
   junction explained by reference bias is `PAN_REF_RESCUED_FALSE_NOVEL`, via either:
   - *Pangenome axis* (`--pangenome-junctions`): **all** of the isoform's novel
     junctions are realizable on a pangenome graph haplotype path (mechanism
     `population_known`) — consistent with reference bias, not proof the splice is
     used. Independent only if the junction set is from population assemblies, which
     PanIsoGuard cannot verify, so it is gated by the **same circularity firewall**
     as the variant axis: `--pangenome-provenance population`/`external` promote,
     while `sample_derived`/`unknown` (the default) are held `AMBIGUOUS`. It is
     evaluated first, so an *independent* graph rescues even when the sample's own
     haplotype provenance is circular-risk. **Validated** on the real HPRC v1.1 chr22
     graph (GATE-1; see [benchmark/pangenome](../benchmark/pangenome)).
   - *Variant axis* (`--reference-haplotype`): the junction is non-canonical on the
     reference but canonical on a personalized haplotype (mechanism
     `variant_created`). **Circularity firewall:** if the haplotype provenance is
     RNA-derived/unknown this is *not* independent evidence, so it is held
     `AMBIGUOUS` with `circularity_flag=true` and never promoted; only
     `wgs`/`external` provenance promotes.
2. Otherwise project (Axis A × Axis B):
   - `SUPPORTED × none → HIGH_CONF_NOVEL`
   - `SUPPORTED × mechanism → MEDIUM_CONF_NOVEL`
   - `PARTIAL × none → MEDIUM_CONF_NOVEL`
   - `PARTIAL × mechanism → LOW_CONF_PARTIAL`
   - `UNSUPPORTED × none → LOW_CONF_PARTIAL`
   - `UNSUPPORTED × mechanism → ARTIFACT`
   - `UNKNOWN` support → `AMBIGUOUS`, except:
     - `UNKNOWN × strong mapping artifact → ARTIFACT` (a decisive mapping signal
       stands without short-read support); and
     - `UNKNOWN × caller-consensus (≥ consensus_min_callers, no mapping artifact) →
       MEDIUM_CONF_NOVEL` (mechanism `none`) or `LOW_CONF_PARTIAL` (any other
       mechanism). Multi-caller agreement corroborates the chain methodologically
       when the short-read axis is silent; a strong mapping artifact still preempts
       it. Disable with `ablate --without consensus`.
- Pass-through: `full-splice_match → HIGH_CONF_KNOWN`; `incomplete-splice_match →
  LOW_CONF_PARTIAL`; non-targeted categories → `AMBIGUOUS`.

**`ARTIFACT` semantics (read carefully).** `ARTIFACT` denotes *(UNSUPPORTED × any
mechanism)* or *(UNKNOWN × strong mapping)* — i.e. a **combined confidence ×
mechanism** signal, **not** a standalone claim that the call is spurious. A novel
isoform with an unusual motif but no short-read backing lands here; users should
cross-check `ARTIFACT` calls against the original reads before hard filtering.

**NIC behaviour.** A `novel_in_catalog` isoform's novel junction is a *novel pairing
of known splice sites*, i.e. a novel intron not present in the catalog, so it IS
counted as a novel junction and DOES receive short-read / BAM / variant evaluation
(it is not auto-`AMBIGUOUS`). The residual `AMBIGUOUS` set is non-targeted categories
plus isoforms whose every junction is an already-known intron (no novel junction to
evaluate). `rule_trace` records which `UNKNOWN` reason applies
(`caller_chain_absent` / `axis_absent` / `no_novel_junctions`).

## Provenance and circularity

Every reclassification records its driving axis and (for the variant axis) a
provenance/circularity flag. A verdict resting solely on variants derived from the
same RNA reads being adjudicated is flagged `circularity_risk` and is **not**
promoted to a hard class. `provenance.log` reports which axes were evaluable.

## Evidence gates

- **GATE-0**: the variant-injection proof-of-concept (does a variant create/destroy
  a canonical GT-AG?) is verified deterministically at the motif level
  (`tests/unit/variant_motif_test.cpp`; see [benchmark/variant_inject](../benchmark/variant_inject)).
- **GATE-1 (met)**: the pangenome reference is the **HPRC v1.1 Minigraph-Cactus
  GRCh38** graph, with HG002/NA24385 confirmed **out-of-graph** (absent from the 47
  assembly samples → non-circular). On chr22, `vg deconstruct` → population deletions →
  graph-supported junctions: the rescue produces **0 false rescues** on real FLAIR novel
  junctions, **fires correctly** when a novel intron is a real population deletion
  (population provenance), and is **held by the circularity firewall** under
  circular-risk provenance — see [benchmark/pangenome](../benchmark/pangenome) and
  [benchmark/results/pangenome](../benchmark/results/pangenome).

## Future work

- **In-process pangenome tier** (gbwtgraph/GBZ, `-DWITH_PANGENOME_LIB`): the
  file-based pangenome axis (`--pangenome-junctions`, graph-supported junctions
  pre-extracted with vg/rpvg) is implemented and emits
  `PAN_REF_RESCUED_FALSE_NOVEL` (mechanism `population_known`). Loading and
  traversing a GBZ graph in-process (instead of consuming a pre-extracted junction
  file) remains behind the experimental build flag.
- **Table-driven projection**: move the (Axis A × Axis B) grid fully into the TOML.
- **SQANTI RF concordance**: optionally annotate `rule_trace` with agreement vs the
  SQANTI3 `filter_result`, reinforcing the orthogonal-attribution framing.
- **Variant-disrupted handling**: `HaplotypeProvider` already computes the disrupted
  case; its adjudication policy is not yet defined.
