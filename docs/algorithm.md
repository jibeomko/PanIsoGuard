# Algorithm

The `adjudicate` subcommand processes one SQANTI3 isoform at a time. The core loop is
in [`src/core/adjudicator.cpp`](../src/core/adjudicator.cpp) (`build_evidence` +
`adjudicate`); the decision step is in [`src/core/rules.cpp`](../src/core/rules.cpp)
(`RuleEngine::evaluate`).

> **Graph view.** Equivalently, each isoform is a *path* through a gene-local splice
> graph (its intron chain), the reference catalog is the *known* graph, and a novel
> junction is an *edge* absent from it; the engine adjudicates the evidence-annotated
> path. This framing — and the per-isoform `graph_trace` it adds to the output — is in
> [method_graph.md](method_graph.md). The algorithm below is unchanged by it.

## Per-isoform adjudication

```text
for each SQANTI3 isoform record r:
    1. Join r (structural category + QC priors) with the caller isoform structure,
       matched by isoform id (id -> IntronChain).
    2. Convert the isoform structure into a normalized 0-based half-open intron chain.
    3. Compare each intron against the reference catalog (Catalog::has_intron,
       strand-aware).
    4. Identify the novel junctions = introns not present in the catalog.
    5. For each novel junction, collect evidence from the axes that are available:
        - short-read SJ.tab support     (SjTable::find_exact + read/motif thresholds)
        - long-read BAM mapping features (BamReader::features_at_junction)
        - SQANTI3 priors                (RTS_stage / all_canonical / perc_A_downstream_TTS)
        - haplotype splice-motif change (HaplotypeProvider::classify)
        - pangenome junction support    (PangenomeJunctions::supports)
    6. Store all signals in one EvidenceVector
       (counts of novel vs supported junctions, worst-case BAM fractions, rescue
        flags + their circularity-risk state).
    7. Apply the RuleEngine (deterministic projection; see below).
    8. Emit a Verdict: confidence class, primary mechanism, circularity flag, and a
       rule_trace recording every rule that fired.
```

Notes:

- Novel-junction axes require the caller chain **and** the reference catalog; without
  them, novelty-support is `UNKNOWN` (held, never guessed).
- A junction is matched by **exact, strand-aware coordinates** — no fuzzy/overlap
  matching — so all readers must agree on the 0-based half-open convention (they do,
  by normalizing at read time).
- BAM features are memoized per `(chrom, strand, intron)` so junctions shared across
  isoforms are probed once.

## Decision projection (RuleEngine.evaluate)

Applied in priority order; the full grid and rationale are in
[decision_engine.md](decision_engine.md):

```text
1. Pass-through categories:
     full-splice_match        -> HIGH_CONF_KNOWN
     incomplete-splice_match  -> LOW_CONF_PARTIAL
     non-targeted (not novel) -> AMBIGUOUS
2. Reference-bias rescue (highest priority for novel isoforms):
     if ALL novel junctions are graph-supported (pangenome axis) -> PAN_REF_RESCUED_FALSE_NOVEL
         (held AMBIGUOUS if the junction-set provenance is circular-risk)
     else if ALL novel junctions are variant-explained (haplotype axis) -> PAN_REF_RESCUED_FALSE_NOVEL
         (held AMBIGUOUS if the haplotype provenance is circular-risk)
3. Otherwise project the 2-axis grid:
     Axis A novelty-support  = short-read corroboration of novel junctions
     Axis B artifact mechanism = mapping > noncanonical > rt_switch > degradation
     (support x mechanism) -> {HIGH/MEDIUM/LOW_CONF_NOVEL, AMBIGUOUS, ARTIFACT}
```

The rescue branches return early and are **all-or-nothing**: a single explained
junction does not rescue a multi-novel-junction isoform whose other junctions are
genuinely new.

## Caller integration (combine)

`combine` is a separate pre-step ([`src/core/consensus.cpp`](../src/core/consensus.cpp)):
each caller's isoforms are reduced to intron-chain fingerprints and merged across
callers (which use different running-ID schemes) into a stable `PIG.NNNNNN` id plus a
caller-support matrix. Monoexonic transcripts are kept unmerged (honest, no
coordinate guessing).
