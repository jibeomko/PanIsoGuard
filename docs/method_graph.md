# A splice-graph view of PanIsoGuard

> **PanIsoGuard represents each candidate isoform as a path through a gene-local
> splice graph and adjudicates that path using orthogonal edge- and path-level
> evidence.**

This document is a *framing*, not a new algorithm. The deterministic rule engine
described in [decision_engine.md](decision_engine.md) and [algorithm.md](algorithm.md)
is unchanged; this page re-states what it already does in graph-theoretic terms,
because the splice-graph view makes the method's structure explicit and gives each
verdict a compact, machine-readable `graph_trace` (see [Output](#output)).

It deliberately stops at **splice graph / path / edge support / graph distance**. It
introduces no learned graph model, message passing, centrality, or PageRank: those add
no biological signal here and would only invite "why is this needed?". The contribution
remains the *deterministic adjudication of an evidence-annotated path*.

## The gene-local splice graph

For a locus, build a directed graph:

- **nodes** — transcript boundaries and splice sites (donor / acceptor positions, plus
  the path source and sink).
- **edges** — **intron edges** (a splice junction connecting a donor to an acceptor)
  and exon segments between them.
- **a transcript is a path** — an ordered source-to-sink walk; equivalently, its
  *ordered intron-edge sequence* (its intron chain). PanIsoGuard already normalizes
  every isoform to this intron chain
  ([`types.hpp`](../include/panisoguard/types.hpp): `Junction`, `IntronChain`), so the
  "path" is exactly the object the engine operates on.

The **reference annotation is the known splice graph**: the union of all annotated
intron edges and splice-site nodes
([`Catalog`](../include/panisoguard/consensus.hpp), built from the reference GTF). An
edge is *known* if `Catalog::has_intron` accepts it (exact, strand-aware), *novel*
otherwise.

| PanIsoGuard concept | Splice-graph object |
|---------------------|---------------------|
| caller isoform | candidate path (ordered intron-edge sequence) |
| reference annotation | known splice graph (known edges + known site nodes) |
| novel junction | **novel edge** (intron not in the known graph) |
| caller consensus (`combine`) | one path supported by several callers (shared fingerprint) |
| short-read `SJ.tab` | **edge support** (an independent observation of an intron edge) |
| BAM spanning reads | **path/edge read-quality** at the junction |
| SQANTI3 QC priors | **path-level priors** (RTS, canonicality, intra-priming) |
| variant / haplotype | an edge that is **non-canonical on the reference graph but canonical on the personal graph** |
| pangenome junction | an edge present on a **population graph** path though absent from the linear reference |
| verdict | a **deterministic projection** of the evidence-annotated path |

## Novelty as graph distance

A candidate path's novelty is its **edge-set distance to the known splice graph** — the
number of its intron edges that are not known edges:

```text
path_distance_to_reference = | edges(path) \ edges(known graph) | = novel_edges
```

The *kind* of novelty is a property of the novel edges' endpoint nodes:

- **known-site recombination** — every novel edge connects two *known* splice-site
  nodes in a new pairing (a novel intron built from known sites). This is the
  `novel_in_catalog` (NIC) case: combinatorial novelty with no new node.
- **novel splice site** — at least one novel edge introduces a *new* donor/acceptor
  node. This is the `novel_not_in_catalog` (NNC) case.

This distinction matters because the short-read axis can corroborate an *edge* but the
combinatorial *path* of a NIC isoform may carry no individually-novel edge to
corroborate — which is exactly why the engine honestly **abstains** (`AMBIGUOUS`) on
NIC-combinatorial novelty rather than guessing
(see [decision_engine.md](decision_engine.md), "NIC behaviour").

## Evidence-annotated edges and paths

Each axis annotates an edge or the whole path with an orthogonal, independently-skippable
signal (an absent input is *not_evaluable*, never silently counted for or against):

- **edge support** — short-read `SJ.tab` exact match (`SjTable::find_exact`): an
  independent observation that the intron edge exists.
- **path/edge read-quality** — BAM spanning-read features (`BamReader`): low-MAPQ /
  supplementary / soft-clip / indel-near fractions on the reads that traverse the edge.
- **path-level priors** — SQANTI3 `all_canonical` / `RTS_stage` /
  `perc_A_downstream_TTS` annotate the path as a whole.
- **personal-graph canonicalization** — `HaplotypeProvider::classify`: a novel edge
  whose motif is non-canonical on the reference but canonical on the sample's
  haplotype, i.e. the edge is real on the *personal* splice graph.
- **population-graph support** — `PangenomeJunctions::supports`: the novel edge lies on
  a haplotype path of a population graph though it is absent from the linear reference.

## Artifacts as graph pathologies

The `ARTIFACT` / held outcomes have clean graph readings:

- **unsupported bridge edge** — a novel edge with no short-read support: a bridge added
  to the graph that nothing independent observes.
- **non-canonical edge without graph rescue** — a non-canonical novel edge that is
  *not* canonicalized on the personal graph nor present on the population graph.
- **low-confidence path** — the reads traversing the edge map poorly (mapping
  mechanism); the path may be an alignment artifact, so upstream concerns dominate.
- **circular evidence edge** — the only support for a rescue edge derives from the
  *same sample's* RNA/graph; the circularity firewall holds it `AMBIGUOUS` rather than
  promoting it (variant and pangenome axes both gated — see
  [decision_engine.md](decision_engine.md)).

The **reference-bias rescue is all-or-nothing on the path**: a path is rescued only if
*every* novel edge is explained (canonical on the personal graph, or present on the
population graph). One explained edge does not rescue a path whose other edges are
genuinely new — this is the early-return rescue in
[`rules.cpp`](../src/core/rules.cpp).

## Output

The verdict is unchanged. Each isoform's `<prefix>.attribution.jsonl` row additionally
carries a `graph_trace` object — a re-expression of the already-computed
`EvidenceVector` under the graph framing (no new computation):

```json
"graph_trace": {
  "novel_edges": 1,
  "path_distance_to_reference": 1,
  "edges_sr_supported": 1,
  "novel_site_type": "novel_splice_site",
  "pangenome_edge_support": false,
  "variant_canonicalized_edge": false
}
```

| field | meaning | source field |
|-------|---------|--------------|
| `novel_edges` | intron edges not in the known graph | `EvidenceVector::n_novel_junctions` |
| `path_distance_to_reference` | edge-set distance to the known graph (= `novel_edges`) | `n_novel_junctions` |
| `edges_sr_supported` | novel edges with short-read support | `n_novel_jx_sr_supported` |
| `novel_site_type` | `known_site_recombination` (NIC) / `novel_splice_site` (NNC) / `none` | `structural_category` |
| `pangenome_edge_support` | all novel edges on a population-graph path | `pangenome_rescue` |
| `variant_canonicalized_edge` | all novel edges canonical on the personal graph | `variant_rescue` |

## See also

- [decision_engine.md](decision_engine.md) — the 2-axis grid, rescue precedence, and the
  circularity firewall (the decision this framing describes).
- [algorithm.md](algorithm.md) — the per-isoform loop that builds and adjudicates the path.
- [architecture.md](architecture.md) — where the splice graph (catalog) and path
  (intron chain) live in the code.
