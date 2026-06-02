# Relationship to SQANTI3 (what PanIsoGuard does *not* recompute)

> **PanIsoGuard does not recompute SQANTI3 QC descriptors. It consumes them as priors
> and integrates them with independent path-level evidence from caller consensus,
> short-read junction support, long-read mapping, personalized haplotypes, and
> pangenome-supported junctions.**

[SQANTI3](https://github.com/ConesaLab/SQANTI3) is a mature transcript-QC tool: its
classification file already provides transcript-level QC features and its junctions
file provides junction-level attributes. PanIsoGuard is a **decision layer on top of**
that output, not a reimplementation of it. To avoid being read as "SQANTI3 rebuilt", it
deliberately **does not** recompute any of SQANTI3's descriptors.

## Overlap map

| Feature | Already in SQANTI3 | PanIsoGuard's stance |
|---------|--------------------|----------------------|
| TSS/TTS support (`dist_to_CAGE_peak`, `within_CAGE_peak`, `dist_to_polyA_site`, `within_polyA_site`, `ratio_TSS`) | yes | **consume, don't recompute** — surfaced verbatim as `bio_flags` |
| polyA motif / intra-priming (`polyA_motif`, `polyA_dist`, `polyA_motif_found`, `perc_A_downstreamTTS`) | yes | **consume** — `perc_A_downstream_TTS` is a degradation prior; the rest are `bio_flags` |
| ORF / coding / NMD (`coding`, `ORF_length`, `CDS_*`, `predicted_NMD`) | yes (TransDecoder2) | **consume** — `predicted_NMD` surfaced as a `bio_flags` annotation |
| splice-site motif / junction QC (`splice_site`, `RTS_junction`, `indel_near_junct`, `junction_category`) | yes (junctions file) | **consume** — `all_canonical` / `RTS_stage` are mechanism priors |
| SQANTI-prior artifact rules | yes | already consumed: `RTS_stage`, `all_canonical`, `perc_A_downstream_TTS` |
| short-read splice-junction support | partial | PanIsoGuard adds **per-edge** `SJ.tab` corroboration of *novel* junctions |
| **multi-caller intron-chain consensus** | no | **PanIsoGuard differentiator** (`combine`) |
| **graph / path-level caller recurrence** | no | **PanIsoGuard differentiator** |
| **pangenome / haplotype reference-bias rescue + circularity firewall** | no | **PanIsoGuard's most independent contribution** |

Sources: SQANTI3 QC output reference
(<https://github-wiki-see.page/m/ConesaLab/SQANTI3/wiki/Understanding-the-output-of-SQANTI3-QC>)
and the SQANTI3 *Nature Methods* paper
(<https://www.nature.com/articles/s41592-024-02229-2>), which describe the QC module's
structural category, TSS/TTS annotation, non-canonical junctions, intra-priming,
RT-switching, CAGE/Quant-seq, and short-read support.

## Where the line is drawn

PanIsoGuard's contribution is the **adjudication of an evidence-annotated path through
a gene-local splice graph** (see [method_graph.md](method_graph.md)) — integrating
orthogonal, largely SQANTI-independent evidence:

- **caller consensus** — the same intron-chain path supported by multiple callers
  (`combine`); SQANTI3 is single-caller.
- **short-read edge support** — per-novel-junction `SJ.tab` corroboration.
- **long-read mapping** — read-level MAPQ / supplementary / soft-clip / indel features
  on the reads that traverse a junction.
- **personalized haplotype** — a novel junction that is non-canonical on the reference
  but canonical on the sample's haplotype (reference bias, not new splicing).
- **pangenome support** — a junction realizable on a population-graph path though absent
  from the linear reference.
- **circularity firewall** — rescues that would rest on the sample's own RNA-derived
  data are held `AMBIGUOUS`, never promoted.

It explicitly avoids re-deriving TSS/TTS, ORF/NMD, polyA motif, or splice motif, and
makes no SQANTI3-style ML/rule QC filter — those would duplicate SQANTI3.

## `bio_flags` (consume, not recompute)

To make the "consume, don't recompute" stance concrete, each
`<prefix>.attribution.jsonl` row carries a `bio_flags` object that **passes through**
the relevant SQANTI3 descriptors verbatim (`null` when the column is absent or `NA`):

```json
"bio_flags": {
  "dist_to_CAGE_peak": 5,
  "dist_to_polyA_site": 8,
  "predicted_NMD": null,
  "polyA_motif_found": null
}
```

These are **biological-plausibility annotations only** — read from the SQANTI3
classification, never computed by PanIsoGuard, and they do **not** enter the rule engine
or change the verdict (which remains the deterministic path-level projection in
[decision_engine.md](decision_engine.md)). They let a downstream user weigh SQANTI3's QC
context alongside PanIsoGuard's independent verdict without PanIsoGuard re-deriving it.
