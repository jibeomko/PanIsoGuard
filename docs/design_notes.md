# Design notes

What PanIsoGuard adds to SQANTI3's output, what it deliberately does not recompute, and the splice-graph way of reading its rules.

PanIsoGuard takes SQANTI3's QC values as they are and combines them with other, independent
evidence: caller agreement, short-read junctions, long-read alignments, the sample's
haplotypes, and pangenome junctions. It does not recompute TSS/TTS, ORF/NMD, polyA or reference
splice motifs, and it is not a SQANTI3-style QC filter (see
[docs/relationship_to_sqanti3.md](relationship_to_sqanti3.md)). What it adds is the
**decision logic**: fixed, readable rules that turn the evidence into a verdict. The thresholds
are in [`config/rules.default.toml`](../config/rules.default.toml) and can be changed without
rebuilding, and every verdict records its `rule_trace` (see
[docs/decision_engine.md](decision_engine.md)).

> PanIsoGuard is an **independent project**; the SQANTI3 authors have not made or endorsed it.
> It only reads SQANTI3's output files. No SQANTI3 code is included or linked, so SQANTI3's
> GPL-3.0 license does not extend to PanIsoGuard's MIT code. If you use SQANTI3 in your
> pipeline, please cite it; see
> [docs/relationship_to_sqanti3.md](relationship_to_sqanti3.md#attribution-licensing--citation).

You can also picture each isoform as a **path through its gene's splice graph**: a novel
junction is an edge that the reference graph does not have. This is only another way to
describe the same rules; there is no graph model or new algorithm.
