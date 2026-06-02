# Function / module I/O

The main modules and their input → output contracts. Types are defined under
[`include/panisoguard/`](../include/panisoguard/).

## Input readers

| Module | Input | Output | Purpose |
|--------|-------|--------|---------|
| SQANTI reader (`sqanti_reader.cpp`) | SQANTI3 `*_classification.txt` | `SqantiTable` (`SqantiRecord[]`) | Load structural category + QC priors (name-indexed columns). |
| GTF/BED reader (`gtf_reader.cpp` / `bed12_reader.cpp`) | caller isoform GTF or BED12 | `isoform_id → IntronChain` | Convert transcript structures into normalized 0-based half-open intron chains. |
| Catalog builder (`gtf_reader.cpp`) | reference GTF | `Catalog` | Build the known-intron / known-transcript catalog (strand-aware `has_intron`). |
| SJ.tab reader (`sj_tab_reader.cpp`) | STAR `SJ.out.tab` | `SjTable` | Short-read junction set with exact-coordinate lookup. |
| Pangenome reader (`pangenome_reader.cpp`) | graph-junction TSV | `PangenomeJunctions` | Junctions present on a pangenome graph haplotype path. |
| FASTA fetcher (`variant_motif.cpp`) | genome + haplotype FASTA (`faidx`) | `FastaFetcher` / `HaplotypeProvider` | Compare splice motifs between reference and a personalized haplotype. |
| BAM reader (`bam_features.cpp`) | indexed BAM/CRAM (htslib) | `BamReader` | Per-junction read-level mapping features. |

## Core

| Module | Input | Output | Purpose |
|--------|-------|--------|---------|
| Fingerprint (`fingerprint.cpp`) | `IntronChain` | 64-bit hash | FNV-1a over chrom+strand+sorted introns; caller-agnostic identity. |
| Consensus / combine (`consensus.cpp`) | multiple callers' chains | `PIG.NNNNNN` ids + caller-support matrix | Integrate isoforms across callers by fingerprint. |
| Adjudicator (`adjudicator.cpp`) | `SqantiTable`, `id→IntronChain`, optional `Catalog`/`SjTable`/`BamReader`/`HaplotypeProvider`/`PangenomeJunctions` | `AdjudicationResult[]` | Per isoform: detect novel junctions, assemble `EvidenceVector`, apply `RuleEngine`. |
| RuleEngine (`rules.cpp`) | `EvidenceVector` | `Verdict` | Deterministic projection → confidence class + primary mechanism + circularity flag + `rule_trace`. |

## Analysis (reviewer-facing)

| Module | Input | Output | Purpose |
|--------|-------|--------|---------|
| Non-redundancy (`analysis.cpp`, `benchmark`) | adjudication results + SQANTI3 `filter_result` | `NonRedundancy` (2×2, Jaccard, McNemar) | Show PanIsoGuard is orthogonal to the SQANTI3 filter, not redundant. |
| Ablation (`analysis.cpp`, `ablate`) | assembled evidence + `RuleEngine` | `AxisAblation` (per-axis class-change counts) | Quantify each axis's marginal contribution. |

## Output writer

| Module | Input | Output | Purpose |
|--------|-------|--------|---------|
| Result writer (`result_writer.cpp`) | `AdjudicationResult[]` + `RunProvenance` | `<prefix>.adjudicated.tsv`, `.attribution.jsonl`, `.provenance.log` | TSV verdicts; per-isoform `rule_trace` + `graph_trace` ([method_graph.md](method_graph.md)) + `bio_flags` (SQANTI3 descriptors passed through, verdict-neutral; [relationship_to_sqanti3.md](relationship_to_sqanti3.md)) in the JSONL; and run provenance. |

## Key data types

| Type | Fields (abridged) |
|------|-------------------|
| `IntronChain` | `chrom`, `strand`, `introns: Junction[]` (0-based half-open) |
| `EvidenceVector` | novel-junction counts, short-read support count, BAM worst-case fractions, SQANTI prior flags, variant/pangenome rescue flags + circularity-risk state |
| `Verdict` | `confidence: ConfidenceClass`, `primary_mechanism: Mechanism`, `circularity_flag: bool`, `rule_trace: string[]` |
| `AdjudicationResult` | `isoform_id`, `chrom`, `strand`, `structural_category`, `evidence`, `verdict` |
