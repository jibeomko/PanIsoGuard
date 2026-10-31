# Input formats

All internal coordinates are **0-based, half-open** introns `[start, end)`. The
readers convert each source so the same intron yields the same coordinates, which
is what makes fingerprint matching and short-read corroboration exact:

- STAR `SJ.tab` (1-based inclusive `s..e`) → `{s-1, e}`
- GTF exons (1-based); intron between exon_i and exon_{i+1} → `{exon_i.end, exon_{i+1}.start-1}`
- BED12 blocks (0-based); intron between blocks → `{block_i.end, block_{i+1}.start}`

## SQANTI3 classification (`--classification`, required)

Tab-separated, with a header. The only mandatory column is `isoform`; everything
else is read by **name** (tolerant of column reordering and appended columns).
Read (as priors, never recomputed): `chrom`, `strand`, `structural_category`, `RTS_stage`,
`all_canonical`, `perc_A_downstream_TTS` and, when present, `polyA_motif_found`. Every
other column is ignored. Category strings are the
SQANTI3 spelling (`full-splice_match`, `incomplete-splice_match`,
`novel_in_catalog`, `novel_not_in_catalog`, …).

## Caller isoforms (`--isoforms-bed` or `--isoforms-gtf`, required)

The structure of the called isoforms, as BED12 or GTF. **Transcript/isoform IDs
must match the SQANTI3 `isoform` IDs** — this is the join key between QC and
structure (e.g. FLAIR's combined IDs `{n}-1_{ENST}_{ENSG}` appear identically in
both files). GTF: `exon` lines grouped by `transcript_id`. BED12: `blockSizes` /
`blockStarts` relative to `chromStart`.

## Reference GTF (`--ref-gtf`, recommended)

A reference annotation (e.g. GENCODE) used to build the known-intron/known-site
catalog. Without it, novel-junction-based axes (short-read, BAM, variant) are not
evaluable and novel isoforms fall to `UNRESOLVED`.

## STAR SJ.tab (`--sj-tab`, optional)

9-column STAR `SJ.out.tab` / merged `combined_SJ.tab`: `chrom, intron_start,
intron_end, strand{0,1,2}, motif{0..6}, annotated, n_uniq, n_multi, max_overhang`.
A novel junction is corroborated if it matches (chrom, strand, coords) exactly with
`n_uniq ≥ sj_min_uniq_reads` and (by default) a canonical motif. A file with no junctions
at all leaves the short-read axis not evaluable (with a warning) rather than marking every
novel junction unsupported.

## BAM/CRAM (`--bam`, optional)

A coordinate-sorted, **indexed** (`.bai`/`.csi`) genome alignment of the long reads
(splice-aware: `N` CIGAR ops at introns). CRAM additionally needs `--reference`.
Read-level features are computed only over reads whose `N` op exactly matches a
novel junction. A read whose aligner recorded the transcript strand (`XS:A` from STAR or
HISAT2, `ts:A` from minimap2, which is relative to the read and flipped for reverse
alignments) counts only for junctions and chains on that strand; a read without the tag
counts for either. A read error in the file (for example a truncated BAM) stops the run
instead of counting as a junction with no reads. Note: a per-sample BAM gives
sample-specific coverage; for a merged isoform set, use a merged BAM for full coverage.

## Personalized haplotype FASTA (`--reference-haplotype`, optional; repeatable)

A `faidx`-indexed personalized genome FASTA (e.g. from `bcftools consensus` on a
sample VCF) plus `--reference` (the unmodified genome). The variant axis compares
the splice motif at each novel junction between the reference and the haplotype(s).
**Caveat:** coordinates must align between reference and haplotype — fine for SNVs;
indels shift downstream coordinates (a documented limitation). A splice dinucleotide that
contains an N, a junction of unknown strand, or a haplotype whose contig is missing from
the FASTA is not evaluable at that junction: `created` needs one readable canonical
haplotype, `disrupted` needs every haplotype readable. Set
`--haplotype-provenance {wgs|external}` for independent genomic data (promotes the
rescue) vs `{rna_derived|unknown}` (held under the circularity firewall).

## Pangenome graph junctions (`--pangenome-junctions`, optional)

A tab-separated file (no header; `#` comment lines allowed) of splice junctions that
are present on at least one haplotype path of a pangenome graph (e.g. HPRC),
pre-extracted with `vg`/`rpvg`. PanIsoGuard does not traverse the graph itself; it
adjudicates using this junction set. Columns:

```text
chrom   intron_start   intron_end   strand   [n_haplotypes]   [haplotype_ids]
```

`intron_start`/`intron_end` are **1-based inclusive** (first/last intronic base, the
same convention as STAR `SJ.tab` and GTF introns); `strand` is `+`/`-` (`1`/`2` also
accepted); the optional 5th column is the number of supporting graph haplotypes
(default `1`, gated by `axis_pangenome.min_haplotypes`; note `1` is a single graph path,
not a population frequency — raise it for stricter support); the optional 6th column is
a comma-separated list of the haplotype/path ids carrying the junction (e.g.
`HG002#1,HG005#2`; duplicate rows for one junction union their ids;
`benchmark/pangenome/derive_junctions.py` writes it from the genotype columns of a
`vg deconstruct` VCF). When **all** of an
isoform's novel-vs-linear-reference junctions are found here — on **one common path**, if
path ids are given — the apparent novelty is consistent
with reference bias → `REFERENCE_BIAS` (mechanism `population_known`); with path ids and
no common path the call is held `UNRESOLVED` (`path_inconsistent`). This
is independent evidence only if the file was extracted from population assemblies;
PanIsoGuard cannot verify that, so the rescue is gated by
`--pangenome-provenance {population|external|sample_derived|unknown}` — `population`/
`external` promote, while the default `unknown` (and `sample_derived`) are held
`UNRESOLVED` under the same circularity firewall as the variant axis. When path ids are
given on some rows only, the rows that have them must share a path; a row without ids
cannot hide a contradiction between the others. Checked on real HPRC v1.1 chr22 graph
junctions ([benchmark/pangenome](../benchmark/pangenome), CTest `integration_pangenome_gate1`).

## Rules TOML (`--config`, optional)

[`config/rules.default.toml`](../config/rules.default.toml) lists every key with its default;
a missing key keeps the default. A key that is present must have the right type and a value
in its range, or the run stops with the key name and the range: integer counts (`sj_min_uniq_reads`,
`consensus_min_callers`, `min_spanning`, `min_full_reads`, `min_haplotypes`) are 1 or more;
`min_mapq` is 0–255; `junction_window_bp` is 0–10000; `max_perc_A_downstream_TTS` is 0–100;
the three mapping fraction gates are 0–2, where a value above 1 can never fire and so
switches that gate off.
