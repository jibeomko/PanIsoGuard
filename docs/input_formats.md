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
Consumed as priors (never recomputed): `chrom`, `strand`, `structural_category`,
`associated_gene/transcript`, `subcategory`, `RTS_stage`, `all_canonical`,
`perc_A_downstream_TTS`, `n_indels_junc`, `dist_to_CAGE_peak`, `dist_to_polyA_site`,
and the filter column (`filter_result` or `filter`). Category strings are the
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
evaluable and novel isoforms fall to `AMBIGUOUS`.

## STAR SJ.tab (`--sj-tab`, optional)

9-column STAR `SJ.out.tab` / merged `combined_SJ.tab`: `chrom, intron_start,
intron_end, strand{0,1,2}, motif{0..6}, annotated, n_uniq, n_multi, max_overhang`.
A novel junction is corroborated if it matches (chrom, strand, coords) exactly with
`n_uniq ≥ sj_min_uniq_reads` and (by default) a canonical motif.

## BAM/CRAM (`--bam`, optional)

A coordinate-sorted, **indexed** (`.bai`/`.csi`) genome alignment of the long reads
(splice-aware: `N` CIGAR ops at introns). CRAM additionally needs `--reference`.
Read-level features are computed only over reads whose `N` op exactly matches a
novel junction. Note: a per-sample BAM gives sample-specific coverage; for a merged
isoform set, use a merged BAM for full coverage.

## Personalized haplotype FASTA (`--reference-haplotype`, optional; repeatable)

A `faidx`-indexed personalized genome FASTA (e.g. from `bcftools consensus` on a
sample VCF) plus `--reference` (the unmodified genome). The variant axis compares
the splice motif at each novel junction between the reference and the haplotype(s).
**Caveat:** coordinates must align between reference and haplotype — fine for SNVs;
indels shift downstream coordinates (a documented limitation). Set
`--haplotype-provenance {wgs|external}` for independent genomic data (promotes the
rescue) vs `{rna_derived|unknown}` (held under the circularity firewall).

## Pangenome graph junctions (`--pangenome-junctions`, optional)

A tab-separated file (no header; `#` comment lines allowed) of splice junctions that
are present on at least one haplotype path of a pangenome graph (e.g. HPRC),
pre-extracted with `vg`/`rpvg`. PanIsoGuard does not traverse the graph itself — it
adjudicates using this junction set (the in-process GBZ tier stays behind
`-DWITH_PANGENOME_LIB`). Columns:

```text
chrom   intron_start   intron_end   strand   [n_haplotypes]
```

`intron_start`/`intron_end` are **1-based inclusive** (first/last intronic base, the
same convention as STAR `SJ.tab` and GTF introns); `strand` is `+`/`-` (`1`/`2` also
accepted); the optional 5th column is the number of supporting graph haplotypes
(default `1`, gated by `axis_pangenome.min_haplotypes`; note `1` is a single graph path,
not a population frequency — raise it for stricter support). When **all** of an isoform's
novel-vs-linear-reference junctions are found here, the apparent novelty is consistent
with reference bias → `PAN_REF_RESCUED_FALSE_NOVEL` (mechanism `population_known`). This
is independent evidence only if the file was extracted from population assemblies;
PanIsoGuard cannot verify that, so the rescue is gated by
`--pangenome-provenance {population|external|sample_derived|unknown}` — `population`/
`external` promote, while the default `unknown` (and `sample_derived`) are held
`AMBIGUOUS` under the same circularity firewall as the variant axis. **Experimental:**
not yet validated on real graph data (GATE-1).
