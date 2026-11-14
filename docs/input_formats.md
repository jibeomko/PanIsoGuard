# Input formats

All internal coordinates are **0-based, half-open** introns `[start, end)`. The
readers convert each source so the same intron yields the same coordinates, which
is what makes fingerprint matching and short-read corroboration exact:

- STAR `SJ.tab` (1-based inclusive `s..e`) → `{s-1, e}`
- GTF exons (1-based); intron between exon_i and exon_{i+1} → `{exon_i.end, exon_{i+1}.start-1}`
- BED12 blocks (0-based); intron between blocks → `{block_i.end, block_{i+1}.start}`

## SQANTI3 classification (`--classification`, required unless `--no-sqanti3`)

Tab-separated, with a header. The only mandatory column is `isoform`; everything
else is read by name (tolerant of column reordering and appended columns).
Read (as priors, never recomputed): `chrom`, `strand`, `structural_category`, `RTS_stage`,
`all_canonical`, `perc_A_downstream_TTS` and, when present, `polyA_motif_found`. Every
other column is ignored. Category strings are the
SQANTI3 spelling (`full-splice_match`, `incomplete-splice_match`,
`novel_in_catalog`, `novel_not_in_catalog`, …).

`--no-sqanti3` replaces the file. The structural category is derived from the intron chain
against `--ref-gtf`: an annotated chain is `full-splice_match`, a contiguous part of one
`incomplete-splice_match`, a new combination of annotated introns `novel_in_catalog`, a
chain with a novel intron in an annotated locus on its strand `novel_not_in_catalog`, one
whose locus is annotated only on the other strand `antisense`, one in an unannotated region
`intergenic`; a mono-exonic isoform is `mono-exon` and not adjudicated. With `--reference`,
`all_canonical` reads every intron's dinucleotides (GT-AG, GC-AG, AT-AC) and
`perc_A_downstream_TTS` the A content of the 20 bases after the 3′ end on the transcript
strand, SQANTI3's definitions; without it both are NA. `RTS_stage` and `polyA_motif_found`
are NA. The provenance log says the rows were computed. On the chr22 five-caller set the
categories equal SQANTI3's for 98.6 % of the multi-exon isoforms; the rest are
retained-intron sub-chains that SQANTI3 files under `incomplete-splice_match` and this mode
under `novel_in_catalog`, so they are adjudicated (benchmark/no_sqanti3).

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

A coordinate-sorted, indexed (`.bai`/`.csi`) genome alignment of the long reads
(splice-aware: `N` CIGAR ops at introns). CRAM additionally needs `--reference`.
Read-level features are computed only over reads whose `N` op exactly matches a
novel junction. A read whose aligner recorded the transcript strand (`XS:A` from STAR or
HISAT2, `ts:A` from minimap2, which is relative to the read and flipped for reverse
alignments) counts only for junctions and chains on that strand; a read without the tag
counts for either. A read error in the file (for example a truncated BAM) stops the run
instead of counting as a junction with no reads. Note: a per-sample BAM gives
sample-specific coverage; for a merged isoform set, use a merged BAM for full coverage.

`--validation-bam PATH` (repeatable) names libraries the caller did not use: other
replicates, the other platform, a second run. Whole-chain reads are then counted in each of
them separately; the isoform's count is the largest, a library supports the chain when it
reaches `chain_min_full_reads` on its own, and the JSONL records the count per library
(`chain_full_reads_by_library`) next to the caller's own BAM count
(`chain_n_full_reads_own`, when `--bam` is also given). The same file given twice, or the
`--bam` file, is refused. The junction-level mapping features still come from `--bam`.

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

## Sample VCF/BCF (`--haplotype-vcf`, optional; repeatable)

An indexed VCF (`bgzip` + `tabix`) or BCF (`bcftools index`) with a genotype column, plus
`--reference`. The alleles of the sample's genotype are the haplotypes: the first allele
builds haplotype 1 and the second haplotype 2, as `bcftools consensus -H 1` and `-H 2` do, so
a run from the VCF gives the verdicts of a run from those FASTAs. Several files add their
haplotypes in order (two haploid assembly-derived VCFs give two haplotypes). With more than
one sample in the file, name one with `--vcf-sample`. Only same-length substitutions (SNVs,
MNPs) are applied at the splice dinucleotides; indels and symbolic alleles are skipped and
counted (`haplotype_vcf_non_snv_skipped_at_motifs` in the provenance log). A missing allele
(`.`) or a genotype with fewer alleles than the sample's ploidy leaves the motif base unknown
(N) on that haplotype, so the motif is not evaluable there, never reference
(`haplotype_vcf_missing_alleles_at_motifs`). The ploidy is the largest genotype in the file,
read once when the file is opened. An unphased heterozygous genotype is applied by allele
order and counted (`haplotype_vcf_unphased_het_applied_at_motifs`); the junction table's
`variant_phase` column says `phased` (phased or homozygous), `unphased`, or `NA`. One created
junction is rescued whichever haplotype carries the allele. Two records inside one motif (a
donor and an acceptor variant, or two bases of one dinucleotide) that place alleles by
haplotype index, heterozygous ones and partly missing ones such as `1|.` or `0|.`, are combined
only when they are phased into one phase set (FORMAT `PS`; phased genotypes without `PS` form
one default set; a `/` before a missing allele is unphased); otherwise no haplotype's motif at
that junction is known and the junction table says `unlinked`. A homozygous record places the
same base on every haplotype and needs no linkage. Two or more created junctions that rest on heterozygous alleles are
rescued only when those alleles share one phase set: unphased alleles, or alleles in different
sets, do not establish a common haplotype and the isoform is held with reason
`haplotype_phase_unknown` (`haplotype_phase_unknown` in the JSONL). In both cases the call is
the same whichever haplotype index the alleles are written on. A record without a `GT` field
leaves its bases unknown on every haplotype. A contig missing from the VCF header is
not evaluable; a contig present without records is reference. `--haplotype-vcf` and
`--reference-haplotype` are not combined in one run.

## Pangenome graph junctions (`--pangenome-junctions`, optional)

A tab-separated file (no header; `#` comment lines allowed) of splice junctions that
are present on at least one haplotype path of a pangenome graph (e.g. HPRC),
pre-extracted with `vg`/`rpvg`. PanIsoGuard does not traverse the graph itself; it
adjudicates using this junction set. A whole-genome GRCh38 set from HPRC v1.1 (deletions
≥ 50 bp, 71,435 junctions) is attached to each release (see the README); `population`
provenance applies to any sample outside the graph's 47 assemblies. Columns:

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
`vg deconstruct` VCF). When all of an
isoform's novel-vs-linear-reference junctions are found here — on one common path, if
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
in its range, a section must be a table, and a key PanIsoGuard does not read (a misspelling)
is an error, or the run stops with the key name and the range: integer counts (`sj_min_uniq_reads`,
`consensus_min_callers`, `min_spanning`, `min_full_reads`, `min_haplotypes`) are 1 or more;
`min_mapq` is 0–255; `junction_window_bp` is 0–10000; `max_perc_A_downstream_TTS` is 0–100;
the three mapping fraction gates are 0–2, where a value above 1 can never fire and so
switches that gate off.
