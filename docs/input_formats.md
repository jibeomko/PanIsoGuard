# Input formats

All internal coordinates are **0-based, half-open** introns `[start, end)`. The
readers convert each source so the same intron yields the same coordinates, which
is what makes fingerprint matching and short-read corroboration exact:

- STAR `SJ.tab` (1-based inclusive `s..e`) gives `{s-1, e}`
- GTF exons (1-based); intron between exon_i and exon_{i+1} gives `{exon_i.end, exon_{i+1}.start-1}`
- BED12 blocks (0-based); intron between blocks gives `{block_i.end, block_{i+1}.start}`

## SQANTI3 classification (`--classification`, required unless `--no-sqanti3`)

Tab-separated, with a header. The only mandatory column is `isoform`; everything
else is read by name (tolerant of column reordering and appended columns).
Read (as priors, never recomputed): `chrom`, `strand`, `structural_category`, `RTS_stage`,
`all_canonical`, `perc_A_downstream_TTS` and, when present, `polyA_motif_found`. Every
other column is ignored. Category strings are the
SQANTI3 spelling (`full-splice_match`, `incomplete-splice_match`,
`novel_in_catalog`, `novel_not_in_catalog` and so on).

`--no-sqanti3` replaces the file. The structural category is derived from the intron chain
against `--ref-gtf`: an annotated chain is `full-splice_match`, a contiguous part of one
`incomplete-splice_match`, a new combination of annotated introns `novel_in_catalog`, a
chain with a novel intron in an annotated locus on its strand `novel_not_in_catalog`, one
whose locus is annotated only on the other strand `antisense`, one in an unannotated region
`intergenic`; a mono-exonic isoform is `mono-exon` and not adjudicated. With `--reference`,
`all_canonical` reads every intron's dinucleotides (GT-AG, GC-AG, AT-AC) and
`perc_A_downstream_TTS` the A content of the 20 bases after the 3' end on the transcript
strand, SQANTI3's definitions; without it both are NA. `RTS_stage` and `polyA_motif_found`
are NA. The provenance log says the rows were computed. On the chr22 five-caller set the
categories equal SQANTI3's for 98.6 % of the multi-exon isoforms; the rest are
retained-intron sub-chains that SQANTI3 files under `incomplete-splice_match` and this mode
under `novel_in_catalog`, so they are adjudicated (benchmark/no_sqanti3).

## Caller isoforms (`--isoforms-bed` or `--isoforms-gtf`, required)

The structure of the called isoforms, as BED12 or GTF. **Transcript/isoform IDs
must match the SQANTI3 `isoform` IDs**: this is the join key between QC and
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
A haplotype FASTA is read at the reference's coordinates, so it can carry substitutions
only: build it from SNVs. A contig whose length differs from the reference's is not read
(that haplotype is not evaluable there; a warning names the file and the contig, and the
provenance log counts them in `reference_haplotype_contig_length_mismatches`); give
insertions and deletions with `--haplotype-vcf`. A splice dinucleotide that
contains an N, a junction of unknown strand, or a haplotype whose contig is missing from
the FASTA is not evaluable at that junction: `created` needs one readable canonical
haplotype, `disrupted` needs every haplotype readable. Set
`--haplotype-provenance {wgs|external}` for independent genomic data (promotes the
rescue) vs `{rna_derived|unknown}` (held under the circularity firewall).

## Sample VCF/BCF (`--haplotype-vcf`, optional; repeatable)

An indexed VCF (`bgzip` + `tabix`) or BCF (`bcftools index`) with a genotype column, plus
`--reference`. The alleles of the sample's genotype are the haplotypes: the first allele
builds haplotype 1 and the second haplotype 2, as `bcftools consensus -H 1` and `-H 2` do, so
for a file of substitutions a run from the VCF gives the verdicts of a run from those FASTAs.
Several files add their haplotypes in order (two haploid assembly-derived VCFs give two
haplotypes). With more than one sample in the file, name one with `--vcf-sample`. The ploidy
is the largest genotype in the file, read once when the file is opened. A contig missing from
the VCF header is not evaluable; a contig present without records is reference.
`--haplotype-vcf` and `--reference-haplotype` are not combined in one run.

- Opening reads the whole file once, every sample's columns included. From a panel of many
  samples, write the sample's own file first (`bcftools view -s NAME -c 1`).
- Every record is applied, whatever its FILTER column says, as `bcftools consensus` does. To
  leave records out, filter the file first (`bcftools view -f PASS`).
- Contig names are those of `--reference` and the isoform file. A VCF or haplotype FASTA that
  names none of the isoforms' contigs (`1` for `chr1`) gets a warning, and the variant axis is
  not evaluable with it.
- A file that cannot be read whole is refused when it is opened: a BGZF file without its
  end-of-file block (cut short) and a record htslib cannot parse. Otherwise the records after
  the break would read as reference.

**What is applied.** Every allele that has a sequence: substitutions (SNVs, MNPs), insertions,
deletions and other replacements. An insertion or deletion is moved towards the intron as far
as the reference sequence allows, so the exon is read as the reference has it wherever that is
possible, and a haplotype's motif is the two bases next to the last exonic base on each side.

- An insertion or deletion inside or next to the motif changes it as a substitution does
  (`haplotype_vcf_indels_applied_at_motifs` in the provenance log).
- An insertion that ends up exactly between the exon and the intron can be intronic (the
  junction reads do not carry it) or part of the exon (they do); the file cannot tell. With
  `--bam` the junction reads decide: with at least `min_spanning` of them, an insertion of that
  length within 3 bp plus its length on the exon side in at most 20 % of the reads makes it
  intronic, in at least 80 % exonic. Between the two, and without `--bam`, that haplotype's
  motif is not evaluable.
- A deletion equal to the intron (wherever the reference lets the same deletion be written)
  means the haplotype lacks those bases: its reads align with a gap there without a splice.
  The junction table says `deleted` and gives the haplotype's motif as `DEL`; the junction
  counts as explained, as a created motif does, and the verdict reason is `haplotype_deletion`
  (`haplotype_vcf_deletions_equal_to_junctions`).
- A symbolic `<DEL>` with `END` or `SVLEN` is a deletion. An allele without a sequence
  (another symbolic allele, a breakend, a `<DEL>` without an end, a `*` whose deletion the
  file does not carry) leaves its bases unknown and is counted
  (`haplotype_vcf_unapplied_alleles_at_motifs`); a `*` inside a deletion or replacement on the
  same haplotype is that deletion.

**What stays unknown.** An unknown base reads N: the motif is not evaluable on that haplotype,
never reference, and a base once unknown is not restored by another record.

- A missing allele (`.`), a genotype with fewer alleles than the sample's ploidy, an allele
  index the record does not define, or a record without `GT`: the record's REF span is unknown
  on the haplotypes concerned (`haplotype_vcf_missing_alleles_at_motifs`).
- A REF that is not the reference sequence: the record describes another genome, and its span
  is unknown on every haplotype (`haplotype_vcf_ref_mismatches_at_motifs`).
- Two records that give one haplotype different content at one base. An allele states the
  sequence over its whole REF span, so the unchanged bases of an MNP count; the spans of both
  records are unknown, in whichever order the rows come. Two substitutions that agree are one
  statement. Two different insertions at one point leave that point and the base on either
  side unknown; the same insertion written twice is one statement. A site written twice with
  the same alleles and another genotype is unknown on every haplotype; an identical repeat is
  one record (`haplotype_vcf_conflicting_records_at_motifs`).

**Phase.** An unphased heterozygous genotype is applied by allele order and counted
(`haplotype_vcf_unphased_het_applied_at_motifs`); its allele order is not part of what it
states (`0/1` and `1/0` are the same record). The junction table's `variant_phase` column says
`phased` (phased or homozygous), `unphased`, `unlinked` or `NA`. One explained junction is
rescued whichever haplotype carries the allele. The file links heterozygous records, and partly
missing ones such as `1|.`, only inside one phase set (FORMAT `PS`; phased genotypes without
`PS` form one default set; an unphased record is a group of its own; a `/` before a missing
allele is unphased). `PS` is read as the Integer the VCF specification defines; where the
header declares it a String (some call sets name their sets, such as `PATMAT`), each distinct
value is a set, and a `PS` of another type is refused. At each junction every placement of the
groups found within 64 bp of either boundary is tried (up to six groups; more are not tried and
the junction is `unlinked`). A group matters when turning it over changes what a haplotype
carries under some placement of the others. With two or more groups that matter, the
per-haplotype bases would be a composition by allele order, so no haplotype's motif at that
junction is known (`unlinked`). A heterozygous record that changes nothing under any placement
does not count, and a homozygous record needs no linkage. Two or more explained junctions that
rest on heterozygous alleles are rescued only when those alleles share one phase set: unphased
alleles, or alleles in different sets, do not establish a common haplotype and the isoform is
held with reason `haplotype_phase_unknown`. In every case the call is the same whichever
haplotype index the alleles are written on, and however an insertion or deletion is placed
inside a repeat.

## Junction reads and the haplotype (`--bam` with either haplotype input)

A junction that only some of the sample's haplotypes explain is used by reads of those
haplotypes. For such a junction the primary alignments that splice exactly there (the first
200) are compared with the sample's other heterozygous single-base sites they cover: sites of
the same phase set for a VCF, positions where the haplotype FASTAs differ otherwise. A read is
assigned to the haplotypes it agrees with at `phase_margin` (2) more sites than with the
others, so one sequencing error cannot place it. The junction table reports `read_phase` with
the two counts (`read_phase_n_motif`, `read_phase_n_other`):

| `read_phase` | meaning |
|---|---|
| `consistent` | at least `phase_min_reads` (3) reads assigned, all to a haplotype that explains the junction |
| `contradicted` | at least `phase_min_reads` assigned, none to a haplotype that explains it: the junction is not counted as explained |
| `mixed` | reads of both kinds |
| `uninformative` | fewer than `phase_min_reads` assigned (no heterozygous site covered, or too few reads) |
| `NA` | every haplotype explains the junction or none does, a haplotype is unknown, or no `--bam` |

Only `contradicted` changes a verdict: the isoform then goes through the normal rules, and its
rule trace says so. Both thresholds are in `[axis_variant]` of the rules TOML.

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
not a population frequency; raise it for stricter support); the optional 6th column is
a comma-separated list of the haplotype/path ids carrying the junction (e.g.
`HG002#1,HG005#2`; duplicate rows for one junction union their ids;
`benchmark/pangenome/derive_junctions.py` writes it from the genotype columns of a
`vg deconstruct` VCF). When all of an
isoform's novel-vs-linear-reference junctions are found here (on one common path, if
path ids are given), the apparent novelty is consistent
with reference bias, so `REFERENCE_BIAS` (mechanism `population_known`); with path ids and
no common path the call is held `UNRESOLVED` (`path_inconsistent`). This
is independent evidence only if the file was extracted from population assemblies;
PanIsoGuard cannot verify that, so the rescue is gated by
`--pangenome-provenance {population|external|sample_derived|unknown}`: `population`/
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
`consensus_min_callers`, `min_spanning`, `min_full_reads`, `min_haplotypes`, `phase_margin`,
`phase_min_reads`) are 1 or more; `min_fraction` is 0-1;
`min_mapq` is 0-255; `junction_window_bp` is 0-10000; `max_perc_A_downstream_TTS` is 0-100;
the three mapping fraction gates are 0-2, where a value above 1 can never fire and so
switches that gate off.
