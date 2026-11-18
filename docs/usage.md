# Usage and output reference

[Back to the README](../README.md#inputs--outputs).

`adjudicate` is the main command. Provide:

- `--classification`: the SQANTI3 classification file, or `--no-sqanti3` to compute
  structural categories from the isoform chains and reference annotation.
- `--isoforms-gtf` or `--isoforms-bed`: the caller's isoforms, with matching isoform IDs
  when using a classification file.
- `--out-prefix`: where to write the results.

With `--no-sqanti3`, `--ref-gtf` is required. Add `--reference` to evaluate splice motifs
and downstream A content; SQANTI3's RT-switch flag is not evaluated in this mode.

When using SQANTI3 output, also give `--ref-gtf`, the annotation you gave SQANTI3. Without it no junction can be
recognized as novel, and every novel isoform ends up `UNRESOLVED`. All other inputs are
optional; each one adds one more kind of evidence (see [Evidence it uses](#evidence-it-uses)).
PanIsoGuard accepts caller output in the documented GTF or BED12 formats. The file names
below are placeholders for your own files.

**Minimal**: SQANTI3 output only. With no read evidence, novel isoforms can be flagged or held
for review, but never confirmed:

```bash
panisoguard adjudicate \
  --classification classification.txt \
  --isoforms-bed   isoforms.bed \
  --ref-gtf        annotation.gtf \
  --out-prefix     out/sample
```

**Recommended**: add the short-read junctions from STAR (`--sj-tab`) and the long-read
alignments (`--bam`). These two are what let PanIsoGuard confirm or reject a novel isoform.
`--reference`, the genome FASTA, is needed for CRAM input and for the haplotype check
below; it also enables motif and downstream A-content checks in `--no-sqanti3` mode.
The indel-near check is measured in excess of each read's own indel rate, so ONT
reads do not trip it by their chemistry alone; the thresholds were set on PacBio HiFi and hold
on ONT reads as they are (`--platform ont`, which turned the indel-near check off, is deprecated
and changes no verdict):

```bash
panisoguard adjudicate \
  --classification classification.txt \
  --isoforms-bed   isoforms.bed \
  --ref-gtf        annotation.gtf \
  --sj-tab         SJ.out.tab \
  --bam            aligned.bam \
  --reference      genome.fa \
  --out-prefix     out/sample
```

**Reference-bias check with the sample's haplotypes**: a splice site can look non-canonical on
the reference genome (for example GT-AC) but be canonical (GT-AG) in the sample, because of a
variant. Give the sample's genotypes as an indexed VCF or BCF with `--haplotype-vcf` (the two
alleles of each genotype are the two haplotypes; several files add their haplotypes), or the
haplotype sequences with `--reference-haplotype` (once per haplotype), and PanIsoGuard reports
such isoforms as reference bias instead of artifacts. Say where the haplotypes came from with
`--haplotype-provenance`. Only `wgs` or `external` (built
from DNA sequencing or another independent source) can lead to the reference-bias verdict. With
`rna_derived` or `unknown` the isoform is held as `UNRESOLVED`, because haplotypes built from the
same RNA reads cannot independently explain those reads. This guard is called the circularity
firewall. This is a genomic explanation for the splice motif; it does not independently
validate the full transcript.

```bash
panisoguard adjudicate \
  --classification classification.txt \
  --isoforms-bed   isoforms.bed \
  --ref-gtf        annotation.gtf \
  --reference      genome.fa \
  --haplotype-vcf  sample.phased.vcf.gz \
  --haplotype-provenance wgs \
  --out-prefix     out/sample
```

From the VCF only same-length substitutions reach the splice motifs; indels are skipped and
counted in the provenance log. A missing allele (`.`, or a genotype shorter than the sample's
ploidy, or a record with no `GT`) leaves that haplotype's base unknown, so the motif there is
not evaluable rather than reference. Two variants inside one motif that place alleles by
haplotype index (heterozygous or partly missing) are combined only when they are phased into
one phase set. An unphased heterozygous genotype is applied by allele order, as `bcftools consensus
-H 1` and `-H 2` would, and the junction table marks it `unphased`. A rescue that needs two or
more junctions on one haplotype requires their alleles to be phased in one phase set (`PS`);
when they are unphased or in different sets the call is held (`haplotype_phase_unknown`), since
the file does not say they share a haplotype. Haplotype FASTAs give the same verdicts when they
are built from SNVs only; a consensus that also applied indels shifts every downstream coordinate
and reads the wrong bases, and a consensus cannot carry phase sets, so with FASTAs the
across-junction claim rests on how they were built.

**Reference-bias check with a pangenome** (tested on HPRC v1.1 chr22): instead of the sample's
own haplotypes, give splice junctions taken from a pangenome graph such as HPRC (extracted
beforehand, for example with vg/rpvg). The rescue requires support for every novel junction;
when path identifiers are available, they must support a consistent path (see
[the input format](input_formats.md#pangenome-graph-junctions---pangenome-junctions-optional)). This is evidence for a
possible reference-bias explanation, not proof that the sample carries the path or expresses
the isoform. The evidence is treated as independent only with a declared independent source, so
`--pangenome-provenance` works like `--haplotype-provenance`: `population` or `external` can
lead to the verdict, while the default `unknown` (or `sample_derived`) holds the isoform as
`UNRESOLVED`.

```bash
panisoguard adjudicate \
  --classification classification.txt \
  --isoforms-bed   isoforms.bed \
  --ref-gtf        annotation.gtf \
  --pangenome-junctions pangenome_junctions.tsv \
  --pangenome-provenance population \
  --out-prefix     out/sample
```

A ready-made junction set for GRCh38 ships with each release: the 71,435 splice junctions
implied by deletions of 50 bp or more in the HPRC v1.1 Minigraph-Cactus graph (47 assemblies;
`vg deconstruct` against the GRCh38 path, then `benchmark/pangenome/derive_junctions.py`; both
strands listed, with the number of carrying haplotypes). It is population data, so
`--pangenome-provenance population` applies, unless the sample is one of the 47 HPRC
assemblies.

```bash
wget -O hprc_junctions.tsv.gz https://github.com/jibeomko/PanIsoGuard/releases/download/v0.1.0/hprc_v1.1_mc_grch38_deletion_junctions_ge50bp.tsv.gz
echo "1503dda3e41e819bd68bc42b4c50241657f40bc82fb2c8514b7da1168cfa3241  hprc_junctions.tsv.gz" | sha256sum -c
gunzip hprc_junctions.tsv.gz    # -> --pangenome-junctions hprc_junctions.tsv
```

**Several callers** (optional): if you ran more than one caller, `combine` matches their
isoforms by intron chain (the same chain of introns counts as the same isoform) and writes a
table of which callers found each one. Give this table to `adjudicate` with `--caller-support`:

```bash
panisoguard combine \
  --gtf flair:flair.gtf \
  --gtf isoquant:isoquant.gtf \
  --gtf bambu:bambu.gtf \
  --ref-gtf annotation.gtf \
  --out caller_support_matrix.tsv
```

`combine` uses exact intron-chain matching. On the five-caller chr22 benchmark it reproduces
`gffcompare -i`'s caller-count distribution and precision-recall curve. Caller agreement is
one evidence axis in the verdict; see
[the comparison with merge tools](relationship_to_merge_tools.md).

## Outputs

`adjudicate` writes four files whose names start with your `--out-prefix`:

| File | What is in it |
|------|----------|
| `<prefix>.adjudicated.tsv`   | One row per isoform: the confidence class, the `verdict` (known / known_fragment / confirmed / reference_bias / unconfirmed / artifact / not_adjudicated) and `verdict_reason` (the one rule that decided it), the short-read support level, the main artifact cause (if any), novel-junction counts, and `genomic_compatibility` (the DNA side on its own: `haplotype` / `pangenome` / `partial` / `none` / `NA`, whatever the verdict) |
| `<prefix>.attribution.jsonl` | One JSON record per isoform: all the evidence that was used and the `rule_trace` (the rules that decided the verdict, in order) |
| `<prefix>.junctions.tsv`     | One row per junction of the adjudicated isoforms: which isoforms use it, whether it is in the reference, its short-read count, its splice motif on the reference and on each haplotype (`created` / `disrupted` / `none`), pangenome support and paths, and the long-read alignment fractions. A reference-bias junction inside an isoform that was not rescued (another junction unexplained) is visible here |
| `<prefix>.provenance.log`    | how the run was done: tool and rule versions, the thresholds, which evidence was used, whether any input could make the evidence circular, and how many isoforms got each class |

With `--write-kept PATH`, `adjudicate` also copies the input isoform GTF (or BED12) keeping only
the isoforms whose verdict is in `--keep-verdicts` (default `known,known_fragment,confirmed,
reference_bias`), so the filtered annotation can go straight to quantification; the provenance
log records both.

Here is what the verdicts look like for four isoforms of the notes' toy gene (real output of one
run with short reads, a long-read BAM and two haplotypes). The outlined evidence is what
decided each verdict.

![Verdict examples: a table of four novel isoforms of a toy gene. For each, the evidence PanIsoGuard checks (short reads, long reads, SQANTI3 QC, genome) is marked as supports, partly, against or nothing found, followed by the verdict and its class, with the evidence that decided it outlined: iso_B confirmed novel (junction in 12 short reads), iso_C unconfirmed (1 of 2 junctions confirmed), iso_G artifact (no short-read support at the junction and an indel beside it), iso_A reference bias (non-canonical only on the reference genome, canonical on the person's haplotype).](figures/verdict_examples.png)

## Confidence classes

Every isoform gets one of eight classes. The `verdict` column groups them into what you do
with the isoform, and `verdict_reason` names the one rule that decided it.

| Class | `verdict` | Meaning |
|---|---|---|
| `KNOWN` | known | Matches a known transcript exactly (SQANTI3 FSM). Not adjudicated. |
| `KNOWN_FRAGMENT` | known_fragment | A truncated known transcript (SQANTI3 ISM) whose whole intron chain at least 3 long reads carry, with no intra-priming sign at the 3' end (`ism_chain_support`). A known structure, not a novel isoform. |
| `NOVEL_CONFIRMED` | confirmed | Short reads confirm every novel junction, and there is no sign of an artifact. |
| `NOVEL_PROBABLE` | confirmed | Short reads confirm every novel junction but there is one artifact sign; or short reads could not be checked and either at least 2 callers found the same isoform (`caller_consensus`) or at least 3 long reads carry its whole intron chain (`chain_support`); both are methodological rather than experimental support, and an RT-switching sign does not demote them (`;rt_switch_sign`). |
| `UNCONFIRMED` | unconfirmed | Weak evidence: only some novel junctions confirmed (`partial_support`), or none but no artifact sign (`unsupported_no_sign`). Also an ISM without whole-chain read support (`ism`, `ism_intrapriming`, `ism_rts`), and, when there are no short reads, a novel isoform whose long-read alignments at the junction look wrong (`mapping_without_short_reads`). |
| `UNRESOLVED` | unconfirmed | Nothing to test: no short reads (`no_short_reads`), no coverage at the locus (`no_coverage`), no novel junction to check (`no_novel_junction`, NIC combinations of known junctions), or a reference-bias explanation held by the circularity firewall (`firewall_held`). SQANTI3 categories PanIsoGuard does not judge (genic, antisense, fusion, intergenic) are here with verdict `not_adjudicated`. |
| `REFERENCE_BIAS` | reference_bias | Haplotype (`haplotype_variant`) or pangenome (`pangenome_path`) evidence supports a reference-bias explanation for every novel junction. This does not independently confirm the full isoform; a population path also does not establish the sample's genotype. |
| `ARTIFACT` | artifact | Probably not real: short reads were checked, none supports any novel junction, the isoform's known junctions (if it has any) do have reads, and there is an artifact sign (`unsupported_<sign>`). Never given without short reads. |

The names before 0.0.6 were `HIGH_CONF_KNOWN`, `HIGH_CONF_NOVEL`, `MEDIUM_CONF_NOVEL`,
`LOW_CONF_PARTIAL`, `AMBIGUOUS`, `PAN_REF_RESCUED_FALSE_NOVEL` and `ARTIFACT`, in that order.

The classes are labels, not a scale. On SQANTI-SIM truth (chr22, 1 M simulated short-read
pairs, with the BAM), the share of genuine isoforms in each class is:

| Class | genuine | false | precision |
|---|---:|---:|---:|
| `NOVEL_CONFIRMED` | 351 | 3 | 0.992 |
| `NOVEL_PROBABLE` | 198 | 1 | 0.995 (151 of them lifted by whole-chain reads) |
| `KNOWN_FRAGMENT` | 153 | 11 | 0.933 (ISM carried end to end by the reads; not a novel call) |
| `UNRESOLVED` | 0 | 6 | 0.000 (without the BAM: 160 / 7, 0.958) |
| `UNCONFIRMED` | 23 | 36 | 0.390 (without the BAM: 169 / 61, 0.735, mostly ISM) |
| `ARTIFACT` | 5 | 83 | 0.057 |

Without the BAM, `UNRESOLVED` is mostly genuine: NIC isoforms with no novel junction to test.
With the BAM (0.0.7), the ones whose whole intron chain is carried by at least 3 long reads
move to `NOVEL_PROBABLE`, the ISM isoforms carried end to end move to `KNOWN_FRAGMENT`, and
`UNRESOLVED` keeps the 6 false ones. What stays `UNCONFIRMED` is well below the 0.83 base
rate. The same table at every short-read depth:
`benchmark/sqanti3_filter_h2h/per_class.py`.

The signs of an artifact are checked in this order, and the first one found is reported:
problems in the long-read alignments at the junction (`mapping_or_repeat`; needs at least 3
spanning reads), a non-canonical splice site (`noncanonical`), reverse-transcriptase template
switching (`rt_switch`), and an A-rich stretch right after the 3' end with no polyA motif
before it (`degradation`). The reference-bias checks come before all of these. Details:
[docs/decision_engine.md](decision_engine.md).

## Evidence it uses

| Evidence | Input | Needed? | What PanIsoGuard checks |
|---|---|---|---|
| SQANTI3 QC | `--classification`, or `--no-sqanti3` | Required (one of the two) | the structural category, and three artifact signs from SQANTI3's columns: non-canonical splice site, RT template switching, A-rich 3' end. With `--no-sqanti3` the category, the motif and the A content are computed from the chains, the annotation and `--reference`; the RT-switch flag is not evaluable |
| Short reads | `--sj-tab` (STAR `SJ.out.tab`) | recommended | whether short reads confirm each novel junction; by default that takes at least 3 uniquely mapped reads and a canonical splice motif |
| Long reads | `--bam` | recommended | whether the reads across each novel junction look misaligned: by default, more than half of them share one problem (low mapping quality, a supplementary alignment, or an indel next to the junction beyond what the reads' own indel rate predicts). The defaults were set on PacBio HiFi. Reads whose aligner put them on the other strand (`XS:A`, `ts:A`) are not counted. |
| Caller agreement | `--caller-support` (from `combine`) | optional | how many callers found the same intron chain; used only when short-read support cannot be checked |
| Whole-chain reads | `--bam` (same file as above) or `--validation-bam` (repeatable: other replicates or platforms, each counted on its own) | With the BAM | how many long reads carry the isoform's entire intron chain; used only when short-read support cannot be checked; it can yield `NOVEL_PROBABLE` (verdict `confirmed`), but not `NOVEL_CONFIRMED`. With `--validation-bam` (a library the caller did not use) the count comes from there and the reason carries `;independent_library` |
| Sample's genome | `--reference` with `--haplotype-vcf` or `--reference-haplotype` | optional | whether a variant makes a non-canonical splice site canonical in the sample (reference bias) |
| Pangenome | `--pangenome-junctions` | optional | whether every novel junction is found in a pangenome (reference bias; tested on HPRC v1.1 chr22) |

## All commands

| Command | What it does |
|---------|---------|
| `adjudicate` | gives each isoform a confidence class, the likely artifact cause, and the reasons (4 output files) |
| `combine`    | matches several callers' isoforms by intron chain and writes which callers found each one |
| `decoys`     | measures, on your own BAM, how often the whole-chain read support would lift a chain the annotation does not know (recombined known transcripts), as an upper bound on its false-lift rate |
| `version`    | shows the version and the linked htslib |

Run `panisoguard <command> --help` for all options. Input file formats:
[docs/input_formats.md](input_formats.md).

### Decoys

```bash
# how often >= 3 whole-chain reads carry a chain the annotation does not know, in this BAM
panisoguard decoys    --ref-gtf ref.gtf --bam aln.bam --out decoys/sample
```

`decoys` joins the 5' introns of one annotated transcript to the 3' introns of another of the
same gene at a shared intron (a novel combination of known junctions, which no short-read test
can reach) and counts the reads that carry each such chain end to end. The fraction with at
least 3 is what the `chain_support` lift would do to chains that are probably wrong, in your
data. In simulation, where nothing expressed is unannotated, wrong combinations reach 3 reads
in 0-0.1 % of cases (and the deleted genuine transcripts in 100 %); on real PacBio and Kinnex
data a third to a half of the recombinants of two expressed isoforms do, with as many reads as
their parents, which is what unannotated combinatorial isoforms look like rather than chimeras
([benchmark/chain_support](../benchmark/chain_support)). The number is an upper bound on the
false-lift rate and a measure of how much of your transcriptome the annotation misses.

## A worked locus

[examples/worked_locus](../examples/worked_locus/README.md) follows one toy gene, from input files to
verdicts: a gene on a 1,700 bp chromosome with three reference transcripts and nine caller
isoforms, with every input a run can take. Adding the inputs one at a time changes the verdict of
`iso_A`, the isoform with one novel junction, from `UNRESOLVED` (nothing to test it with) to
`ARTIFACT` (short reads do not support its non-canonical junction) to `REFERENCE_BIAS` (the
junction is canonical on the sample's haplotype, with independent provenance). Two scripts check the
example and run in CI: `check_worked_locus.py` recomputes every verdict with the Python
standard library and compares it with the binary, and `rules_replay.py` is the whole rulebook of
`RuleEngine::evaluate()` in 50 lines of Python, replayed against the binary's `attribution.jsonl`.
