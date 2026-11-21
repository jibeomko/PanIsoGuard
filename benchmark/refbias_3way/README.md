# Reference bias three ways: GRCh38 pipeline vs personalized-genome pipeline vs PanIsoGuard

Status: run on HG03516 and HG02717 (arms A, B1, B2, C and, since 0.1.1, C with the assembly's
calls as VCFs and with the 1000 Genomes short-read calls of the same people). Numbers:
[Results](#results-hg03516-and-hg02717) below and
[../results/refbias_3way/metrics.json](../results/refbias_3way/metrics.json).

## Question

PanIsoGuard's reference-bias rescue is so far validated against its own criterion (a
splice-site SNV that is non-canonical on GRCh38 and canonical on the person's haplotype)
plus one non-circular phasing check ([refbias_phasing](../refbias_phasing): 27 of 28
heterozygous junctions). The only head-to-head table ([sqanti3_filter_h2h](../sqanti3_filter_h2h))
is on SQANTI-SIM, which is simulated from GRCh38 and therefore contains no reference bias
by construction, so the rescue cannot show there.

This benchmark gives the rescue an independent standard: the novel isoforms that a
full personalized-genome re-analysis of the same reads retains or discards. It asks

1. how much of what the personalized pipeline changes, the post-hoc rescue recovers
   (sensitivity, precision), and at what fraction of the compute;
2. which reference-bias mechanisms the rescue cannot see by construction (indels, SVs,
   mapping shifts, motif-*disrupting* variants), and how large each is.

Expected shape of the result: with SNV-consensus FASTAs (arm C) the rescue should recover most
of the **splice-dinucleotide SNV** class (that is what a FASTA in reference coordinates can
carry) and none of the rest; with the assembly's calls as VCFs (arm C-VCF, 0.1.1) it also reads
insertions and deletions at the motif and deletions equal to a junction. Arm B1 below isolates the first;
the difference B2 - B1 measures the second. A sensitivity well below 100 % against the full
personalized pipeline is the expected, reportable scope of the method, not a failure.

## Arms

All arms use the same reads (HG03516 Kinnex Iso-Seq flnc, 11.34 M aligned reads), the same
caller (IsoQuant 3.13, `pacbio_ccs`), the same annotation (GENCODE v49) and the same QC
(SQANTI3 6.0.1 + its rules filter). Only the genome the reads are aligned to differs.

| Arm | Genome | Annotation | Extra steps | What it measures |
|---|---|---|---|---|
| A | GRCh38 primary | GENCODE v49 | none | the reference-only pipeline (today's default) |
| B1 | SNV-personalized diploid: `mat_consensus.fa` + `pat_consensus.fa` (GRCh38 backbone with the haplotype's assembly SNVs applied; GRCh38 coordinates) | GENCODE v49, contigs renamed per haplotype | none (same coordinates) | SNV-mediated reference bias only, the component PanIsoGuard models |
| B2 | HPRC R2 phased assembly, diploid: `mat.fa.gz` + `pat.fa.gz` (72 + 74 contigs, chromosome-scale) | GENCODE v49 lifted per haplotype with liftoff | lift junctions back to GRCh38 (`paftools.js liftover` on the asm5 PAF) + read-level cross-check | the full personalized pipeline: SNVs, indels, SVs, assembly-specific sequence |
| C | GRCh38 (arm A's alignment and models) | GENCODE v49 | `panisoguard adjudicate --bam rna.aln.bam --reference GRCh38 --reference-haplotype mat_consensus.fa --reference-haplotype pat_consensus.fa --haplotype-provenance wgs` (no short reads exist for this sample, so the short-read axis is off) | PanIsoGuard post hoc on the reference-only pipeline |
| C-1kGP | as C | GENCODE v49 | the same with one diploid `--haplotype-vcf`: the person's genotypes in the 1000 Genomes high-coverage phased panel (30x Illumina; SNVs, indels, SVs), as phased and with the phase removed | the post-hoc run with the variant calls a short-read WGS study has, independent of the assembly |
| C-VCF | as C | GENCODE v49 | the same with `--haplotype-vcf mat.hap.vcf.gz --haplotype-vcf pat.hap.vcf.gz` instead of the FASTAs: every `paftools.js call` record of each assembly haplotype against GRCh38 (3.8-4.0 M per haplotype: substitutions, insertions, deletions), one haploid VCF per haplotype | the same post-hoc run with the indels the consensus FASTAs leave out |

Arm A alignment and IsoQuant come from [hg03516_refbias](../hg03516_refbias) (68,011
models, 37,568 novel); `A_sqanti` runs SQANTI3 and its rules filter on all of them.

Diploid rather than per-haplotype alignment for B1/B2 because that is what a personalized
pipeline does (each read picks its haplotype). Fallback if IsoQuant or SQANTI3 object to a
6.2 Gb diploid reference: run per haplotype and take the union (a junction counts as
canonical/known if it is on either haplotype, which is also the rescue's rule).

## Unit and labels

The unit is a junction (intron, GRCh38 coordinates, strand) used by at least one
caller-reported novel model (NIC/NNC) in arm A: `J_A` (98,672 for HG03516). Isoform-level
counts (what a user keeps after the SQANTI3 filter in each arm) are reported alongside.

For each `j` in `J_A`, arm A gives its status on GRCh38 (`known`, `novel-canonical`,
`novel-noncanonical`) and the filter outcome of the models that use it.

For arm B (B1 and B2), the label comes from the reads, not from coordinate liftover
alone, because IsoQuant model coordinates can be displaced from where the reads splice
(the P2RX5 case in [hg03516_refbias](../hg03516_refbias)). Take the reads that span `j`
exactly in arm A; find where the same reads splice on the personal genome; lift that
position back to GRCh38 (B2 only; B1 is already in GRCh38 coordinates). Each read falls in
one class, and the junction takes the majority class:

| class | meaning |
|---|---|
| `same-canonical` | Same junction, now canonical in the person: the non-canonical flag was reference bias |
| `same-noncanonical` | Same junction, still non-canonical: not reference bias |
| `shifted-known` | the reads splice at a different, GENCODE-annotated junction in the person: the novelty was a mapping artifact of reference bias (indel, SV) |
| `shifted-novel` | different junction, still novel: ambiguous |
| `unspliced` | the reads align through the position without an intron in the person (B2: mostly short introns that minimap2 renders as a deletion on the assembly, see Results) |
| `no-reads` | No read spans the junction exactly in arm A (model coordinate displaced from the reads) |
| `unassessable` | read unaligned, or region not lifted |

Reference-bias truth: `T = same-canonical (non-canonical on GRCh38) + shifted-known +
unspliced`, restricted to `j` that arm A calls novel; `T_core` leaves `unspliced` out (see
Results for why). `T_B1` (SNV only) and `T_B2` (full).

PanIsoGuard's sets, both under `wgs` provenance: `R_rule` = novel junctions that meet the
rescue's own criterion (non-canonical on GRCh38, canonical on a consensus haplotype), i.e.
what the variant axis sees; `R` = novel junctions inside isoforms that got the
`REFERENCE_BIAS` verdict (the verdict needs every novel junction of the isoform explained,
and only NIC/NNC isoforms are adjudicated). `R_unknown` (provenance `unknown`) must be
empty: the firewall check. For arm C-VCF the sets are read from the tool's own junction table:
`V_rule` = novel junctions it reports as `created` or `deleted` and that the junction reads do
not contradict (`read_phase`), `V` = those inside `REFERENCE_BIAS` isoforms.

## Metrics

- sensitivity `|R and T| / |T|` and precision `|R and T| / |R|`, against `T_B1` and `T_B2`;
  bootstrap CIs over junctions.
- `T_B2 \ R` by mechanism: dinucleotide SNV missed (implementation), SNV outside the
  dinucleotide, indel at the site, SV / assembly-specific exon, alignment shift. This is the
  scope statement.
- `R \ T_B2`: every case inspected by hand (candidate false rescues).
- the blind spot in the other direction: arm-A `novel-canonical` junctions that become
  non-canonical or absent in the person (the scan's DISRUPTED class). PanIsoGuard does not
  model these; report their count.
- isoform level: novel models kept by the SQANTI3 rules filter in A, in B1, in B2, and kept
  by PanIsoGuard in C (`NOVEL_CONFIRMED` / `NOVEL_PROBABLE` + `REFERENCE_BIAS`); in
  particular "discarded in A, kept in B2" vs "rescued in C".
- cost per arm: wall time, CPU-hours, peak RSS, disk; C counted as adjudicate only (arm A
  is sunk cost for both A and C).

## Inputs on disk (HG03516, `/mnt/Data/PanIsoGuard_bench/pig_hg03516`, 38 GB)

| file | what | made by |
|---|---|---|
| `flnc.bam` (13.5 GB) | Kinnex Iso-Seq flnc, HPRC R2 | `dl.sh` |
| `rna.aln.bam` | minimap2 `splice:hq` vs GRCh38; 11,342,107 primary alignments; 66 min wall / 17.4 CPU-h / 23 GB at 16 threads | `rna.sh` |
| `isoquant_out/iq/iq.transcript_models.gtf` | IsoQuant 3.13; 47 min at 32 threads | `isoquant.sh` |
| `mat.fa.gz`, `pat.fa.gz` | HPRC R2 v1.1.0 phased assembly (`HG03516#2#` and `#1#` prefixes) | `dl.sh` |
| `{mat,pat}.vcf`, `.snv.vcf.gz` | paftools.js call vs GRCh38 (3.2 M SNVs, 0.4 M 1-bp indels, 0.6 M longer indels, 18 k SVs ≥ 50 bp per haplotype) | `geno.sh` |
| `{mat,pat}_consensus.fa` | bcftools consensus, SNV only, length-matched to GRCh38 | `geno.sh` |
| `scan/`, `refbias_unique.tsv` | the 40-junction rescue subset and its per-chromosome adjudicate inputs | `scan_refbias.sh` |

HG02717 (`/mnt/Data/PanIsoGuard_bench/pig_hg02717`) has the same set (69,228 models, 24,941 novel; 46
rescue junctions). Run HG03516 first, replicate on HG02717 once the scripts are fixed.

Reference: `/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa`,
`gencode.v49.annotation.gtf`.

## Prerequisites and known snags

- liftoff: `conda create -n liftoff --override-channels -c conda-forge -c bioconda liftoff=1.6.3`
  (bioconda has 1.6.3; it pulls minimap2); `run.sh` expects it at `$CONDA/liftoff/bin/liftoff`.
- SQANTI3 (`~/miniconda3/envs/sqanti3`) works only with that env's `bin` first in
  `PATH`: the wrapper's `#!/usr/bin/env python` otherwise picks another python and fails on
  `BCBio`. `export PATH=$HOME/miniconda3/envs/sqanti3/bin:$PATH` as the h2h script does.
- The asm5 PAF (`B2_paf`) is kept under `paf/`; `paftools.js liftover <aln.paf>
  <query.bed>` lifts BED in query (assembly) coordinates to the target (GRCh38) and
  writes `contig_start_end[_t5][_t3]` names for partial lifts, which the classifier treats
  as locus-only.
- IsoQuant on the diploid references needs ~2× the GRCh38 run's memory; the machine has
  32 cores / 123 GB, so run arms sequentially, not in parallel with each other.
- Keep IsoQuant's `--complete_genedb` off for the lifted annotation unless liftoff output
  is checked to be complete (liftoff drops genes it cannot place; count them, and mark their
  regions `unassessable`).

## Results (HG03516 and HG02717)

Arm A reports 21,099 (HG03516) and 16,014 (HG02717) junctions used by novel models, 2,959
and 2,233 of them non-canonical on GRCh38. Where the same reads splice on the personal genome:

| label (B arm) | HG03516 B1 | HG03516 B2 | HG02717 B1 | HG02717 B2 |
|---|---:|---:|---:|---:|
| Same-canonical, canonical on GRCh38 (no reference bias) | 17,673 | 16,974 | 13,554 | 13,117 |
| Same-noncanonical, non-canonical on GRCh38 (not reference bias) | 2,727 | 2,503 | 2,058 | 1,877 |
| Same-canonical, non-canonical on GRCh38 (splice-dinucleotide SNV; the rescue's target) | 32 | 35 | 39 | 43 |
| shifted-known (reads use an annotated junction in the person) | 7 | 12 | 11 | 10 |
| Same-noncanonical, canonical on GRCh38 (motif *disrupted* in the person) | 23 | 30 | 7 | 6 |
| shifted-novel | 169 | 255 | 83 | 171 |
| unspliced | 49 | 499 | 21 | 271 |
| unassessable / no-reads | 399 / 20 | 771 / 20 | 223 / 18 | 501 / 18 |

`T_core` (SNV class + shifted-known) = 39 / 47 (HG03516 B1 / B2) and 50 / 53 (HG02717). B2
adds the SNV-class and shifted-known junctions that the SNV-only genome cannot show: indels
and assembly-specific sequence at the splice site.

`unspliced` is left out of `T_core`. B2's 499 (HG03516) are short (median 82 bp, 55 % ≤ 100
bp; median 4 reads): on the diploid assembly the reads run through the position without an
intron. 308 of them are `same-canonical` on the SNV genome (B1), which has GRCh38's sequence
there and so shows the same gap as GRCh38. The class holds two things: a deletion in the
person that the GRCh38 pipeline reports as an intron (221 of the 499 equal a deletion in the
assembly's own calls; the VCF arm below), and short introns that minimap2 renders as a deletion
on the assembly. `T_core` keeps to the classes whose reading does not depend on that split.
HG02717 shows the same (271 on B2, 21 on B1; 130 equal a deletion).

### Rule level and verdict level

| | HG03516 B1 | HG03516 B2 | HG02717 B1 | HG02717 B2 |
|---|---:|---:|---:|---:|
| `R_rule` (junctions meeting the rescue criterion) | 35 | 35 | 42 | 42 |
| `R_rule and T` | 31 | 31 | 37 | 38 |
| rule precision | 0.886 | 0.886 | 0.881 | 0.905 |
| rule sensitivity, SNV class | 31 / 32 = 0.969 | 31 / 35 = 0.886 | 37 / 39 = 0.949 | 38 / 43 = 0.884 |
| rule sensitivity, `T_core` | 0.795 | 0.660 | 0.740 | 0.717 |
| `R` (junctions in `REFERENCE_BIAS` isoforms) | 13 (18 isoforms) | 13 | 22 (28 isoforms) | 22 |
| `R and T` | 12 | 12 | 21 | 19 |
| verdict precision (95 % CI) | 0.923 (0.77-1.00) | 0.923 | 0.955 (0.86-1.00) | 0.864 (0.73-1.00) |
| verdict sensitivity, SNV class (95 % CI) | 0.375 (0.22-0.56) | 0.343 (0.20-0.51) | 0.538 (0.38-0.69) | 0.442 (0.30-0.58) |

The motif criterion itself finds 31 of 32 (HG03516) and 37 of 39 (HG02717) SNV-class
junctions on B1. The `R_rule` junctions outside `T` (4 and 5) are `same-noncanonical` in the
person: the haplotype with the canonical motif exists, but the reads come from the other
one (the rescue is allele-agnostic; [refbias_phasing](../refbias_phasing) found
the same 1-in-10). Shifted-known junctions (indel / SV mapping shifts) are not caught: the
FASTA route cannot carry indels (the VCF arm below reads them at the splice site, not where
they move a junction to another position).

The drop from `R_rule` to `R` (35 to 13 and 42 to 22) is not the criterion. Of the SNV-class
junctions in `T_B1` without a rescue verdict (20 and 18): 15 and 10 sit in isoforms with
another novel junction the haplotypes do not explain (the verdict needs all of them), 4 and
6 are in isoforms outside NIC/NNC (not adjudicated), 1 and 2 are criterion misses. B2 adds
3 criterion misses in each individual, the non-SNV variants. `R minus T` is one junction in
HG03516 (chr6:136,665,560-136,668,208) and one locus in HG02717 (three junctions sharing the
acceptor chr6:130,022,220), all `same-noncanonical` in the person: candidate false rescues,
1 of 13 and 1-3 of 22.

### With the assembly's calls as VCFs: insertions, deletions and the reads' haplotype (0.1.1)

Arm C-VCF, against the assembly arm B2 (B1 carries GRCh38's sequence wherever the person has
an insertion or deletion, so it cannot label these):

| | HG03516 | HG02717 |
|---|---:|---:|
| `T` on B2 (SNV class + shifted-known + unspliced) | 35 + 12 + 499 = 546 | 43 + 10 + 271 = 324 |
| junctions the haplotypes explain: `created` / `deleted` | 37 / 229 | 45 / 135 |
| of them contradicted by the junction reads | 1 | 1 |
| `V_rule` | 265 | 179 |
| `V_rule and T` | 255 | 171 |
| precision (95 % CI) | 0.962 (0.94-0.98) | 0.955 (0.92-0.98) |
| sensitivity over `T` (95 % CI) | 0.467 (0.42-0.51) | 0.528 (0.48-0.59) |
| `created`: in the SNV class / all | 33 / 36 = 0.917 | 41 / 45 = 0.911 |
| SNV-class junctions found | 33 / 35 = 0.943 | 41 / 43 = 0.953 |
| `deleted`: `unspliced` on B2 / all | 221 / 229 = 0.965 | 130 / 134 = 0.970 |
| `unspliced` junctions found | 221 / 499 = 0.443 | 130 / 271 = 0.480 |
| `V` (junctions in `REFERENCE_BIAS` isoforms) | 107 (144 isoforms) | 81 (108 isoforms) |
| `V and T` | 105 = 0.981 | 78 = 0.963 |

- **The motif class.** Two more SNV-class junctions per individual than the FASTA route
  (33 and 41 against 31 and 38 on B2): one made canonical by a 1-2 bp insertion or deletion
  inside the motif, and chr19:41,565,056-41,565,563, where both haplotypes carry a 4-bp
  insertion exactly at the acceptor and the 491 and 1,186 junction reads show that it is
  intronic. The misses left (2 and 2) have no variant at the splice site in the call set. The
  junctions `created` outside the class are `same-noncanonical` in the person (3 and 3) and
  one `shifted-novel`; none has enough assigned reads to be contradicted.
- **The reads' haplotype.** HG03516 chr3:184,709,997-184,735,943 is canonical on the maternal
  haplotype only and all 3 assigned junction reads are paternal: `contradicted`, so the
  junction is no longer counted as reference bias (its isoform is `genic`, a category that is
  not adjudicated, so no verdict changes). That is the one junction
  [refbias_phasing](../refbias_phasing) found contradicted, an alignment artifact, and
  the precision of `created` rises from 0.886 to 0.917. Of the 37 and 45 `created` junctions,
  13 and 16 are `consistent`, 8 and 14 `uninformative`, 15 and 15 explained on both haplotypes.
- **Deletions equal to a junction.** 229 and 135 novel junctions of arm A equal a deletion one
  of the assembly haplotypes carries against GRCh38. On the assembly, 221 and 130 of them are
  `unspliced`: the reads run through without an intron, because the person's genome does not
  have those bases. The GRCh38 pipeline reports each as a novel intron. Where the deletion is
  heterozygous and at least 3 junction reads can be assigned, they come from the deletion
  haplotype only at 20 of 23 (HG03516) and 21 of 24 (HG02717) junctions; 3 and 2 are `mixed`
  and 1 (HG02717, in the IGH locus) is `contradicted`. This revises the reading of
  `unspliced` above: 44-48 % of that class is a deletion in the person, not alignment
  instability, and B1's `same-canonical` label for them says only that B1 has GRCh38's sequence
  there. The other 278 and 141 `unspliced` junctions have no equal deletion in the call set.
- **Verdicts.** 144 and 108 isoforms are `REFERENCE_BIAS` (18 and 28 with the FASTAs): 126 and
  80 by `haplotype_deletion`, 18 and 28 by `haplotype_variant`. SQANTI3's rules filter in arm A
  discards 142 and 104 of them (all of the `haplotype_variant` ones; 124 of 126 and 76 of 80 of
  the deletions).
- **Not seen.** `shifted-known` (12 and 10: the reads use an annotated junction at another
  position in the person) stays out of reach: 1 of them equals a deletion, the rest need the
  alignment on the personal genome.
- **Cost.** One run, genome-wide with `--bam` and the two VCFs: 3.6 and 5.9 min, 0.6 GB.

### With a call set that does not come from the assembly: 1000 Genomes short reads (0.1.1)

In the arm above the variant calls and the standard are the same assembly. Both individuals are
also in the 1000 Genomes high-coverage panel: SNVs, indels and SVs called from 30x Illumina
reads of 3,202 samples and phased with the pedigree (Byrska-Bishop et al., Cell 2022). That is
the input a study with short-read WGS of its sample has, and it owes nothing to the assembly.
`fetch_1kgp.sh` reads each person's genotypes through the remote index (in the regions of the
isoforms that have a novel junction; the panel is 30 GB) into one diploid VCF, and `C_1kgp`
runs it twice: as phased, and with every `|` written `/`, as an unphased call set gives them.
On chr22, 98.9 % and 98.3 % of each call set's SNV alleles are also in that person's assembly
calls (the scorer's same-person check). Standard as before: where the same reads splice on the
assembly (B2).

| | HG03516: assembly calls | 1000 Genomes, phased | phase removed | HG02717: assembly calls | 1000 Genomes, phased | phase removed |
|---|---:|---:|---:|---:|---:|---:|
| `V_rule` (junctions explained, not contradicted by the reads) | 265 | 103 | 97 | 179 | 90 | 87 |
| `V_rule and T` | 255 | 96 | 89 | 171 | 85 | 82 |
| precision (95 % CI) | 0.962 (0.94-0.98) | 0.932 (0.88-0.97) | 0.918 (0.87-0.97) | 0.955 (0.92-0.98) | 0.944 (0.89-0.99) | 0.943 (0.90-0.99) |
| sensitivity over `T` (95 % CI) | 47 % (42-51 %) | 18 % (14-21 %) | 16 % (13-19 %) | 53 % (48-59 %) | 26 % (22-31 %) | 25 % (21-30 %) |
| splice-site class found | 33 / 35 | 31 / 35 | 31 / 35 | 41 / 43 | 38 / 43 | 37 / 43 |
| `unspliced` found | 221 / 499 | 65 / 499 | 58 / 499 | 130 / 271 | 47 / 271 | 45 / 271 |
| novel junctions left `unlinked` | 0 | 0 | 49 | 0 | 0 | 25 |
| `REFERENCE_BIAS` isoforms (variant / deletion) | 144 (18 / 126) | 55 (18 / 37) | 48 (18 / 30) | 108 (28 / 80) | 65 (25 / 40) | 64 (25 / 39) |

- **What it calls is as often right.** 96 of 103 and 85 of 90 junctions are confirmed on the
  assembly (precision 0.93 and 0.94, against 0.96 and 0.96 with the assembly's own calls; the
  intervals overlap). The junctions outside `T` are the same kind as before:
  `same-noncanonical` in the person, or not assessable on the assembly.
- **The splice-site class is nearly all there.** 31 of 35 and 38 of 43 are found
  (33 and 41 with the assembly's calls).
- **Most deletions are not.** Of the junctions that equal a deletion in the assembly's calls,
  64 of 222 and 47 of 130 are found, and the share falls with length:

| length of the junction | HG03516: found | no deletion in the call set | a deletion with other coordinates, overlap ≥ 80 % | other | HG02717: found | no deletion in the call set | other coordinates, overlap ≥ 80 % | other |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| under 50 bp | 52 / 89 | 30 | 1 | 6 | 40 / 65 | 18 | 2 | 5 |
| 50-99 bp | 10 / 45 | 29 | 3 | 3 | 6 / 24 | 12 | 5 | 1 |
| 100-299 bp | 1 / 38 | 30 | 6 | 1 | 1 / 24 | 20 | 2 | 1 |
| 300 bp or more | 1 / 50 | 18 | 28 | 3 | 0 / 17 | 4 | 10 | 3 |

  Under 50 bp the tool finds 6 in 10. Above that the call set mostly has no deletion there
  (107 and 54 junctions over all lengths). Where it has one, usually a `<DEL>` from an SV
  caller, its coordinates differ from the junction's (38 and 19 overlap it by 80 % or more,
  some by 2 bp): the tool asks for a deletion equal to the intron, and it gives no credit to
  one that is nearly so. So overall sensitivity is 18 % and 26 % of `T`, against 47 % and 53 %.
- **At 13 and 9 of the missed junctions the motif is not evaluable.** 8 and 7 have a `<DEL>` of
  the call set at the junction that is not equal to it: a deletion that starts at the exon
  boundary and is longer than the 64 bp read on either side leaves no base to read on that
  haplotype. At the other 5 and 2 the call set writes the region as a cluster of small records.
- **Without phase little changes.** With every genotype unphased, 89 of 97 and 82 of 87 are
  confirmed. 49 and 25 novel junctions are `unlinked` (two heterozygous records near the
  junction whose relation the file does not give), which costs 15 and 10 of the assembly arm's
  junctions; 1 isoform is held for `haplotype_phase_unknown`. The read-haplotype check needs
  phased sites, so it reports `uninformative` throughout and the one contradicted junction is
  counted.
- **Verdicts.** 55 and 65 isoforms are `REFERENCE_BIAS` (48 and 64 unphased); SQANTI3's rules
  filter discards 55 and 62 of them.

So with the variant calls a WGS study has, the axis keeps its precision and its splice-site
class, and finds about a third of the deletions an assembly shows. Two things would raise that
and are not done: accepting a deletion that is nearly, not exactly, the intron (it would have to
be weighed against the exactness the verdict rests on), and reading a haplotype's motif across
a long deletion that is not the intron. The two runs of `C_1kgp` took 15 min (HG03516) and
20 min (HG02717) in all, 0.6 GB.

### Isoform level, firewall, and what the rescue changes for a user

SQANTI3's rules filter in arm A calls `Artifact` 9,886 of HG03516's adjudicated novel models
and `Isoform` 9,253 (HG02717: 7,202 / 10,067). PanIsoGuard's `REFERENCE_BIAS` isoforms (18
and 28) are all among the `Artifact`s in both individuals: without the rescue a user
following the filter discards them. The rest (HG03516 / HG02717, 0.1.1 binary, ruleset 0.0.8):
`NOVEL_PROBABLE` 7,068 / 8,364, every one with ≥ 3 whole-chain long reads (0.0.7 chain
support; no short reads exist for these samples, so nothing is `CONFIRMED`): 5,775 / 7,007
are SQANTI3 `Isoform`s and 1,293 / 1,357 are SQANTI3 `Artifact`s whose only sign is
`RTS_stage` (`chain_support;rt_switch_sign`; held `UNCONFIRMED` before ruleset 0.0.5), which
is the real-data size of the RT-switch decision; `KNOWN_FRAGMENT` 3,925 / 1,854 (ISM carried
end to end with no intra-priming sign; 314 / 177 of them SQANTI3 `Artifact`s); `UNCONFIRMED`
11,134 / 5,725 (ISM with the intra-priming sign 6,612 / 2,273; chain-supported isoforms with
an intra-priming sign 2,892 / 1,978 or a non-canonical motif 630 / 748, both still held; the
mapping sign without short reads 738 / 601, `ARTIFACT` before 0.0.5); `UNRESOLVED` 15,414 /
8,969; `ARTIFACT` 0 / 0, because no short reads were given.

Firewall: 18 and 28 rescued under `wgs` provenance, 0 under `unknown`, all held `UNRESOLVED`
with `circularity_flag`.

### Cost (32 cores, 123 GB; `/usr/bin/time -v`)

| arm | step | HG03516 wall | peak RSS | HG02717 wall | peak RSS |
|---|---|---:|---:|---:|---:|
| A | minimap2 splice:hq (existing) | 66 min | 23 GB | 2 h 8 min | 23 GB |
| A | IsoQuant (existing) | 47 min | n/a | n/a | n/a |
| A | SQANTI3 + rules filter, all models | 5.8 min | 4.9 GB | 4.2 min | 4.9 GB |
| B1 | build diploid SNV reference + GTF | 2.1 min | n/a | 1.1 min | n/a |
| B1 | minimap2 splice:hq, 6.2 Gb diploid | 2 h 13 min | 37 GB | 4 h 31 min | 37 GB |
| B1 | IsoQuant | 8 h 55 min | 15 GB | 2 h 22 min | 15 GB |
| B1 | SQANTI3 | 5.7 min | 9.6 GB | 8.0 min | 9.6 GB |
| B2 | asm5 PAF × 2 haplotypes | 1 h 41 min | 27 GB | 32 min | 29 GB |
| B2 | liftoff × 2 | 1 h 37 min | 78 GB | 1 h 8 min | 78 GB |
| B2 | build diploid + GTF | 2.9 min | n/a | 1.1 min | n/a |
| B2 | minimap2 splice:hq | 4 h 21 min | 38 GB | 3 h 32 min | 38 GB |
| B2 | IsoQuant | 22 min (resumed after an OOM kill next to B1) | 5.9 GB | 1 h 37 min | 17 GB |
| B2 | SQANTI3 | 8.9 min | 9.5 GB | 3.7 min | 9.4 GB |
| C | `adjudicate` × 3 provenance runs, genome-wide with `--bam` | 10.9 min | 0.6 GB | 18.0 min | 0.6 GB |
| C-VCF | `adjudicate`, one run with `--bam` and the two haploid VCFs | 3.6 min | 0.6 GB | 5.9 min | 0.6 GB |

HG03516's B1 and B2 ran partly concurrently, so its walls are upper bounds; HG02717 ran
sequentially (B1 7.0 h, B2 6.9 h end to end). The personalized pipelines cost 7-12 h and
38-78 GB on top of arm A; the post-hoc rescue costs 3-6 minutes per run (the three provenance
runs above; 1.8-1.9 min each before the 0.0.7 whole-chain read pass) and under 1 GB, and recovers
31 of 32 and 37 of 39 SNV-class junctions at the rule level (12 of 32 and 21 of 39 as
isoform verdicts, because of the all-junctions rule).

B2's SQANTI3 step failed silently in the first HG03516 run (liftoff writes CDS lines without
frames, which SQANTI3's `gtfToGenePred` rejects); the scorer never reads B2's SQANTI3
output, so no number depended on it, but the cost row did. `run.sh` now feeds SQANTI3 the
gene/transcript/exon lines only, and both individuals' rows above are from that fixed step.

## Deliverables

`run.sh` (all arms; `SAMPLE=`, `WORK=` and tool paths as env vars), `classify_junctions.py`
(read-level classes; `--same-coordinates` for B1, `--paf` liftover for B2), `score.py`
(metrics, bootstrap, cost table; `--vcf-arm armC_vcf/wgs` scores the VCF arm from the tool's
junction table and `--call-set NAME=PREFIX[,VCF]` a run with another call set;
`--label <ID> --emit-metrics` adds the individual to the tracked envelope checked by
`benchmark/results/collect.py --check`), `fetch_1kgp.sh` and `novel_isoform_regions.py` (the
1000 Genomes genotypes of the regions that hold a novel junction).
