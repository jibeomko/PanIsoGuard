# Reference bias three ways: GRCh38 pipeline vs personalized-genome pipeline vs PanIsoGuard

**Status: run on HG03516 and HG02717 (arms A, B1, B2, C; PanIsoGuard 0.0.6).** Numbers:
[Results](#results-hg03516-and-hg02717) below and
[../results/refbias_3way/metrics.json](../results/refbias_3way/metrics.json).

## Question

PanIsoGuard's reference-bias rescue is so far validated against its own criterion (a
splice-site SNV that is non-canonical on GRCh38 and canonical on the person's haplotype)
plus one non-circular phasing check ([refbias_phasing](../refbias_phasing): 27 of 28
heterozygous junctions). The only head-to-head table ([sqanti3_filter_h2h](../sqanti3_filter_h2h))
is on SQANTI-SIM, which is simulated from GRCh38 and therefore contains no reference bias
by construction, so the rescue cannot show there.

This benchmark gives the rescue an **independent standard**: the novel isoforms that a
full personalized-genome re-analysis of the same reads retains or discards. It asks

1. how much of what the personalized pipeline changes, the post-hoc rescue recovers
   (sensitivity, precision), and at what fraction of the compute;
2. which reference-bias mechanisms the rescue cannot see by construction (indels, SVs,
   mapping shifts, motif-*disrupting* variants), and how large each is.

Expected shape of the result: the rescue should recover most of the **splice-dinucleotide
SNV** class (that is what it models) and none of the rest. Arm B1 below isolates the first;
the difference B2 − B1 measures the second. A sensitivity well below 100 % against the full
personalized pipeline is the expected, reportable scope of the method, not a failure.

## Arms

All arms use the same reads (HG03516 Kinnex Iso-Seq flnc, 11.34 M aligned reads), the same
caller (IsoQuant 3.13, `pacbio_ccs`), the same annotation (GENCODE v49) and the same QC
(SQANTI3 6.0.1 + its rules filter). Only the genome the reads are aligned to differs.

| Arm | Genome | Annotation | Extra steps | What it measures |
|---|---|---|---|---|
| **A** | GRCh38 primary | GENCODE v49 | none | the reference-only pipeline (today's default) |
| **B1** | SNV-personalized diploid: `mat_consensus.fa` + `pat_consensus.fa` (GRCh38 backbone with the haplotype's assembly SNVs applied; GRCh38 coordinates) | GENCODE v49, contigs renamed per haplotype | none (same coordinates) | SNV-mediated reference bias only, the component PanIsoGuard models |
| **B2** | HPRC R2 phased assembly, diploid: `mat.fa.gz` + `pat.fa.gz` (72 + 74 contigs, chromosome-scale) | GENCODE v49 lifted per haplotype with liftoff | lift junctions back to GRCh38 (`paftools.js liftover` on the asm5 PAF) + read-level cross-check | the full personalized pipeline: SNVs, indels, SVs, assembly-specific sequence |
| **C** | GRCh38 (arm A's alignment and models) | GENCODE v49 | `panisoguard adjudicate --bam rna.aln.bam --reference GRCh38 --reference-haplotype mat_consensus.fa --reference-haplotype pat_consensus.fa --haplotype-provenance wgs` (no short reads exist for this sample, so the short-read axis is off) | PanIsoGuard post hoc on the reference-only pipeline |

Arm A alignment and IsoQuant come from [hg03516_refbias](../hg03516_refbias) (68,011
models, 37,568 novel); `A_sqanti` runs SQANTI3 and its rules filter on all of them.

Diploid rather than per-haplotype alignment for B1/B2 because that is what a personalized
pipeline does (each read picks its haplotype). Fallback if IsoQuant or SQANTI3 object to a
6.2 Gb diploid reference: run per haplotype and take the union (a junction counts as
canonical/known if it is on either haplotype, which is also the rescue's rule).

## Unit and labels

The unit is a **junction** (intron, GRCh38 coordinates, strand) used by at least one
caller-reported novel model (NIC/NNC) in arm A: `J_A` (98,672 for HG03516). Isoform-level
counts (what a user keeps after the SQANTI3 filter in each arm) are reported alongside.

For each `j ∈ J_A`, arm A gives its status on GRCh38 (`known`, `novel-canonical`,
`novel-noncanonical`) and the filter outcome of the models that use it.

For arm B (B1 and B2), the label comes from the **reads**, not from coordinate liftover
alone, because IsoQuant model coordinates can be displaced from where the reads splice
(the P2RX5 case in [hg03516_refbias](../hg03516_refbias)). Take the reads that span `j`
exactly in arm A; find where the same reads splice on the personal genome; lift that
position back to GRCh38 (B2 only; B1 is already in GRCh38 coordinates). Each read falls in
one class, and the junction takes the majority class:

| class | meaning |
|---|---|
| `same-canonical` | same junction, now canonical in the person: the non-canonical flag was reference bias |
| `same-noncanonical` | same junction, still non-canonical: not reference bias |
| `shifted-known` | the reads splice at a different, GENCODE-annotated junction in the person: the novelty was a mapping artifact of reference bias (indel, SV) |
| `shifted-novel` | different junction, still novel: ambiguous |
| `unspliced` | the reads align through the position without an intron in the person (B2: mostly short introns that minimap2 renders as a deletion on the assembly, see Results) |
| `no-reads` | no read spans the junction exactly in arm A (model coordinate displaced from the reads) |
| `unassessable` | read unaligned, or region not lifted |

Reference-bias truth: `T = same-canonical (non-canonical on GRCh38) ∪ shifted-known ∪
unspliced`, restricted to `j` that arm A calls novel; `T_core` leaves `unspliced` out (see
Results for why). `T_B1` (SNV only) and `T_B2` (full).

PanIsoGuard's sets, both under `wgs` provenance: `R_rule` = novel junctions that meet the
rescue's own criterion (non-canonical on GRCh38, canonical on a consensus haplotype), i.e.
what the variant axis sees; `R` = novel junctions inside isoforms that got the
`REFERENCE_BIAS` verdict (the verdict needs every novel junction of the isoform explained,
and only NIC/NNC isoforms are adjudicated). `R_unknown` (provenance `unknown`) must be
empty: the firewall check.

## Metrics

- sensitivity `|R ∩ T| / |T|` and precision `|R ∩ T| / |R|`, against `T_B1` and `T_B2`;
  bootstrap CIs over junctions.
- `T_B2 \ R` by mechanism: dinucleotide SNV missed (implementation), SNV outside the
  dinucleotide, indel at the site, SV / assembly-specific exon, alignment shift. This is the
  scope statement.
- `R \ T_B2`: every case inspected by hand (candidate false rescues).
- the blind spot in the other direction: arm-A `novel-canonical` junctions that become
  non-canonical or absent in the person (the scan's DISRUPTED class). PanIsoGuard does not
  model these; report their count.
- isoform level: novel models kept by the SQANTI3 rules filter in A, in B1, in B2, and kept
  by PanIsoGuard in C (`HIGH/NOVEL_PROBABLE` + `REFERENCE_BIAS`); in
  particular "discarded in A, kept in B2" vs "rescued in C".
- cost per arm: wall time, CPU-hours, peak RSS, disk; C counted as adjudicate only (arm A
  is sunk cost for both A and C).

## Inputs on disk (HG03516, `/mnt/Data/pig_hg03516`, 38 GB)

| file | what | made by |
|---|---|---|
| `flnc.bam` (13.5 GB) | Kinnex Iso-Seq flnc, HPRC R2 | `dl.sh` |
| `rna.aln.bam` | minimap2 `splice:hq` vs GRCh38; 11,342,107 primary alignments; 66 min wall / 17.4 CPU-h / 23 GB at 16 threads | `rna.sh` |
| `isoquant_out/iq/iq.transcript_models.gtf` | IsoQuant 3.13; 47 min at 32 threads | `isoquant.sh` |
| `mat.fa.gz`, `pat.fa.gz` | HPRC R2 v1.1.0 phased assembly (`HG03516#2#…`, `#1#…`) | `dl.sh` |
| `{mat,pat}.vcf`, `.snv.vcf.gz` | paftools.js call vs GRCh38 (3.2 M SNVs, 0.4 M 1-bp indels, 0.6 M longer indels, 18 k SVs ≥ 50 bp per haplotype) | `geno.sh` |
| `{mat,pat}_consensus.fa` | bcftools consensus, SNV only, length-matched to GRCh38 | `geno.sh` |
| `scan/`, `refbias_unique.tsv` | the 40-junction rescue subset and its per-chromosome adjudicate inputs | `scan_refbias.sh` |

HG02717 (`/mnt/Data/pig_hg02717`) has the same set (69,228 models, 24,941 novel; 46
rescue junctions). Run HG03516 first, replicate on HG02717 once the scripts are fixed.

Reference: `/mnt/Data/LR-RNAseq/reference/GRCh38.primary_assembly.genome.fa`,
`gencode.v49.annotation.gtf`.

## Prerequisites and known snags

- **liftoff**: `conda create -n liftoff --override-channels -c conda-forge -c bioconda liftoff=1.6.3`
  (bioconda has 1.6.3; it pulls minimap2); `run.sh` expects it at `$CONDA/liftoff/bin/liftoff`.
- **SQANTI3** (`~/miniconda3/envs/sqanti3`) works only with that env's `bin` first in
  `PATH`: the wrapper's `#!/usr/bin/env python` otherwise picks another python and fails on
  `BCBio`. `export PATH=$HOME/miniconda3/envs/sqanti3/bin:$PATH` as the h2h script does.
- **The asm5 PAF** (`B2_paf`) is kept under `paf/`; `paftools.js liftover <aln.paf>
  <query.bed>` lifts BED in **query (assembly)** coordinates to the target (GRCh38) and
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
| same-canonical, canonical on GRCh38 (no reference bias) | 17,673 | 16,974 | 13,554 | 13,117 |
| same-noncanonical, non-canonical on GRCh38 (not reference bias) | 2,727 | 2,503 | 2,058 | 1,877 |
| **same-canonical, non-canonical on GRCh38 (splice-dinucleotide SNV; the rescue's target)** | **32** | **35** | **39** | **43** |
| **shifted-known** (reads use an annotated junction in the person) | **7** | **12** | **11** | **10** |
| same-noncanonical, canonical on GRCh38 (motif *disrupted* in the person) | 23 | 30 | 7 | 6 |
| shifted-novel | 169 | 255 | 83 | 171 |
| unspliced | 49 | 499 | 21 | 271 |
| unassessable / no-reads | 399 / 20 | 771 / 20 | 223 / 18 | 501 / 18 |

`T_core` (SNV class + shifted-known) = 39 / 47 (HG03516 B1 / B2) and 50 / 53 (HG02717). B2
adds the SNV-class and shifted-known junctions that the SNV-only genome cannot show: indels
and assembly-specific sequence at the splice site.

**`unspliced` is left out of `T_core`.** B2's 499 (HG03516) are short introns (median 82 bp,
55 % ≤ 100 bp; median 4 reads) at which minimap2 places a deletion instead of an intron on
the diploid assembly; 308 of them are `same-canonical` on the SNV genome (B1, same
coordinates), so the reads do splice there. That is alignment instability of short introns
on a diploid assembly, not reference bias, and no post-hoc method could or should "rescue"
it. HG02717 shows the same (271 on B2, 21 on B1).

### Rule level and verdict level

| | HG03516 B1 | HG03516 B2 | HG02717 B1 | HG02717 B2 |
|---|---:|---:|---:|---:|
| `R_rule` (junctions meeting the rescue criterion) | 35 | 35 | 42 | 42 |
| `R_rule ∩ T` | 31 | 31 | 37 | 38 |
| rule precision | 0.886 | 0.886 | 0.881 | 0.905 |
| rule sensitivity, SNV class | **31 / 32 = 0.969** | 31 / 35 = 0.886 | **37 / 39 = 0.949** | 38 / 43 = 0.884 |
| rule sensitivity, `T_core` | 0.795 | 0.660 | 0.740 | 0.717 |
| `R` (junctions in `REFERENCE_BIAS` isoforms) | 13 (18 isoforms) | 13 | 22 (28 isoforms) | 22 |
| `R ∩ T` | 12 | 12 | 21 | 19 |
| verdict precision (95 % CI) | 0.923 (0.77–1.00) | 0.923 | 0.955 (0.86–1.00) | 0.864 (0.73–1.00) |
| verdict sensitivity, SNV class (95 % CI) | 0.375 (0.22–0.56) | 0.343 (0.20–0.51) | 0.538 (0.38–0.69) | 0.442 (0.30–0.58) |

The motif criterion itself finds 31 of 32 (HG03516) and 37 of 39 (HG02717) SNV-class
junctions on B1. The `R_rule` junctions outside `T` (4 and 5) are `same-noncanonical` in the
person: the haplotype with the canonical motif exists, but the reads come from the other
one (the rescue is allele-agnostic by design; [refbias_phasing](../refbias_phasing) found
the same 1-in-10). Shifted-known junctions (indel / SV mapping shifts) are never caught: the
method does not model them.

The drop from `R_rule` to `R` (35 → 13 and 42 → 22) is not the criterion. Of the SNV-class
junctions in `T_B1` without a rescue verdict (20 and 18): 15 and 10 sit in isoforms with
another novel junction the haplotypes do not explain (the verdict needs all of them), 4 and
6 are in isoforms outside NIC/NNC (not adjudicated), 1 and 2 are criterion misses. B2 adds
3 criterion misses in each individual, the non-SNV variants. `R ∖ T` is one junction in
HG03516 (chr6:136,665,560–136,668,208) and one locus in HG02717 (three junctions sharing the
acceptor chr6:130,022,220), all `same-noncanonical` in the person: candidate false rescues,
1 of 13 and 1–3 of 22.

### Isoform level, firewall, and what the rescue changes for a user

SQANTI3's rules filter in arm A calls `Artifact` 9,886 of HG03516's adjudicated novel models
and `Isoform` 9,253. PanIsoGuard's `REFERENCE_BIAS` isoforms (18 and 28) are all among the
`Artifact`s in both individuals: without the rescue a user following the filter discards
them. The rest: `ARTIFACT` 741 / 605 and `UNRESOLVED` 18,380 / 20,055 (no short reads exist
for these samples, so nothing is confirmed).

Firewall: 18 and 28 rescued under `wgs` provenance, 0 under `unknown`, all held `UNRESOLVED`
with `circularity_flag`.

### Cost (32 cores, 123 GB; `/usr/bin/time -v`)

| arm | step | HG03516 wall | peak RSS | HG02717 wall | peak RSS |
|---|---|---:|---:|---:|---:|
| A | minimap2 splice:hq (existing) | 66 min | 23 GB | 2 h 8 min | 23 GB |
| A | IsoQuant (existing) | 47 min | — | — | — |
| A | SQANTI3 + rules filter, all models | 5.8 min | 4.9 GB | 4.2 min | 4.9 GB |
| B1 | build diploid SNV reference + GTF | 2.1 min | — | 1.1 min | — |
| B1 | minimap2 splice:hq, 6.2 Gb diploid | 2 h 13 min | 37 GB | 4 h 31 min | 37 GB |
| B1 | IsoQuant | 8 h 55 min | 15 GB | 2 h 22 min | 15 GB |
| B1 | SQANTI3 | 5.7 min | 9.6 GB | 8.0 min | 9.6 GB |
| B2 | asm5 PAF × 2 haplotypes | 1 h 41 min | 27 GB | 32 min | 29 GB |
| B2 | liftoff × 2 | 1 h 37 min | 78 GB | 1 h 8 min | 78 GB |
| B2 | build diploid + GTF | 2.9 min | — | 1.1 min | — |
| B2 | minimap2 splice:hq | 4 h 21 min | 38 GB | 3 h 32 min | 38 GB |
| B2 | IsoQuant | 22 min (resumed after an OOM kill next to B1) | 5.9 GB | 1 h 37 min | 17 GB |
| B2 | SQANTI3 | 8.9 min | 9.5 GB | 3.7 min | 9.4 GB |
| **C** | **`adjudicate` × 3 provenance runs, genome-wide with `--bam`** | **5.4 min** | **0.4 GB** | **5.7 min** | **0.5 GB** |

HG03516's B1 and B2 ran partly concurrently, so its walls are upper bounds; HG02717 ran
sequentially (B1 7.0 h, B2 6.9 h end to end). The personalized pipelines cost 7–12 h and
38–78 GB on top of arm A; the post-hoc rescue costs 5–6 minutes and under 1 GB, and recovers
31 of 32 and 37 of 39 SNV-class junctions at the rule level (12 of 32 and 21 of 39 as
isoform verdicts, because of the all-junctions rule).

B2's SQANTI3 step failed silently in the first HG03516 run (liftoff writes CDS lines without
frames, which SQANTI3's `gtfToGenePred` rejects); the scorer never reads B2's SQANTI3
output, so no number depended on it, but the cost row did. `run.sh` now feeds SQANTI3 the
gene/transcript/exon lines only, and both individuals' rows above are from that fixed step.

## Deliverables

`run.sh` (all arms; `SAMPLE=`, `WORK=` and tool paths as env vars), `classify_junctions.py`
(read-level classes; `--same-coordinates` for B1, `--paf` liftover for B2), `score.py`
(metrics, bootstrap, cost table; `--label <ID> --emit-metrics` adds the individual to the
tracked envelope checked by `benchmark/results/collect.py --check`).
