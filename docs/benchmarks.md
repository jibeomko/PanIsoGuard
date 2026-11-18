# Benchmarks and limitations

[Back to the README](../README.md#benchmark--limitations).

## When to use PanIsoGuard

Use PanIsoGuard to inspect the evidence and reasoning behind novel-isoform verdicts,
including possible reference bias and candidates held for review. The comparisons below
show the precision-recall trade-offs under each tested input configuration.

- **You ran more than one caller** (FLAIR, IsoQuant, Bambu, ESPRESSO, TALON and others) and they
  disagree about the novel isoforms. PanIsoGuard matches isoforms across callers by their
  intron chains and records how many callers found each one. On the five-caller chr22 simulation, novel chains unique to one
  caller have much lower precision than consensus chains under the tested matching rules ([benchmark/multicaller](../benchmark/multicaller)).
- **You want a clear class and reason for each novel isoform**, not just a filtered GTF. Each
  isoform gets one of 8 classes and a `rule_trace` that scripts can read. You can keep `NOVEL_CONFIRMED` and `NOVEL_PROBABLE` and
  review the rest, instead of looking at reads one by one.
- **You have the sample's haplotypes or a pangenome** and want to investigate possible
  *reference bias*: splice motifs made canonical by sample variants, deletions in the sample
  that a reference pipeline reports as introns, or junctions supported by population paths.
  Population support does not establish the sample's genotype or expression. A provenance gate
  prevents RNA-derived or unknown-source genomic inputs from providing an independent rescue
  explanation. Against a full personalized-genome re-analysis of the same reads (HG03516 and
  HG02717, [benchmark/refbias_3way](../benchmark/refbias_3way)), with the assembly's variant
  calls given as VCFs, 255 of 265 and 171 of 179 junctions PanIsoGuard attributes to the
  genome are confirmed on the person's assembly (precision 0.96). That is 47 % and 53 % of the
  junctions the assembly accounts for: 33 of 35 and 41 of 43 made canonical by a variant at the
  splice site, and 221 and 130 that equal a deletion in the person. One run takes 4-6 minutes
  instead of the 7-11 hours of the personalized pipeline. Consensus FASTAs carry substitutions
  only and reach the splice-site SNV class alone (31 of 32 and 37 of 39 on an SNV-personalized
  genome). As isoform verdicts fewer are recovered, because a verdict needs every novel
  junction of the isoform explained (144 and 108 isoforms). The reads that use a junction
  explained on one haplotype come from that haplotype at 27 of 28 testable junctions
  ([benchmark/refbias_phasing](../benchmark/refbias_phasing)); `adjudicate --bam` makes the
  same check per junction and does not count a contradicted one. Homozygous junctions cannot
  be checked this way. A junction that an indel or structural variant moves to another
  position it does not see.

PanIsoGuard evaluates existing caller output. It can consume SQANTI3's QC results or compute
structural categories, splice motifs and downstream A content with `--no-sqanti3`; that mode
does not evaluate SQANTI3's RT-switch flag. See the
[SQANTI3 comparison](relationship_to_sqanti3.md) and the
[comparison with merge tools](relationship_to_merge_tools.md) for the scope of each tool.

## How it compares

Test data: simulated reads with known answers (SQANTI-SIM on GENCODE v49 chr22, FLAIR +
SQANTI3), 870 isoforms that do not fully match a known transcript (730 real, 140 false). All
methods got the same inputs. "1M" means 1 million simulated Illumina read pairs aligned with
STAR. AUPRC measures how well a method ranks real isoforms above false ones (1 is perfect).
Full tables, bootstrap confidence intervals and caveats are in
[benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h).

| short reads | method | precision | recall | F1 | AUPRC |
|---|---|---:|---:|---:|---:|
| none | SQANTI3 rules filter | 0.918 | 0.886 | 0.902 | 0.909 |
| none | PanIsoGuard | n/a (0 calls) | 0.000 | n/a | 0.824 |
| none | PanIsoGuard + BAM | 0.993 | 0.748 | 0.853 | 0.953 |
| 1M | SQANTI3 rules filter | 0.918 | 0.890 | 0.904 | 0.909 |
| 1M | Every junction ≥ 3 short reads | 0.976 | 0.945 | 0.960 | 0.968 |
| 1M | PanIsoGuard | 0.990 | 0.545 | 0.703 | 0.970 |
| 1M | PanIsoGuard + BAM | 0.993 | 0.752 | 0.856 | 0.975 |

- In this comparison, PanIsoGuard is the most precise, but it keeps fewer of the real isoforms (lower recall).
  With short reads alone it holds back on NIC isoforms that only combine known junctions
  (there is no new junction to confirm) and on ISM isoforms. With the BAM (0.0.7), NIC
  isoforms whose whole intron chain is carried by at least 3 long reads are kept as
  `NOVEL_PROBABLE`; that is the recall difference between the last two rows (151 real
  isoforms gained, 0 false).
- For ranking (AUPRC) it ties the one-line short-read rule, with or without the BAM
  (every difference is within the bootstrap interval).
- Without short reads it confirms nothing unless you give it the BAM or caller agreement.
  With the BAM it keeps 546 of the 730 real isoforms at precision 0.993 and ranks above
  SQANTI3's rules filter (0.953 vs 0.909); the long reads are the caller's own, so these
  calls are `NOVEL_PROBABLE`, never `NOVEL_CONFIRMED`.
- The same thresholds on two chromosomes they never saw (chr20, chr19) give the same picture:
  without short reads AUPRC 0.962 and 0.957 against 0.946 and 0.928 for SQANTI3's rules filter,
  both bootstrap intervals above zero; with 1 M pairs a tie with the one-line rule
  ([benchmark/sqanti_sim_heldout](../benchmark/sqanti_sim_heldout)).
- The same thresholds on the same chr22 design re-simulated as ONT cDNA (NanoSim) keep the
  ordering: without short reads precision 0.76 against 0.18 for SQANTI3's rules filter (AUPRC
  0.549 vs 0.172), with 1 M pairs 0.88 and a tie with the one-line rule on AUPRC; the `hifi` and
  `ont` presets give the same class to every isoform, because the indel-near sign is read in
  excess of the reads' own indel rate ([benchmark/sqanti_sim_ont](../benchmark/sqanti_sim_ont)).
- On the other callers' output (IsoQuant, Bambu, ESPRESSO, TALON, same alignment, no short
  reads) it raises precision every time, by 0.01-0.02 where the caller is already strict and by
  0.32 on TALON, where SQANTI3's rules filter adds nothing; it is not a substitute for a caller's
  own read-count filter but stacks on it ([benchmark/caller_filters](../benchmark/caller_filters)).
- It decides fewer isoforms and is wrong less often on those it decides: at 1 M pairs with
  the BAM it gives a decisive class to 74 % of the isoforms and 1.4 % of those are wrong;
  the one-line rule decides all of them and 6.6 % are wrong
  ([benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h), risk-coverage).
- On real data with experimentally validated isoforms (LRGASP WTC11, PacBio and ONT), none
  of the validated isoforms the caller reconstructed is called an artifact, and without
  short reads the whole-chain reads keep all 17 (PacBio) and 16 (ONT) as
  `NOVEL_PROBABLE`, where SQANTI3's default rules filter discards 2 and 4 of them for its
  RT-switch flag; with the other platform's reads as an independent library
  (`--validation-bam`) 13 of 17 and 14 of 16 are still carried by at least 3 reads
  ([benchmark/lrgasp_wtc11](../benchmark/lrgasp_wtc11)).
- Real-data negative controls are limited: failure of amplicon support does not prove a
  transcript is absent. The 26 amplicon-unsupported LRGASP chains, adjudicated
  as would-be isoforms with no caller in the loop, are confirmed 0 times on PacBio and 3 times on
  ONT, where 8-60 reads carry them end to end; SQANTI3's rules filter keeps 18 of the 26
  ([benchmark/lrgasp_targets](../benchmark/lrgasp_targets)). Of the 13 novel isoforms IsoQuant
  reports in TUSCO's 46 single-isoform genes, PanIsoGuard confirms 3 and SQANTI3's filter keeps 8
  ([benchmark/tusco_negatives](../benchmark/tusco_negatives)).

These comparisons support a conservative operating point with explicit abstention.
With short reads, AUPRC is similar to the simple junction-support rule within bootstrap
uncertainty; the benefits include per-call reasoning and genomic reference-bias assessment.

## Validation

[docs/validation.md](validation.md) lists what has been checked and what is planned
(SQANTI-SIM, HG002/HPRC, LRGASP).

**Thresholds.** On SQANTI-SIM (GENCODE v49 chr22), trying a grid of threshold values changes
AUPRC very little (0.970-0.971), and the default is within 0.0001 of the best
([benchmark/results/sqanti_sim/sweep.tsv](../benchmark/results/sqanti_sim/sweep.tsv)); the
whole-chain read threshold is flat from 2 to 10 ([benchmark/chain_support](../benchmark/chain_support)).
The same thresholds, untouched, hold on chr20 and chr19
([benchmark/sqanti_sim_heldout](../benchmark/sqanti_sim_heldout)) and on the chr22 design re-simulated
as ONT reads ([benchmark/sqanti_sim_ont](../benchmark/sqanti_sim_ont)).

**Pangenome.** The pangenome reference-bias check was tested on the real HPRC v1.1 chr22 graph
(GATE-1): no false rescues on real FLAIR novel junctions, correct rescues on real population
deletions, and the circularity guard held when the junction file was marked as not independent
([benchmark/pangenome](../benchmark/pangenome)). It has been tested on chr22 only; reading the
graph (GBZ) directly is future work.

The confidence classes are labels, not an ordered scale or probabilities. Their observed
precision depends on which inputs are supplied; see the
[per-class results](usage.md#confidence-classes). Motif-defined rescue counts check the
implementation of the rescue rule, not independent biological accuracy. The
[phasing analysis](../benchmark/refbias_phasing) tests 28 heterozygous junctions from two
individuals, with 27 consistent and one contradicted; homozygous and insufficiently covered
junctions cannot be checked this way. See [validation](validation.md) for the distinction
between implementation checks and independent evidence.

## Runtime & memory

One thread, one `adjudicate` run with the long-read BAM, measured with `/usr/bin/time`
(GENCODE v49 unless noted):

| Data set | Isoforms | BAM | Other inputs | Wall time | Peak RAM |
|---|---:|---:|---|---:|---:|
| SQANTI-SIM chr22, FLAIR | 1,579 | 14 MB | SJ.tab | 0.5 s | 18 MB |
| GM12878 ONT direct RNA, 285 novel chains | 285 | 0.27 GB | | 4 s | 0.41 GB |
| LRGASP WTC11 PacBio, IsoQuant | 78,256 | 3.5 GB | SJ.tab | 2.2 min | 0.59 GB |
| LRGASP WTC11 ONT, IsoQuant | 39,542 | 12.4 GB | SJ.tab | 5.9 min | 0.44 GB |
| HG03516 Kinnex, IsoQuant | 68,011 | 5.1 GB | Two haplotype FASTAs | 3.3 min | 0.47 GB |
| HG03516 Kinnex, IsoQuant | 68,011 | 5.1 GB | Two haploid VCFs (8.0 M records) | 3.6 min | 0.6 GB |
| `decoys` on the same BAMs | 139,003 decoys | | | 1-12 min | 0.41 GB |

Building the GENCODE v49 catalog takes about 3 s of that; the rest is reading the BAM at each
junction and chain. The callers and the alignment take most of the total time of a study.
The measured adjudication runs in the table used under 0.6 GB of peak RAM; runtime and memory
on other inputs depend on their size and coverage.
