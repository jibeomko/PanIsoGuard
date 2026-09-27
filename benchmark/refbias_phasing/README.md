# Phasing check of the reference-bias rescue (non-circular)

**Question:** when PanIsoGuard rescues a novel junction as reference bias, is the junction really
used by the haplotype on which it is canonical? The earlier benchmarks cannot say: their truth
label (`CREATED`: non-canonical on GRCh38, canonical on the individual's haplotype) is the rescue
rule itself ([docs/validation.md](../../docs/validation.md)). Reading the motif-creating allele in
the junction reads is impossible too, because that base is intronic and spliced out.

**Test.** For a *heterozygous* reference-bias junction (canonical on one haplotype H only), the
long reads that use the junction should carry H's alleles at the gene's other heterozygous SNVs,
while the other haplotype O is still expressed at the locus (reads not using the junction carry
O alleles). Phase comes from the trio-phased HiFi assembly (pat/mat), not from the RNA.

**Answer: yes, in 27 of 28 testable junctions.**

## Result

Same work dirs as [hg03516_refbias](../hg03516_refbias) / [refbias_cohort](../refbias_cohort)
(PacBio Kinnex Iso-Seq, minimap2 to GRCh38; pat/mat SNVs from `paftools call`). The junction lists
include the one junction per individual that the cohort analysis dropped for zero exact-coordinate
reads (chr17:3681978-3685514, homozygous in both, so untestable either way).

| | HG03516 | HG02717 | total |
|---|---:|---:|---:|
| reference-bias junctions | 41 | 47 | 88 |
| homozygous (cannot be phased) | 18 | 16 | 34 |
| heterozygous | 23 | 31 | 54 |
| &nbsp;&nbsp;too few phased junction reads (< 3) | 9 | 14 | 23 |
| &nbsp;&nbsp;other haplotype not seen at the locus | 2 | 1 | 3 |
| &nbsp;&nbsp;**tested** | **12** | **16** | **28** |
| &nbsp;&nbsp;&nbsp;&nbsp;consistent (all phased junction reads from H) | 11 | 16 | **27** |
| &nbsp;&nbsp;&nbsp;&nbsp;contradicted (junction reads from O) | 1 | 0 | **1** |
| junction reads, tested junctions: H / O | 172 / 3 | 400 / 0 | 572 / 3 |
| other reads at the same loci: H / O | 669 / 654 | 1839 / 1078 | |

Both haplotypes are expressed at the tested loci (last row), yet 572 of 575 phased junction reads
come from the haplotype that carries the canonical motif. Per-junction counts and one-sided Fisher
p-values are in [../results/refbias_phasing/metrics.json](../results/refbias_phasing/metrics.json).

## The one contradiction

HG03516 chr3:184709997-184735943 (+). GRCh38 reads GT…AC; the maternal haplotype AT…AC (canonical,
U12 type), the paternal AC…AC. All 3 phased junction reads are paternal, and all 3 carry a 2-bp
insertion 3 bp upstream of the donor (`57M2I3M25946N…`), so the junction coordinate looks like an
alignment artifact rather than maternal splicing. PanIsoGuard's BAM axis would flag these reads
(indel-near fraction 1.0 > 0.5), but the reference-bias rescue is checked before the mapping
mechanism, so it is rescued anyway. Letting a mapping flag block the rescue is *not* a clean fix:
HG02717 chr2:88861556-88947305 also has indel-near 1.0 on all 81 junction reads and passes the
phasing test (11/11 phased reads from H).

Three heterozygous junctions below the 3-read threshold lean the same way: HG03516
chr6:136665560 (1 H, 1 O), chr7:45791846 (0 H, 2 O) and HG02717 chr4:56466745 (0 H, 2 O).

## Caveats

- Only heterozygous junctions can be tested; the 34 homozygous ones remain unchecked.
- A SNV missing from one haplotype's VCF is read as the reference allele there; about 9% of GRCh38
  is not covered by exactly one contig per haplotype, so a few "heterozygous" sites may be wrong.
  A read is assigned only if its alleles favour one haplotype by at least 2 sites, which keeps a
  single bad site or sequencing error from flipping it.
- "Consistent" shows the junction is used from the haplotype that makes it canonical. It does not
  show whether the resulting isoform is functional, and it does not test the pangenome axis.
- The thresholds >= 3 phased junction reads and >= 2 other-haplotype locus reads were fixed before
  the first run. The read-assignment rule was changed once after it: the first version kept only
  sites where the RNA showed both alleles, and at HG03516 chr7:144278475 that left a single site
  with an error-prone neighbour, which made the junction look contradicted (2 H, 2 O). Using all
  assembly-het sites with the margin-2 rule gives 4 H, 0 O there. The verdict order was also
  changed then, so that junction reads from O count as a contradiction even when O is otherwise
  absent at the locus (that is what makes chr3:184709997 contradicted rather than uninformative).

## Reproduce

```bash
# work dirs from benchmark/hg03516_refbias/run.sh (one per individual); needs pysam + scipy
python3 phase_check.py HG03516:/path/to/hg03516 HG02717:/path/to/hg02717 \
  --emit-metrics ../results/refbias_phasing/metrics.json
```

Runs in a few seconds from the aligned BAM and the two SNV VCFs.
