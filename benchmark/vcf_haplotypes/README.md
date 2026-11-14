# The variant axis from a VCF equals the variant axis from consensus FASTAs

`adjudicate --haplotype-vcf` (0.1.0) reads a sample's haplotypes from an indexed VCF or BCF
instead of consensus FASTAs built with `bcftools consensus`. The two inputs must give the same
verdicts and the same per-junction motifs, since the VCF route is defined as the FASTA route's
semantics (the genotype's first allele is haplotype 1, the second haplotype 2, as `-H 1` / `-H 2`
build them) applied at the splice dinucleotides only. This protocol checks that on every variant-axis
run already on disk.

| data set | runs | isoforms | same class and reason | `REFERENCE_BIAS` (FASTA / VCF) | junctions | same variant call | same haplotype-1 motif |
|---|---:|---:|---:|---:|---:|---:|---:|
| HG002 chr22 (GIAB v4.2.1 biallelic SNVs) | 1 | 8 | 8 | 1 / 1 | 8 | 8 | 8 |
| HG03516, every chromosome, pat + mat | 46 | 412 | 412 | 59 / 59 | 412 | 412 | 412 |
| HG02717, every chromosome, pat + mat | 46 | 347 | 347 | 63 / 63 | 347 | 347 | 347 |

767 of 767 isoforms and 767 of 767 junctions agree. The VCF route marks 302 junctions `phased`
(a variant touched the motif, from a phased or homozygous genotype) and none `unphased`; no indel or
symbolic allele sat at a motif in these files (SNV-only VCFs), so the two counters the provenance
log carries for them are 0 here.

Why the VCF route exists: a consensus FASTA that also applies indels is shifted relative to the
reference, and the splice coordinates of a GRCh38 alignment then read the wrong bases on it. On
HG002 chr22 a consensus that applied the multiallelic indels too was 50 bp longer and read wrong
bases at 7 junctions; the VCF route cannot do that, because it applies same-length substitutions
only and skips (and counts) everything else. The `REFERENCE_BIAS` counts here sum over the two
parental haplotypes of each scan and so count a homozygous junction twice; the per-individual
reference-bias numbers are in [refbias_cohort](../refbias_cohort).

## Reproduce

```bash
bash run.sh      # re-runs the VCF route next to the existing FASTA-route outputs, then compare_routes.py
```

`run.sh` uses the outputs of [hg002](../hg002) (`vt.*`, `chr22.snps.vcf.gz`) and of
[hg03516_refbias](../hg03516_refbias) for both individuals (`scan/<chrom>.<hap>/adj_wgs.*`,
`<hap>.snv.vcf.gz`) and writes `vtvcf.*` / `adjvcf_wgs.*` beside them. Envelope:
`benchmark/results/vcf_haplotypes/metrics.json`.
