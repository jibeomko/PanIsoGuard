# Reference-bias rescue — cohort yield across four GIAB individuals

Extends the single-sample [hg002](../hg002) whole-genome rescue to a **four-individual
GIAB cohort** — HG001 (NA12878), HG002 (NA24385), HG003 (NA24149), HG004 (NA24143) — to
show the reference-bias rescue is a **consistent, perfectly-specific guardrail across
personalized human genomes**, not a one-off. Annotation + variant scan only (no RNA-seq,
no caller runs); all GIAB v4.2.1 / GENCODE v49 data is public.

For each individual: per chromosome, build the SNV-consensus haplotype, scan every GENCODE
intron for a **reference-bias** splice site (`CREATED`: non-canonical on the linear
reference, canonical on that individual's haplotype), and adjudicate the variant axis
(`--haplotype-provenance wgs` for the rescue, `unknown` for the firewall).

## Result (whole-genome, chr1–22)

| individual | reference-bias junctions (`CREATED`) | rescued (`wgs`) | false rescues | firewall-held (circular) |
|------------|:---:|:---:|:---:|:---:|
| HG001 | 37 | 37 | 0 | 37 |
| HG002 | 33 | 33 | 0 | 33 |
| HG003 | 31 | 31 | 0 | 31 |
| HG004 | 36 | 36 | 0 | 36 |
| **cohort** | **137** | **137 / 137 (100 %)** | **0** | **137 / 137** |

**Reading.** Each personalized genome carries **31–37** genuine reference-bias splice
junctions — junctions a *reference-only* pipeline would wrongly report as novel. PanIsoGuard's
variant axis **exonerates every one (137/137, 100 % sensitivity)** with **0 false rescues**,
and the circularity firewall **holds all 137** when the provenance is circular-risk. So the
rescue is a **consistent, perfectly-specific safeguard across individuals** — a few dozen
false-novel calls correctly reclassified per personalized genome — whose practical value
rises with how far a sample diverges from the linear reference (non-reference individuals,
personalized-genome / pangenome contexts).

This is the **impact case** for the rescue (PanIsoGuard's genuine differentiator): rare per
intron (~1 in 18k) but recurrent per genome and caught with perfect specificity. Recorded
in [../results/giab_cohort_rescue](../results/giab_cohort_rescue).

## Reproduce

```bash
# download the four GIAB v4.2.1 GRCh38 benchmark VCFs, pre-split GENCODE by chromosome,
# fill in the paths at the top of run_cohort.sh, then:
./run_cohort.sh
```
