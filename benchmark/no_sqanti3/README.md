# `--no-sqanti3` against SQANTI3

`adjudicate --no-sqanti3` computes the classification rows PanIsoGuard otherwise reads from
SQANTI3: the structural category from the intron chain against the reference GTF, the splice
motif and the 3′ A content from the genome. This protocol adjudicates the five chr22
SQANTI-SIM callers twice, from SQANTI3 6.0.1's classification and from the computed one, in
the two settings of the main text, and compares the categories, the motif prior, the
confidence classes, and both modes' accuracy against the truth.

## Result, no short reads + BAM

| caller | multi-exon isoforms | categories equal | `all_canonical` equal | classes equal | SQANTI3 mode P / R | computed mode P / R |
|---|---:|---:|---:|---:|---:|---:|
| FLAIR | 1530 | 98.3 % | 1530 / 1530 | 98.4 % | 0.993 / 0.748 | 0.993 / 0.778 |
| IsoQuant | 1441 | 98.3 % | 1441 / 1441 | 98.3 % | 0.991 / 0.708 | 0.991 / 0.75 |
| Bambu | 552 | 99.6 % | 552 / 552 | 99.6 % | 0.764 / 0.951 | 0.761 / 0.951 |
| ESPRESSO | 1500 | 100.0 % | 1500 / 1500 | 100.0 % | 0.987 / 0.939 | 0.987 / 0.939 |
| TALON | 2703 | 98.0 % | 2703 / 2703 | 98.1 % | 0.716 / 0.709 | 0.702 / 0.749 |

Over the 7726 multi-exon isoforms the categories agree on 98.6 %. Mono-exonic
isoforms are classified too (inside an exon of a mono-exonic transcript, of a multi-exonic
one, overlapping an exon, inside a transcript span) but are not adjudicated in either mode.
The splice-motif prior agrees everywhere both modes evaluate it; the RT-switch flag is not
computed, so the `rt_switch` sign is absent in the computed mode.

With 1 M short-read pairs + BAM:

| caller | classes equal | SQANTI3 mode P / R | computed mode P / R |
|---|---:|---:|---:|
| FLAIR | 96.1 % | 0.993 / 0.752 | 0.993 / 0.784 |
| IsoQuant | 96.5 % | 0.991 / 0.717 | 0.991 / 0.758 |
| Bambu | 94.9 % | 0.75 / 0.961 | 0.747 / 0.961 |
| ESPRESSO | 98.0 % | 0.987 / 0.944 | 0.987 / 0.944 |
| TALON | 96.7 % | 0.712 / 0.704 | 0.698 / 0.746 |

## Where the categories differ

| SQANTI3 | computed | isoforms |
|---|---|---:|
| `incomplete-splice_match` | `novel_in_catalog` | 100 |
| `genic` | `novel_not_in_catalog` | 5 |
| `antisense` | `novel_not_in_catalog` | 2 |
| `fusion` | `novel_not_in_catalog` | 1 |

Every `incomplete-splice_match` the computed mode calls `novel_in_catalog` is an isoform
SQANTI3 itself files under the `intron_retention` subcategory: its junctions are a
consecutive part of an annotated chain, but its terminal exon spans one more intron of that
transcript. SQANTI3 keeps some of these under ISM and others under NIC depending on the
reference transcript it associates; the computed mode treats every retained intron as a
novel structure. Those isoforms are then adjudicated instead of being held as fragments,
which is where the computed mode's higher recall comes from (FLAIR 0.748 to 0.778, IsoQuant
0.708 to 0.750, TALON 0.709 to 0.749 at a precision within 0.014 of the SQANTI3 mode). The
few `genic`, `antisense` and `fusion` multi-exon isoforms SQANTI3 does not adjudicate become
`novel_not_in_catalog` here when their locus is annotated on their strand, and are
adjudicated.

## Reproduce

```bash
WORK=/path/to/sqanti_sim/work MC=/path/to/multicaller/data GENOME=/path/to/chr22.fa \
  SJ=/path/to/h2h_sr/sr_1000000/SJ.out.tab OUT=/tmp/no_sqanti3 bash run.sh
```

Inputs are the public/simulated chr22 data of [../sqanti_sim](../sqanti_sim) and
[../multicaller](../multicaller). Result envelope:
[../results/no_sqanti3/metrics.json](../results/no_sqanti3/metrics.json).
