# TUSCO single-isoform genes as real-data negatives

The one thing no real data set in this repository gave until now is a negative: an isoform
known to be false. TUSCO (Liu et al., *Nat Commun* 2026, doi:10.1038/s41467-026-72089-1) curates 46 human genes
that express one isoform across tissues (the "universal" panel; list in
[data/tusco_human.tsv](data/tusco_human.tsv), from ConesaLab/SQANTI3). By TUSCO's own scoring, a
transcript model in one of these genes that is not the reference transcript is a false
positive. So every NIC / NNC isoform a caller reports in a TUSCO gene is a presumed false
novel isoform, and what a filter does with it is a real-data precision check.

| data set | isoforms in TUSCO genes | novel (NIC / NNC) | PanIsoGuard confident | PanIsoGuard artifact | SQANTI3 rules filter keeps |
|---|---:|---:|---:|---:|---:|
| LRGASP WTC11 PacBio, + Illumina | 48 | 4 | **0** | 3 | 1 |
| LRGASP WTC11 PacBio, no short reads | 48 | 4 | **0** | 0 | 1 |
| LRGASP WTC11 ONT, + Illumina | 25 | 0 | – | – | – |
| HG03516 Kinnex (no short reads) | 26 | 2 | **0** | 0 | 2 |
| HG02717 Kinnex (no short reads) | 37 | 7 | 3 | 0 | 5 |
| total without short reads (one run per data set) | | **13** | **3** | 0 | **8** |

- Thirteen presumed-false novel isoforms in total: PanIsoGuard confirms 3 (23 %), SQANTI3's
  default rules filter keeps 8 (62 %). With short reads (PacBio) PanIsoGuard confirms none and
  rejects 3 of 4 as `ARTIFACT`; without short reads it holds the same four.
- The three it confirms are HG02717 chains with 3, 3 and 5 whole-chain reads lifted by
  `chain_support`, right at the ≥ 3 threshold; no sign, canonical motifs, kept by SQANTI3 too.
  TUSCO is a strong prior, not a truth: a universal single-isoform gene can still express an
  unannotated isoform in one lymphoblastoid line, so these three are suspicious, not proven
  wrong. They are also the case for a reader raising `chain_min_full_reads` on a deep sample.
- The non-canonical ones (two NNC models in ENSG00000185085 on both PacBio and HG02717, one
  with the RTS flag) are `ARTIFACT` with short reads and `UNCONFIRMED` without; SQANTI3's filter
  rejects them too. The novel-in-catalog model in ENSG00000135045 that SQANTI3 keeps on every
  data set has no novel junction to test and stays `UNRESOLVED`.
- Small by construction (46 genes), so it is reported as a check that points the same
  way as the simulation, not as a precision estimate.

## Reproduce

```bash
python3 score_tusco.py --tusco data/tusco_human.tsv \
  --run lrgasp_pb_sr=/path/pb/sq/iq_classification.txt,/path/pb/adj_sr.adjudicated.tsv,/path/pb/sq_rules/rules_RulesFilter_classification.txt \
  --run hg02717_kinnex=/path/armA/sq/iq_classification.txt,/path/armC/wgs.adjudicated.tsv,/path/armA/sq/rules/rules_RulesFilter_classification.txt \
  --emit-metrics ../results/tusco_negatives/metrics.json
```

Inputs are the tracked runs of [../lrgasp_wtc11](../lrgasp_wtc11) (`run.sh sqanti_filter adjudicate`)
and [../refbias_3way](../refbias_3way) (arm A SQANTI3 + rules filter, arm C `wgs`). The gene
list is redistributed unchanged from SQANTI3 (GPL-3.0) with attribution; the TUSCO portal is
<https://tusco.uv.es>.
