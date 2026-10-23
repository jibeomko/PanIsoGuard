# verdict_distribution: what a real run returns

Descriptive, not an accuracy claim: the verdict-class distribution PanIsoGuard writes on real data,
over all isoforms and over the novel ones (SQANTI3 NIC + NNC), with the share of probable calls that
carry the RT-switch annotation and the reasons behind held and undecided calls. The inputs are the
adjudicated tables that `benchmark/lrgasp_wtc11` (WTC11 PacBio and ONT, with and without Illumina)
and `benchmark/refbias_3way` (HG03516 and HG02717 Kinnex with the personal haplotypes, no short
reads) already wrote; nothing is re-run.

```
W=/mnt/Data   # the work roots of the two protocols
benchmark/verdict_distribution/score_distribution.py \
  --run lrgasp_pb_nosr=$W/pig_lrgasp/pb/adj_nosr.adjudicated.tsv --run lrgasp_pb_sr=$W/pig_lrgasp/pb/adj_sr.adjudicated.tsv \
  --run lrgasp_ont_nosr=$W/pig_lrgasp/ont/adj_nosr.adjudicated.tsv --run lrgasp_ont_sr=$W/pig_lrgasp/ont/adj_sr.adjudicated.tsv \
  --run hg03516_kinnex=$W/pig_hg03516/refbias_3way/armC/wgs.adjudicated.tsv --run hg02717_kinnex=$W/pig_hg02717/refbias_3way/armC/wgs.adjudicated.tsv \
  --emit-metrics benchmark/results/verdict_distribution/metrics.json
```

Envelope: `benchmark/results/verdict_distribution/metrics.json`.
