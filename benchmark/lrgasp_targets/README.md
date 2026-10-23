# lrgasp_targets: the LRGASP validation targets as would-be isoforms

[../lrgasp_wtc11](../lrgasp_wtc11) scores the amplicon-validated LRGASP WTC11 targets that IsoQuant
reconstructed, which leaves the unsupported targets (the only amplicon-tested negatives) mostly
untested because the caller did not report them. This protocol removes the caller from the loop:
every multi-exon target chain of both Supplementary Data 15 sheets (`targets.py` of lrgasp_wtc11
writes them as GTFs) is classified by SQANTI3 against GENCODE v49 and adjudicated by PanIsoGuard
with the PacBio and the ONT library, with and without the Illumina junctions, as if a caller had
reported it. For the chains novel against v49: supported = positives, unsupported = negatives.
SQANTI3's default rules filter runs on the same classification.

```
bash run.sh /mnt/Data/pig_lrgasp                      # -> /mnt/Data/pig_lrgasp/targets_run
score_targets.py --lrgasp /mnt/Data/pig_lrgasp --emit-metrics ../results/lrgasp_targets/metrics.json
```

Envelope: `benchmark/results/lrgasp_targets/metrics.json`.
