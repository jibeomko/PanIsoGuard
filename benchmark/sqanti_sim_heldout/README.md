# sqanti_sim_heldout: the chr22 protocol on chromosomes it never saw

PanIsoGuard's thresholds (whole-chain reads ≥ 3, the mapping gates, the support axis) were set on
the SQANTI-SIM chr22 truth set of [../sqanti_sim](../sqanti_sim) and compared there against the
SQANTI3 filters and a one-line short-read rule ([../sqanti3_filter_h2h](../sqanti3_filter_h2h)).
This protocol repeats both, unchanged, on other chromosomes of GENCODE v49 with nothing re-tuned:
the same SQANTI-SIM design (`equal`, 1,500 deleted transcripts: 300 NIC, 300 NNC, 200 ISM; 100,000
PBSIM3 reads; seed 1), FLAIR, SQANTI3 6.0.1, and the two short-read settings of the main text
(none; 1 M wgsim pairs → STAR 2-pass). Scoring and bootstrap CIs are imported from the chr22
head-to-head so the numbers are directly comparable.

```
bash run.sh chr20          # -> /mnt/Data/pig_sqantisim_chr20/{work,h2h}
bash run.sh chr19
emit_heldout.py chr20=/mnt/Data/pig_sqantisim_chr20 chr19=/mnt/Data/pig_sqantisim_chr19 \
    --emit-metrics ../results/sqanti_sim_heldout/metrics.json
```

Envelope: `benchmark/results/sqanti_sim_heldout/metrics.json`.
