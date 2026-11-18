# The chr22 truth set re-simulated as ONT cDNA

Every threshold PanIsoGuard ships was set on PacBio HiFi-like reads (the [sqanti_sim](../sqanti_sim)
chr22 truth set, PBSIM3) and the ONT evidence so far was real data without a truth set
([gm12878_realdata](../gm12878_realdata), the ONT library of [lrgasp_wtc11](../lrgasp_wtc11)). This
protocol keeps the chr22 design unchanged (the same reduced annotation, the same 1,500 simulated
transcripts with the same requested counts, so the same hidden novel transcripts) and re-simulates the
reads as ONT cDNA with NanoSim (SQANTI-SIM's pre-trained human WTC11 cDNA guppy model, 100,000
reads, seed 1), aligns them with `minimap2 -ax splice`, and runs the rest of the
[sqanti3_filter_h2h](../sqanti3_filter_h2h) protocol on them: FLAIR, SQANTI3, the SQANTI3 rules /
ML / strict filters, 1 M simulated short-read pairs through STAR, and PanIsoGuard without the BAM,
with the BAM under the default (hifi) preset and with `--platform ont`. Nothing was re-tuned.

## What the ONT reads do to the caller

FLAIR reports 4,787 isoforms from the ONT reads (1,579 from the PacBio reads of the same design); of
the 4,182 non-FSM ones 320 are genuine novel isoforms and 3,862 false, so the genuine share is 8 %
against 84 % on PacBio. The filters are judged on that set, as in the PacBio protocol.

## Results (0.1.1 binary, ruleset 0.0.8)

| short reads | method | TP | FP | precision | recall | AUPRC | coverage | selective risk | AUPRC difference, 95 % CI |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| none | SQANTI3 rules filter | 286 | 1,275 | 0.183 | 0.894 | 0.172 | 1.000 | 0.313 | |
| none | SQANTI3 ML filter | 108 | 478 | 0.184 | 0.338 | 0.136 | 1.000 | 0.165 | |
| none | PanIsoGuard | 0 | 0 | n/a | 0.000 | 0.066 | 0.000 | n/a | vs rules [-0.120, -0.092] |
| none | PanIsoGuard + BAM | 222 | 70 | 0.760 | 0.694 | 0.549 | 0.070 | 0.240 | vs rules [+0.325, +0.430] |
| none | PanIsoGuard + BAM, `--platform ont` | 222 | 70 | 0.760 | 0.694 | 0.549 | 0.070 | 0.240 | vs rules [+0.325, +0.430] |
| 1 M pairs | SQANTI3 rules filter | 288 | 1,275 | 0.184 | 0.900 | 0.173 | 1.000 | 0.313 | |
| 1 M pairs | SQANTI3 rules + min_cov ≥ 3 | 269 | 98 | 0.733 | 0.841 | 0.628 | 1.000 | 0.036 | |
| 1 M pairs | SQANTI3 ML filter | 127 | 69 | 0.648 | 0.397 | 0.541 | 1.000 | 0.063 | |
| 1 M pairs | every junction ≥ 3 short reads | 299 | 114 | 0.724 | 0.934 | 0.681 | 1.000 | 0.032 | |
| 1 M pairs | PanIsoGuard | 178 | 24 | 0.881 | 0.556 | 0.538 | 0.166 | 0.037 | vs the one-line rule [-0.205, -0.080] |
| 1 M pairs | PanIsoGuard + BAM | 234 | 33 | 0.876 | 0.731 | 0.673 | 0.187 | 0.047 | vs the one-line rule [-0.069, +0.053] |
| 1 M pairs | PanIsoGuard + BAM, `--platform ont` | 234 | 33 | 0.876 | 0.731 | 0.673 | 0.187 | 0.047 | vs the one-line rule [-0.069, +0.053] |

Metrics as in [sqanti3_filter_h2h](../sqanti3_filter_h2h): positives are `NOVEL_CONFIRMED` /
`NOVEL_PROBABLE` (PanIsoGuard) or `Isoform` (SQANTI3 filters); coverage is the share of isoforms given a
decisive class and selective risk the wrong share among them; the CIs are 2,000 bootstrap resamples of
the non-FSM isoforms.

What holds from PacBio to ONT:

- The ordering. Without short reads PanIsoGuard with the BAM is the better filter (precision 0.76
  against 0.18 for SQANTI3's rules, AUPRC +0.33 to +0.43 above them); with short reads it is the most
  precise method (0.88) and ties the one-line short-read rule on AUPRC (0.673 against 0.681, the
  interval covering zero), as on PacBio.
- The presets. `--platform hifi` and `--platform ont` give the same class to every one of the 4,787
  isoforms, with and without short reads. Over the 635 novel junctions with at least 3 spanning reads
  the raw indel-near fraction exceeds the 0.5 gate at 354 (56 %), the chemistry alone (mean background
  0.70 per junction); the excess over each read's own indel rate (ruleset 0.0.7) exceeds it at none.
  The ONT preset is therefore inert here too, as on real GM12878 ONT reads.

What does not hold is the absolute level: precision 0.76 and 0.88 against 0.99 on PacBio. The ONT
FLAIR set is 92 % false novels, and the whole-chain read support lifts 70 of the 3,862 false chains
(1.8 %) next to 222 of the 320 genuine ones; on PacBio the same lift takes 4 of 140. The decoy
measurement of that lift on real ONT reads is in [chain_support](../chain_support).

## Reproduce

```bash
bash run.sh [/mnt/Data/PanIsoGuard_bench/pig_sqantisim_ont]   # PB=<the sqanti_sim root> for the chr22 design
```

`run.sh` copies the chr22 design from the PacBio run, simulates with `sqanti-sim.py sim --ont
--read_type cDNA` (NanoSim; `nanosim_compat/` makes SQANTI-SIM's scikit-learn-1.2 model pickles load
under a newer scikit-learn), and ends with `emit_ont.py`, which writes
`benchmark/results/sqanti_sim_ont/metrics.json`.
