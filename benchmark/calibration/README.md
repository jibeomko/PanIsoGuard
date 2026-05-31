# Confidence-class calibration (ECE / Brier)

PanIsoGuard emits ordinal confidence **classes**, not native probabilities.
`calibrate.py` maps each class to an empirical `P(genuine)` learned on a held-out
*calibration* half of a labelled set, then reports the Brier score and Expected
Calibration Error (ECE) on the *test* half — so the calibration is never evaluated
on the data that defined it.

```bash
python calibrate.py <truth.tsv: isoform<TAB>genuine|false> <prefix.adjudicated.tsv>
```

## Results

**controlled_truth** (clean short-read axis, GENCODE v49 chr22; truth = hidden vs
fabricated):

```
class -> P(genuine):  HIGH_CONF_NOVEL 1.000 | AMBIGUOUS 0.996 | LOW_CONF_PARTIAL 0.000
test n=1648   Brier=0.0000   ECE=0.0013
```

**end-to-end** (real FLAIR + SQANTI3 output; truth = intron-chain match to hidden
transcripts):

```
class -> P(genuine):  HIGH_CONF_NOVEL 1.000 | AMBIGUOUS 0.625 | LOW_CONF_PARTIAL 0.000 | ARTIFACT 0.000
test n=143   Brier=0.0121   ECE=0.0114
```

The empirical genuine-rate per class is consistent between the calibration and test
halves (low ECE/Brier), i.e. the confidence classes are **well-calibrated**. Note
that `AMBIGUOUS` is not a 50/50 "unknown": in these sets it is dominated by genuine
NIC-like isoforms the engine honestly held (no novel junction to corroborate), so
its empirical `P(genuine)` is high — consistent with the documented behaviour.

(Small test-N on the end-to-end set; the controlled set is the high-N check.)
