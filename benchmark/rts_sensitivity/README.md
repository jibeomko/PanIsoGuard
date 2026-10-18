# RT-switch policy: what the ruleset-0.0.5 choice moves

Since ruleset 0.0.5 SQANTI3's `RTS_stage` no longer demotes a caller-consensus or whole-chain
lift to `UNCONFIRMED`; the lift stands and the reason carries `;rt_switch_sign`. The 0.0.4
behaviour ("hold") can be read off the same outputs, so `rts_sensitivity.py` tabulates, per
tracked run, how many calls the choice decides and what they are where a truth or an external
filter exists. Why the two policies are both defensible, and why this is the one tracked
decision without a direct real-data test, is spelled out below.

| run | lifts carrying the RTS sign | truth / external view | annotate → hold |
|---|---:|---|---|
| SQANTI-SIM chr22, 1 M pairs + BAM | 26 | 26 genuine / 0 false | P 0.993 / R 0.745 / AUPRC 0.974 → 0.992 / 0.710 / 0.972 |
| SQANTI-SIM chr22, no short reads + BAM | 57 | 57 / 0 | 0.993 / 0.748 / 0.953 → 0.992 / 0.670 / 0.939 |
| end2end chr22 | 4 | 4 / 0 | genuine `NOVEL_PROBABLE` 19/20 → 15/20 |
| LRGASP WTC11 PacBio, + Illumina / none | 910 / 1,321 | validated targets among them: MTHFS; CENPK, MTHFS | the held ones would be `UNCONFIRMED`, never `ARTIFACT` |
| LRGASP WTC11 ONT, + Illumina / none | 121 / 222 | MTHFS, USF1; CENPK, MTHFS, RND2, USF1 | |
| HG03516 Kinnex (arm C) | 1,297 of 7,081 `NOVEL_PROBABLE` | all 1,297 are SQANTI3-filter `Artifact`s | |
| HG02717 Kinnex (arm C) | 1,357 of 8,368 | all `Artifact`s | |
| GM12878 ONT direct RNA | 0 | | |

`KNOWN_FRAGMENT` has the same annotation for ISM (`ism_chain_support;rts_sign`): 8 on
SQANTI-SIM, 222 / 52 on LRGASP PacBio / ONT, 311 / 171 on the two Kinnex samples.

## Reading

- **What the evidence covers.** The simulation produces no RT switching, so its 87 / 0 says only
  that the flag lands on genuine isoforms (it is computed from the reference sequence around
  the junction, so in a simulation every flag is a false sign). LRGASP adds the one real-data
  measurement: `RTS_stage` is TRUE on 12–18 % of the experimentally validated isoforms IsoQuant
  called, and the four validated novel isoforms that carry it (MTHFS, CENPK, USF1, RND2) are the
  ones the policy keeps, and all four are carried by ≥ 3 reads in the other platform's
  library too (CENPK 4 / 16, MTHFS 49 / 57, RND2 38 / 76, USF1 7 / 1 PacBio / ONT reads; USF1
  is the one with a single PacBio read).
- **The artifact side is still unmeasured.** [../lrgasp_wtc11](../lrgasp_wtc11)
  (`cross_library.py`) counts each lifted PacBio isoform's reads in the ONT library of the same
  cells and the reverse. PacBio lifts carrying the sign reproduce in ONT at about half the rate
  of the unflagged lifts in every bin of PacBio depth (own count 3–5: 4.5 % vs 12.6 %; 6–10:
  7.9 % vs 20.1 %; 11–20: 15.1 % vs 26.7 %; 21–50: 26.6 % vs 40.1 %; > 50: 40.7 % vs 51.1 %),
  which at first looks like half of them being library-specific artifacts. The control says
  otherwise: the annotated transcripts carrying the flag show the same deficit (FSM 6 / 11 /
  21 / 29 / 55 % flagged vs 18 / 29 / 39 / 50 / 72 % unflagged). The flag marks a
  repeat-flanked junction, the chain count demands exact coordinates, and ONT reads reach
  them less often across a repeat, so the deficit measures ONT alignment, not RT switching;
  from ONT to PacBio the flagged chains reproduce as well as or better than the unflagged
  ones. What remains unmeasured is how often a real RT-switching chimera that ≥ 3 reads carry
  end to end is flagged, which is what "hold" would protect against. Template switching is
  sequence-driven, so it can recur across molecules; read count alone does not rule it out.
- **What the choice moves.** On real Iso-Seq about one `NOVEL_PROBABLE` in six (1,297 / 7,081;
  1,357 / 8,368) rests on it, and every one of those is an isoform SQANTI3's rules filter
  discards. Under "hold" they are `UNCONFIRMED` (review), not `ARTIFACT`.
- **The trade, and why annotate.** The one real-data number is the false-sign rate: the flag
  is wrong about one time in six on isoforms known to exist, and a sign that weak should not
  by itself move a read-supported chain from "keep" to "review". Under "annotate" the
  confident set keeps the flagged lifts (about one `NOVEL_PROBABLE` in six on real Iso-Seq);
  under "hold" they go to review and the validated flagged isoforms go with them. Whichever
  policy is the default, the information is not lost: the reason string keeps the sign, so the
  other policy is one filter on `;rt_switch_sign` away, and a paper table can report both. The
  intra-priming sign still holds (no real-data measurement either way, and thousands of
  chain-supported isoforms per sample).

## Reproduce

```bash
~/miniconda3/envs/sqsim/bin/python rts_sensitivity.py --sqsim-work /path/to/sqanti_sim/work \
  --h2h /path/to/h2h_work --e2e /path/to/e2e_work --lrgasp /mnt/Data/pig_lrgasp \
  --threeway hg03516=/path/to/pig_hg03516/refbias_3way hg02717=/path/to/pig_hg02717/refbias_3way \
  --gm12878 /path/to/gm12878/rb_ont.adjudicated.tsv --emit-metrics ../results/rts_sensitivity/metrics.json
```
