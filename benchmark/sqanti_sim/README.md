# SQANTI-SIM truth-based benchmark (per-class P/R + AUPRC)

The canonical controlled-novelty simulator [SQANTI-SIM](https://github.com/ConesaLab/SQANTI-SIM)
generates a labelled truth set: a chosen set of real transcripts is **deleted** from the
annotation the caller sees, so a caller must rediscover them as novel. This gives a
ground-truth P/R + AUPRC for PanIsoGuard's novel-isoform adjudication on a *real* caller
pipeline — complementing the in-house [controlled_truth](../controlled_truth) and
[end2end](../end2end) protocols with the standard tool.

## Pipeline (GENCODE v49 chr22)

```
SQANTI-SIM classif chr22.gtf                       # categorize transcripts
SQANTI-SIM design equal -nt 1500 --NIC 300 --NNC 300 --ISM 200
                                                   # delete 769 multi-exon -> genuine-novel truth
SQANTI-SIM sim --pb --pbsim (QSHMM-RSII)           # 98k PacBio HiFi reads
minimap2 -ax splice:hq -> FLAIR collapse           # 1579 isoforms (vs the REDUCED annotation)
SQANTI3 qc (vs reduced)
prep_truth_sj.py   # truth SJ.tab from the simulated transcripts' introns + chain truth map
panisoguard adjudicate --sj-tab truth.SJ.tab
score.py           # match FLAIR isoforms to truth by intron chain -> P/R + AUPRC
```

Truth labels (matched by exact intron chain): **genuine_novel** = chain of a deleted
simulated transcript; **known** = chain still in the reduced annotation; **false_novel** =
multi-exon chain matching no simulated transcript (a FLAIR/alignment artifact). The truth
SJ.tab carries the introns of all simulated (expressed) transcripts, so a genuine novel's
junctions are short-read-corroborated while an artifact's are not.

## Result (1500 simulated transcripts; 730 multi-exon genuine-novel, 652 known, 148 false)

Confusion (truth × PanIsoGuard confidence class):

```
truth          HIGH_KNOWN  HIGH_NOVEL  MED_NOVEL  LOW_PARTIAL  AMBIGUOUS  ARTIFACT
genuine_novel          0         385         45         154        146         0
known                652           0          0           0          0         0
false_novel            8           4          1          61          6        68
```

Genuine-novel detection (predicted genuine = `HIGH/MEDIUM_CONF_NOVEL`):

```
precision = 0.989   recall = 0.589   specificity(false rejected) = 0.966
AUPRC (genuine vs false) = 0.971   (positive rate 0.831 = a random ranking)
recall by category:  NNC 297/297 = 1.000 | NIC 133/279 = 0.477 | ISM 0/154 = 0.000
known -> HIGH_CONF_KNOWN: 652/652 = 1.000
```

(0.0.4: `PARTIAL` support no longer promotes to `MEDIUM_CONF_NOVEL`; before that fix the
same run read 4 false `MEDIUM`, precision 0.982, specificity 0.946, AUPRC 0.970.)

> **Read this with [sqanti3_filter_h2h](../sqanti3_filter_h2h).** The truth SJ.tab is an
> oracle (every expressed junction, 50 reads), and against the real baselines on this same
> truth set PanIsoGuard is **not** more accurate: the one-line rule "every junction has ≥ 3
> short reads" ranks as well, and at the default operating point has far higher recall.

## Interpretation

- **Precision 0.989, specificity 0.966, AUPRC 0.971** — when PanIsoGuard calls an isoform
  genuinely novel it is almost always right, and it rejects artifacts.
- **Zero genuine_novel → ARTIFACT**: it never wrongly rejected a real novel isoform.
- The **moderate recall (0.589) is abstention, not misclassification — but it has a cost**
  (the SQANTI3 rules filter and the one-line short-read rule keep most of these isoforms at
  ~0.93–0.99 precision; see [sqanti3_filter_h2h](../sqanti3_filter_h2h)):
  - **NNC** (novel-not-in-catalog) recall is **perfect (1.000)** — these carry genuinely
    novel junctions that the short-read axis corroborates.
  - **NIC** (novel-in-catalog) recall 0.477 — NIC novelty is a *combination* of known
    junctions; the isoforms with no individually-novel junction are held `AMBIGUOUS`
    (the engine has nothing to corroborate, so it does not guess — documented behavior).
  - **ISM** (incomplete-splice_match) "recall" 0.000 — ISMs are truncations of known
    transcripts *relative to the reduced annotation* and are assigned `LOW_CONF_PARTIAL`,
    not a confident novel call — yet each of these 154 chains exactly matches a deleted
    (genuine) transcript, so this is lost recall, not a labelling artifact.

So SQANTI-SIM shows the design intent working — confirm what an orthogonal axis can
verify (NNC), hold what it cannot (NIC-combinatorial / ISM-partial), distrust artifacts —
and also its price: the abstention costs recall that simpler filters keep.

## Reproduce

`run.sh` documents the exact commands; `prep_truth_sj.py` + `score.py` are the
truth/SJ.tab builder and scorer. Tools: SQANTI-SIM (lean env: biopython, bcbio-gff,
pysam, pbsim3, pbccs, gffread, minimap2, samtools, scikit-learn), FLAIR 3.x, SQANTI3 5.x.
