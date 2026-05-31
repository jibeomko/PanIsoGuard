# Controlled-truth E1 (incomplete-reference protocol)

A fast, fully-local correctness check of PanIsoGuard's **adjudication logic** — no
read simulator or upstream caller required. Truth is generated directly from a
reference annotation, so every isoform's genuine/false label is known.

## Method

From one chromosome of a reference GTF (e.g. GENCODE):

- **Genuine novels** = ~20% of multi-exon transcripts *hidden* from the catalog
  handed to PanIsoGuard. They are real, so their novel junctions appear in the
  short-read `SJ.tab`.
- **False novels** = kept transcripts with one internal exon boundary shifted,
  creating a junction that is neither annotated nor read-supported.

A correct adjudicator should classify genuine novels as `HIGH/MEDIUM_CONF_NOVEL`
and fabricated ones as `LOW_CONF_PARTIAL` / `ARTIFACT`, and may honestly **abstain**
(`AMBIGUOUS`) on genuine isoforms whose every junction is already known (NIC-like:
no novel junction to corroborate).

## Run

```bash
awk -F'\t' '$1=="chr22"' gencode.vNN.annotation.gtf > chr22.gtf
python gen.py chr22.gtf work/
panisoguard adjudicate \
  --classification work/classification.tsv --isoforms-gtf work/caller.gtf \
  --ref-gtf work/reduced_catalog.gtf --sj-tab work/real.SJ.tab --out-prefix work/e1
python score.py work/truth.tsv work/e1.adjudicated.tsv
```

## Result (GENCODE v49 chr22, 11,030 multi-exon transcripts → 2,206 genuine + 1,070 false)

```
         pred_genuine  pred_false  abstain
genuine         1220           0       986
false              0        1068         2

Decisive calls (n=2288): precision(genuine)=1.000 specificity(false)=1.000 F1=1.000 accuracy=1.000
Abstain (AMBIGUOUS): genuine 986/2206 (44.7%), false 2/1070 (0.2%)
```

**Zero misclassifications on decisive calls**: every short-read-corroborated genuine
novel was kept and every read-unsupported fabricated junction was rejected. The
44.7% genuine abstention is by design — those are NIC-like isoforms whose every
junction is already annotated, so there is no novel junction to corroborate; the
engine holds them `AMBIGUOUS` rather than guessing.

## Scope (honest)

This isolates the **short-read corroboration axis** (the primary novelty signal)
against controlled truth; it does not exercise the BAM/variant axes or the upstream
callers. The full read-simulation E1 (PBSIM3/NanoSim → minimap2 → FLAIR → SQANTI3
→ PanIsoGuard), which also validates the simulator and the callers, is the
complementary end-to-end validation (see [../variant_inject](../variant_inject) and
[../../docs/validation.md](../../docs/validation.md)).
