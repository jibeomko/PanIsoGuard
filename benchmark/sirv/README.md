# SIRV spike-in control

The recognized spike-in standard as an independent control transcriptome (7 SIRV
loci, ~70 isoforms; Lexogen SIRV-Set4). Same incomplete-reference protocol as
[end2end](../end2end) but on the clean, dense SIRV reference and multi-contig.

## Run

```bash
# download SIRV-Set4 from Lexogen, then (genome = SIRV_isoforms_multi-fasta...fasta,
# annotation = ..._annotation_C_...gtf):
samtools faidx sirv.fa
python prep.py sirv.fa annotation.gtf work/
# PBSIM3 -> minimap2 (to sirv.fa) -> FLAIR collapse -> SQANTI3 -> panisoguard adjudicate
# then score by contig-aware intron-chain match to the hidden transcripts
```

## Result (SIRV-Set4, 61 multi-exon transcripts, 13 hidden; 3660 reads)

FLAIR exactly reconstructed 12/13 hidden transcripts and produced 228 isoforms. On
the SQANTI-novel set, scored by intron-chain truth:

```
            pred_genuine  pred_false  abstain
genuine            4           0          0
false              3          51          0

specificity(false)=0.944   genuine-recall(decisive)=1.000   precision(genuine)=0.571
```

**Honest interpretation.** Specificity (0.944) is lower than the GENCODE end-to-end
run (0.996) because SIRVs are deliberately dense, heavily-overlapping isoforms that
share splice junctions. The 3 false positives are FLAIR mis-collapses that form a
*novel combination of individually real, short-read-supported junctions* — the
short-read axis confirms each junction but cannot, on its own, judge whether the
combination is genuine. This is the documented limitation that multi-caller
consensus and read-level full-length support are designed to address (see
[../multicaller](../multicaller)); on the dense SIRV control it is the expected
hard case, and PanIsoGuard still removed 51/54 artifacts and rejected no genuine
novel.
