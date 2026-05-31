# GATE-0: variant-injection proof-of-concept

GATE-0 verifies that PanIsoGuard's variant (reference-bias) axis correctly
recognises a junction whose canonical splice motif exists only on the sample's
haplotype — i.e. a "novel" call that is really reference bias, not new splicing.

## Two layers

**1. Motif-reconstruction core (verified here, no external tools).**
The decisive question — *does a variant create/destroy a canonical GT-AG motif?* —
is implemented in `src/evidence/variant_motif.cpp` and verified deterministically
in `tests/unit/variant_motif_test.cpp` against a controlled fixture
(`tests/data/tiny/mini_ref.fa` vs `mini_hap.fa`) where a single SNV (G→T) turns a
non-canonical `GG-AG` donor into a canonical `GT-AG`. The provider returns
`kCreated`, and the rule engine maps it to `PAN_REF_RESCUED_FALSE_NOVEL` (or holds
it `AMBIGUOUS` under circular-risk provenance). This is 100% recovery of the
injected motif at the mechanism level.

**2. Read-simulation recovery (requires external tools; not run in this env).**
The full harness injects N variant-created and N variant-disrupted canonical
motifs at known GENCODE v49 loci, simulates variant-allele long reads, aligns them
to the UNMODIFIED GRCh38, runs the FLAIR→SQANTI pipeline, and confirms ≥95% of the
injected junctions are recovered (plus a negative control: the same loci without
injection do not produce the motif). This requires `pbsim3`/`NanoSim`, `minimap2`,
and `bcftools` (for `bcftools consensus` to build the personalized haplotype FASTA).
Status in the current environment: `minimap2` present; `pbsim3`, `NanoSim`, and
`bcftools` absent — install via:

```
conda install -c bioconda pbsim3 nanosim bcftools
```

then implement `run_truth.sh` (inject → consensus → simulate → align → call →
score recovery). Until then, GATE-0 is satisfied at the motif-reconstruction layer
(layer 1), which is the part PanIsoGuard itself implements; layer 2 validates the
simulation harness that produces truth labels, not the tool.
