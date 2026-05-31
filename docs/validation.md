# Validation

## What is verified

**Unit suite** (`ctest`, 28 cases / 159 assertions, htslib-only build):

- intron-chain fingerprint (order-independence, strand/chrom sensitivity)
- interval index RAII guard (build/query phase split)
- readers on controlled fixtures: SJ.tab (incl. strand discrimination), GTF +
  catalog, BED12, SQANTI3 (name-indexed, real category strings)
- cross-reader coordinate convergence (a BED isoform's fingerprint matches the
  GTF catalog chain)
- rule projection over the full 2-axis grid, incl. the mapping override and the
  variant circularity firewall
- BAM read-level features on a synthetic indexed BAM (spanning count, low-MAPQ,
  supplementary, soft-clip, indel-near)
- **GATE-0 motif core**: a controlled SNV fixture (`GG-AG` reference → `GT-AG`
  haplotype) is recognised as variant-created → `PAN_REF_RESCUED_FALSE_NOVEL`
- non-redundancy statistics (McNemar continuity-corrected χ², erfc p-value,
  Jaccard) and per-axis ablation re-evaluation

**Real-data plumbing / scale** (one cohort, 188,912 isoforms, GRCh38 + GENCODE v49):

- 100% SQANTI↔caller isoform-ID join
- `adjudicate` runs in seconds (priors+SJ) to ~1.6 min (with a 4.5 GB BAM)
- variant axis with reference == haplotype yields **0 false rescues** (the axis
  fires only on genuine reference-vs-haplotype motif differences)
- `combine` integrates two callers' differently-named running IDs by fingerprint
- two independent adversarial code reviews completed; findings fixed

## What is NOT yet established (and the plan)

The current `benchmark` / `ablate` outputs report **non-redundancy and per-axis
contribution**, i.e. *disagreement with SQANTI3 and which axis drives which calls*
— **not correctness**. Whether a reclassification is *right* needs ground truth.
That requires external tools not present in this environment
(`pbsim3`/`NanoSim`/`bcftools`/SQANTI-SIM); the plan, in priority order:

1. **SQANTI-SIM (E1/E2)** — simulate reads from a reference with known transcripts
   hidden to create labelled NIC/NNC/ISM truth; report per-class precision/recall,
   AUPRC for false-novel detection, and confidence calibration (ECE/Brier). This
   also lets `benchmark`/`ablate` report AUPRC / ΔAUPRC instead of only contingency
   counts.
2. **Variant-injection recovery (full GATE-0)** — inject variant-created canonical
   motifs, simulate variant-allele reads, align to the unmodified reference, and
   confirm ≥95% recovery + a negative control
   (see [../benchmark/variant_inject](../benchmark/variant_inject)).
3. **HG002 + HPRC v1.1** — matched-assembly truth for the variant/reference-bias
   axis on an out-of-graph individual.
4. **LRGASP (+ SIRV)** — multi-caller consensus and spike-in artifacts.

Thresholds in `config/rules.default.toml` are conservative defaults to be
**calibrated** by step 1 before any accuracy claim is made.
