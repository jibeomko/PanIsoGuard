# `panisoguard-report` — PDF report for adjudication output

A SQANTI3-style multi-page PDF that summarizes a `panisoguard adjudicate` run at a
glance. It is a **companion tool**: the C++ core stays dependency-free (htslib only),
and the report — like SQANTI3's own report step — is an optional Python post-process.
It **recomputes nothing**; it only visualizes the verdicts and the evidence the engine
already wrote.

## Install

```bash
pip install ./python          # or: pip install -e ./python  (dev)
# only extra dependency is matplotlib
```

## Use

```bash
panisoguard adjudicate --classification cls.txt --isoforms-gtf iso.gtf \
    --ref-gtf ref.gtf --out-prefix run        # produces run.attribution.jsonl + .provenance.log
panisoguard-report --prefix run               # writes run.report.pdf
# or, without installing:
python -m panisoguard.report --prefix run --out run.report.pdf
```

It reads `<prefix>.attribution.jsonl` (per-isoform verdict + evidence + rule_trace) and
`<prefix>.provenance.log` (run metadata + the authoritative class tally).

## What the report shows

| Page | Content |
|------|---------|
| 1. Executive summary | which evidence axes were in play, the headline novel-isoform verdict (trustworthy / needs-a-look / artifact / rescued), the confidence-class distribution, and the **circularity-firewall scorecard** (rescues held vs promoted). |
| 2. Two-axis decision space | the centerpiece **Axis A (novelty-support) × Axis B (mechanism)** projection grid — each cell is an isoform count, the border is the class the rulebook projects to — plus per-category novelty support, mechanism distribution, and short-read corroboration. |
| 3. Mechanism diagnostics | mechanism → resulting class, BAM mapping-evidence quality (gating vs reported-only), structural-signal availability (three-valued), and a register of mapping-down-weighted calls. |
| 4. Rescue & firewall | reference-bias rescue flow (variant / pangenome), promoted-vs-held outcomes, and the held-rescue ledger — the defensibility thesis. |
| 5. Consensus *(only if `--caller-support` was used)* | caller-support distribution, caller-count → class, and the consensus-driven promotions of UNKNOWN-support novels (capped at MEDIUM, never HIGH). |
| 6. SQANTI3 priors *(only if supplied)* | CAGE / polyA / NMD / polyA-motif pass-through — **banner: not used in the verdict**. |
| 7. Appendix | full provenance, most-frequent rule-trace lines, and a defensibility checklist (incl. a provenance-↔-JSONL tally reconciliation). |

Conditional panels degrade to a visible grey-hatch "axis not evaluable" callout rather
than an empty plot, so the *absence* of evidence is itself legible to a reviewer. A
single semantic color ramp is reused on every page (green→red confidence, blue =
reference-bias reclassification, grey = honest abstention, teal/grey-hatch = axis
on/absent).
