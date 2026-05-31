# The PanIsoGuard decision engine

This document describes the adjudication logic — the methodological contribution.
It is intentionally readable so that any verdict can be audited without reading
C++: all thresholds and the label projection live in
[`config/rules.default.toml`](../config/rules.default.toml), and the engine emits
a per-isoform `rule_trace`.

## Framing constraints (non-negotiable)

- PanIsoGuard reports **orthogonal attribution + confidence + provenance**.
- It **never** claims to beat the SQANTI3 random-forest filter or FLAIR2. The
  benchmark module reports *non-redundancy* (McNemar / 2×2 / Jaccard), not
  superiority.
- Validation uses **public + simulated** data (SQANTI-SIM primary; HG002/HPRC v1.1
  with HG002 confirmed out-of-graph; LRGASP/SIRV). A clinical cohort is at most an
  optional real-world demo / engineering scale test, never the source of
  mechanistic ground truth.

## Two-axis evidence space

Each novel isoform is scored on two axes; the engine then projects the (binned)
pair to one of the seven legacy confidence labels via the config's `projection`
table.

- **Axis A — novelty-support**: short-read `SJ.tab` corroboration, multi-caller
  fingerprint consensus, read support. Bins: `SUPPORTED` / `PARTIAL` / `UNSUPPORTED`.
- **Axis B — artifact-mechanism**: `none` / `reference_bias` / `variant_created` /
  `mapping_or_repeat` / `degradation` / `population_known`. Each mechanism gets a
  normalized `[0,1]` sub-score scaled by a config weight.

A mechanism whose input is absent (e.g. the variant axis with no haplotype FASTA)
is `not_evaluable` — a soft-skip, never silently counted as evidence for or
against, and recorded as such in the trace and provenance.

## Provenance and circularity

Every reclassification records its driving axis and a provenance tag
(`reference_only`, `haplotype_fasta`, `external_pangenome`, `population_freq`).
A verdict that rests solely on variants derived from the *same* RNA reads being
adjudicated is flagged `circularity_risk` and is **not** promoted to a hard
mechanistic class.

## Evidence gates (precede headline validation)

- **GATE-0**: the SQANTI-SIM variant-injection fork is built and verified
  (≥ 95 % motif recovery + negative control) before the variant axis (E3) is
  treated as a must-have deliverable.
- **GATE-1**: the pangenome reference is pinned to **HPRC v1.1 (freeze1)** with
  HG002 confirmed out-of-graph before the real-data rescue experiment (E6).
