# The PanIsoGuard decision engine

This document describes the adjudication logic — the methodological contribution.
It is intentionally readable so any verdict can be audited without reading C++:
all thresholds live in [`config/rules.default.toml`](../config/rules.default.toml)
and the engine emits a per-isoform `rule_trace`.

## Framing constraints (non-negotiable)

- PanIsoGuard reports **orthogonal attribution + confidence + provenance**.
- It **never** claims to beat the SQANTI3 random-forest filter or FLAIR2; the
  head-to-head is in [benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h).
- It **does not recompute SQANTI3 QC descriptors** (TSS/TTS, ORF/NMD, polyA/splice
  motif): those are consumed as priors, while the
  verdict rests on independent path-level evidence
  ([relationship_to_sqanti3.md](relationship_to_sqanti3.md)).
- Validation uses **public + simulated** data (SQANTI-SIM; HG002/HPRC; LRGASP). A
  clinical cohort is at most a real-world demo / engineering scale test, never the
  source of mechanistic ground truth.

## Two-axis evidence space

Each *novel* isoform (SQANTI `novel_in_catalog` / `novel_not_in_catalog`) is scored
on two axes; the binned pair is projected to one confidence class.

- **Axis A — novelty support** (from short-read `SJ.tab` corroboration of the
  isoform's novel junctions): `SUPPORTED` / `PARTIAL` / `UNSUPPORTED`, or `UNKNOWN`
  when the axis is not evaluable. `UNSUPPORTED` (0 of N novel junctions corroborated)
  counts as evidence *against* only when the locus is covered: if none of the isoform's
  **known** junctions has short-read support either, the sample is too shallow there and
  the support stays `UNKNOWN` (reason `no_coverage`; since 0.0.6). When short-read
  corroboration is `UNKNOWN`, **multi-caller consensus** (`--caller-support`, from a
  `panisoguard combine` matrix) can stand in as a *weaker, methodological* form of
  support: a novel chain independently recovered by ≥ `consensus_min_callers` callers
  is corroborated, but only enough to reach `NOVEL_PROBABLE`, never `NOVEL_CONFIRMED`
  (caller agreement is not orthogonal experimental evidence).
- **Axis B — artifact mechanism**: `none` / `mapping_or_repeat` / `noncanonical` /
  `rt_switch` / `degradation` / `variant_created`. An axis whose input is absent is
  `not_evaluable` — a soft-skip, never silently counted for or against.

## Evidence axes and where they come from

| Mechanism | Source | Evaluable when |
|-----------|--------|----------------|
| short-read support | STAR `SJ.tab` exact junction match (strand-aware) vs the reference catalog | `--sj-tab` + `--ref-gtf` |
| mapping_or_repeat | BAM reads spanning the junction — the **low-MAPQ / supplementary / indel-near** fractions gate this mechanism (> 0.5 each), taken only over novel junctions with ≥ `min_spanning` (3) reads: with one or two reads a single indel is a fraction of 1.0, and 45 % of HG03516's mapping-artifact calls rested on < 3 reads before 0.0.6. `--platform ont` turns the indel-near gate off (it is HiFi-calibrated and fires on ~90 % of ONT junctions). | `--bam` |
| noncanonical / rt_switch | SQANTI3 priors `all_canonical` / `RTS_stage`; `noncanonical` also when a novel junction is canonical on the reference but non-canonical on **every** haplotype of the sample (`evidence.n_novel_jx_variant_disrupted`, reason suffix `;haplotype_disrupted`; needs `wgs`/`external` haplotype provenance) | always (from classification); the haplotype form with `--reference-haplotype` |
| degradation (intra-priming) | SQANTI3 `perc_A_downstream_TTS ≥ 60` **and** no polyA motif upstream of the 3' end (`polyA_motif_found`, present when SQANTI3 ran with `--polyA_motif_list`; without it perc_A alone decides and `provenance.log` says so). On HG03516, 70–80 % of perc_A-flagged ends lack any of the 12 polyA hexamers (clean ends: 11–17 %), so the sign is real, but the 20–30 % with a motif are kept | always; the motif column when supplied |
| variant_created | reference vs personalized haplotype FASTA splice-motif comparison | `--reference` + `--reference-haplotype` |
| caller consensus | a `panisoguard combine` matrix: how many independent callers recovered the isoform's intron chain (`n_callers`) | `--caller-support` |

## Mechanism priority and its rationale

When several mechanisms fire, one is reported as `primary_mechanism`; the order is
**mapping > noncanonical > rt_switch > degradation** (and the variant rescue, below,
preempts all of them). Rationale: **alignment reliability is the most upstream
concern** — if the spanning reads do not map confidently, the junction may be a
mapping artifact and downstream motif/QC signals are moot; non-canonical motif is the
next strongest structural signal; RT-switching and intra-priming degradation are
softer biochemical signals. Only the first-matching mechanism is reported, both as
`primary_mechanism` and in `rule_trace` (the check stops at the first match); the other
flags stay visible in the `evidence` block of `attribution.jsonl`.

## Projection to confidence classes (7 reachable)

1. **Reference-bias rescue (highest priority).** A novel-vs-linear-reference
   junction explained by reference bias is `REFERENCE_BIAS`, via either:
   - *Pangenome axis* (`--pangenome-junctions`): **all** of the isoform's novel
     junctions are realizable on a pangenome graph haplotype path (mechanism
     `population_known`) — consistent with reference bias, not proof the splice is
     used. Independent only if the junction set is from population assemblies, which
     PanIsoGuard cannot verify, so it is gated by the **same circularity firewall**
     as the variant axis: `--pangenome-provenance population`/`external` promote,
     while `sample_derived`/`unknown` (the default) are held `UNRESOLVED`. It is
     evaluated first, so an *independent* graph rescues even when the sample's own
     haplotype provenance is circular-risk. **Validated** on the real HPRC v1.1 chr22
     graph (GATE-1; see [benchmark/pangenome](../benchmark/pangenome)).
   - *Variant axis* (`--reference-haplotype`): the junction is non-canonical on the
     reference but canonical on a personalized haplotype (mechanism
     `variant_created`). **Circularity firewall:** if the haplotype provenance is
     RNA-derived/unknown this is *not* independent evidence, so it is held
     `UNRESOLVED` with `circularity_flag=true` and never promoted; only
     `wgs`/`external` provenance promotes.

   Both rescues need **all** of the isoform's novel junctions explained, **by one
   haplotype or one graph path**: a transcript is one molecule from one allele, so
   junction 1 fixed on haplotype A plus junction 2 fixed on haplotype B explains
   nothing. The adjudicator intersects the per-junction haplotype sets
   (`evidence.haplotype_consistent`); when every junction is created but the
   intersection is empty the call is held `UNRESOLVED` with reason
   `haplotype_inconsistent`. The pangenome file may list the path ids per junction
   (optional 6th column); when it does, the same intersection applies
   (`evidence.pangenome_path_consistent`, reason `path_inconsistent`), and when it does
   not the rescue rests on the counts alone (`pangenome_path_consistent = null`). The
   per-junction count is in `evidence.n_novel_jx_variant_created` /
   `n_novel_jx_pangenome`. A rescue is not blocked by a BAM mapping sign (that blocked
   a correct rescue in [benchmark/refbias_phasing](../benchmark/refbias_phasing)), but
   the sign is appended to the reason (`haplotype_variant;mapping_sign`) so the call
   can be reviewed. The rescue is allele-agnostic: canonical on *either* haplotype
   suffices, and the reads may come from the other one (about 1 in 10 rescued
   junctions in HG03516; see [benchmark/refbias_3way](../benchmark/refbias_3way)).
2. Otherwise project (Axis A × Axis B):
   - `SUPPORTED × none → NOVEL_CONFIRMED`
   - `SUPPORTED × mechanism → NOVEL_PROBABLE`
   - `PARTIAL × any → UNCONFIRMED` (an uncorroborated novel junction leaves the
     chain unconfirmed; before 0.0.4 `PARTIAL × none` promoted to `NOVEL_PROBABLE`,
     which produced false positives — see
     [benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h))
   - `UNSUPPORTED × none → UNCONFIRMED`
   - `UNSUPPORTED × mechanism → ARTIFACT`
   - `UNKNOWN` support → `UNRESOLVED`, except:
     - `UNKNOWN × strong mapping artifact → UNCONFIRMED` (reason
       `mapping_without_short_reads`): the mapping sign preempts consensus and whole-chain
       support, but without short reads it holds the call instead of rejecting it. Before
       ruleset 0.0.5 this cell was `ARTIFACT`; on SQANTI-SIM 11 of the 31 isoforms it
       rejected were genuine, and on real ONT data the HiFi-calibrated gate fires on most
       junctions; and
     - `UNKNOWN × caller-consensus (≥ consensus_min_callers, no mapping artifact) →
       NOVEL_PROBABLE` (mechanism `none`, or an RT-switch sign, noted as
       `caller_consensus;rt_switch_sign`) or `UNCONFIRMED` (`noncanonical` or
       `degradation`, reason `consensus_with_<mechanism>`). Multi-caller agreement corroborates the chain methodologically
       when the short-read axis is silent; a strong mapping artifact still preempts
       it. It applies only to a chain that was examined: if the caller's chain did
       not join (`caller_chain_absent`) or no `--ref-gtf` was given, the count is
       reported in `evidence.n_callers` but the call stays `UNRESOLVED`. The matrix
       is joined by intron-chain fingerprint (the `fingerprint` column `combine`
       writes since 0.0.6), not by native id, because two callers may reuse one name
       for different chains.
     - `UNKNOWN × whole-chain long-read support (≥ chain_min_full_reads, no mapping
       artifact, no consensus) → NOVEL_PROBABLE` (mechanism `none`, or an RT-switch sign,
       `chain_support;rt_switch_sign`) or `UNCONFIRMED` (`noncanonical` or `degradation`,
       `chain_support_with_<mechanism>`), since 0.0.7. The RT-switch sign stopped demoting
       methodological support in ruleset 0.0.5: on SQANTI-SIM every isoform it held was
       genuine (26 / 0 with short reads, 57 / 0 without), and SQANTI3's `RTS_stage` is TRUE
       on 12–18 % of LRGASP's experimentally validated isoforms (what the choice moves per
       dataset: [benchmark/rts_sensitivity](../benchmark/rts_sensitivity)). The intra-priming sign still
       holds: the simulator does not model it (4 / 0 there says nothing), and on real Iso-Seq
       it marks 2,000–3,000 chain-supported isoforms per sample that SQANTI3's filter also
       rejects ([benchmark/refbias_3way](../benchmark/refbias_3way)). With `--bam`, the adjudicator counts the
       primary reads whose alignment carries the isoform's entire intron chain
       (`evidence.chain_n_full_reads`; `chain_n_partial_reads` are shorter reads
       consistent with it). These are the caller's own reads, so like consensus it is
       methodological support and never reaches `NOVEL_CONFIRMED`; it never overrides a
       short-read verdict (`UNSUPPORTED`, `PARTIAL`). On SQANTI-SIM chr22 it lifts
       146 genuine and 0 false FLAIR isoforms with 1 M short-read pairs (546 / 4 without), but the false
       chains of Bambu, ESPRESSO and TALON are read-supported too (alignment-level
       artifacts), so it separates caller-level from alignment-level false chains only
       ([benchmark/chain_support](../benchmark/chain_support)). Coordinates are exact: on
       real ONT direct-RNA 94.5 % of splice sites among near-matching reads are exact and a
       ± 5 bp tolerance would add 7 % of chains, so none is shipped
       ([benchmark/gm12878_realdata](../benchmark/gm12878_realdata) § 4). With
       `--validation-bam` (a second library the caller did not see) the reads are counted
       there instead; the evidence is then independent (`evidence.chain_independent`,
       reason `chain_support;independent_library`), but the class is still `PROBABLE`. On
       LRGASP WTC11, with the other platform's library as the independent one, 13 of 17
       (PacBio → ONT) and 14 of 16 (ONT → PacBio) validated novel isoforms are carried by
       ≥ 3 independent reads (the misses have 0–2), and within genes the other library
       covers the novel lifts reproduce like known transcripts; but a second library attests
       expression, not correct splicing, so the evidence is recorded rather than promoted ([benchmark/lrgasp_wtc11](../benchmark/lrgasp_wtc11)).
       Same chain+catalog requirement as consensus.
       `panisoguard decoys` measures the lift's false-fire rate on a user's own BAM: it
       recombines each gene's annotated transcripts at a shared intron into
       novel-in-catalog chains and reports the fraction that ≥ `chain_min_full_reads` reads
       carry end to end (an upper bound; some decoys are real, unannotated isoforms).
- Pass-through: `full-splice_match → KNOWN`; `incomplete-splice_match →
  KNOWN_FRAGMENT` (verdict `known_fragment`, reason `ism_chain_support`) when
  ≥ `chain_min_full_reads` long reads carry its whole intron chain and the 3' end shows
  no intra-priming sign (an RTS sign is noted, `;rts_sign`), since ruleset 0.0.5;
  otherwise `UNCONFIRMED` (a truncated known transcript; the intra-priming / RTS sign is
  reported in the reason, `ism_intrapriming` / `ism_rts`, but does not change the class,
  because the splice structure is known and only the 3' end is unreliable); non-targeted
  categories (genic, antisense, fusion, intergenic, genic_intron) → `UNRESOLVED` with
  verdict `not_adjudicated`.

## `verdict` and `verdict_reason` (since 0.0.6)

The class is the detailed label; two more columns give the user's decision unit and the
one rule that decided it:

| `verdict` | classes | meaning |
|---|---|---|
| `known` | KNOWN | matches a known transcript |
| `known_fragment` | KNOWN_FRAGMENT | a truncated known transcript whose whole intron chain the long reads carry (since ruleset 0.0.5); a known structure, not a novel isoform |
| `confirmed` | NOVEL_CONFIRMED, NOVEL_PROBABLE | keep (PROBABLE via `caller_consensus` or `chain_support` is methodological, not experimental, support) |
| `reference_bias` | REFERENCE_BIAS | the junction is real; it looked novel or non-canonical because of the reference |
| `unconfirmed` | UNCONFIRMED, UNRESOLVED | review |
| `artifact` | ARTIFACT | discard |
| `not_adjudicated` | UNRESOLVED | a SQANTI3 category PanIsoGuard does not judge |

`verdict_reason` codes: `full_splice_match`, `ism`, `ism_intrapriming`, `ism_rts`,
`ism_chain_support[;rts_sign]`,
`category_<category>`, `pangenome_path`, `haplotype_variant` (either `;mapping_sign`),
`firewall_held`, `haplotype_inconsistent`, `path_inconsistent`, `all_junctions_supported`
(a `noncanonical` reason carries `;haplotype_disrupted` when the sign came from the
haplotypes rather than SQANTI3),
`supported_with_<mechanism>`,
`partial_support[_with_<mechanism>]`, `unsupported_no_sign`, `unsupported_<mechanism>`,
`mapping_without_short_reads`, `caller_consensus[;rt_switch_sign]`,
`consensus_with_<mechanism>`, `chain_support[;rt_switch_sign]`,
`chain_support_with_<mechanism>` (each `;independent_library` with `--validation-bam`), and the
`UNKNOWN`-support reasons `no_short_reads`, `no_coverage`, `no_novel_junction`,
`caller_chain_absent`. The classes are **not** an ordered scale: on SQANTI-SIM truth
`UNRESOLVED` (mostly NIC with no novel junction to test) has precision 0.96 while
`UNCONFIRMED` novel isoforms have 0.36–0.56 (see
[benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h), `per_class.py`).

**`ARTIFACT` semantics (read carefully).** `ARTIFACT` denotes *(UNSUPPORTED × any
mechanism)* — short reads were checked at a covered locus, none supports the novel
junction, and an artifact sign is present — i.e. a **combined confidence × mechanism**
signal, **not** a standalone claim that the call is spurious. Since ruleset 0.0.5 no
isoform is `ARTIFACT` without short reads (the mapping sign alone holds it
`UNCONFIRMED`). A novel isoform with an unusual motif but no short-read backing lands
here; users should cross-check `ARTIFACT` calls against the original reads before hard
filtering.

**NIC behaviour.** A `novel_in_catalog` isoform's novel junction is a *novel pairing
of known splice sites*, i.e. a novel intron not present in the catalog, so it IS
counted as a novel junction and DOES receive short-read / BAM / variant evaluation
(it is not auto-`UNRESOLVED`). The residual `UNRESOLVED` set is non-targeted categories
plus isoforms whose every junction is an already-known intron (no novel junction to
evaluate). `rule_trace` and `verdict_reason` record which `UNKNOWN` reason applies
(`caller_chain_absent` / `no_short_reads` / `no_novel_junction` / `no_coverage`).

## Provenance and circularity

Every reclassification records its driving axis and (for the variant axis) a
provenance/circularity flag. A verdict resting solely on variants derived from the
same RNA reads being adjudicated is flagged (`circularity_flag = true`) and is **not**
promoted to a hard class. `provenance.log` reports which axes were evaluable.

## Evidence gates

- **GATE-0**: the variant-injection proof-of-concept (does a variant create/destroy
  a canonical GT-AG?) is verified deterministically at the motif level
  (`tests/unit/variant_motif_test.cpp`; see [benchmark/variant_inject](../benchmark/variant_inject)).
- **GATE-1 (met)**: the pangenome reference is the **HPRC v1.1 Minigraph-Cactus
  GRCh38** graph, with HG002/NA24385 confirmed **out-of-graph** (absent from the 47
  assembly samples → non-circular). On chr22, `vg deconstruct` → population deletions →
  graph-supported junctions: the rescue produces **0 false rescues** on real FLAIR novel
  junctions, **fires correctly** when a novel intron is a real population deletion
  (population provenance), and is **held by the circularity firewall** under
  circular-risk provenance — see [benchmark/pangenome](../benchmark/pangenome) and
  [benchmark/results/pangenome](../benchmark/results/pangenome).

## Future work

- **In-process pangenome tier**: the pangenome axis consumes graph-supported junctions
  pre-extracted with vg/rpvg (`--pangenome-junctions`) and emits `REFERENCE_BIAS`
  (mechanism `population_known`); traversing a GBZ graph in-process is not implemented.
- **Table-driven projection**: move the (Axis A × Axis B) grid fully into the TOML.
- **SQANTI RF concordance**: optionally annotate `rule_trace` with agreement vs the
  SQANTI3 `filter_result`, reinforcing the orthogonal-attribution framing.
- **Variant-disrupted junctions as a rescue-side signal**: since 0.0.7 a motif that is
  canonical on the reference only is a `noncanonical` sign (below); whether such
  junctions should also block a `KNOWN` pass-through is not decided.
