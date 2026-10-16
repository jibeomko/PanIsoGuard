#pragma once

#include <string>
#include <vector>

namespace panisoguard {

// Axis A: how well the isoform's novel junctions are independently supported.
enum class NoveltySupport { kUnknown, kUnsupported, kPartial, kSupported };

// Axis B: the dominant artifact mechanism attributed to a novel call. Tier-0
// (SQANTI priors + short-read corroboration) can reach none/noncanonical/
// rt_switch/degradation; mapping/variant/population require later evidence tiers.
enum class Mechanism {
  kNone,
  kNoncanonical,
  kRtSwitch,
  kDegradation,
  kMapping,
  kVariant,
  kPopulationKnown,
  kUnknown,
};

// Final per-isoform confidence class (the projection of the 2-axis grid). Names say the
// verdict, not the evidence; the evidence is in Verdict::verdict_reason and rule_trace.
enum class ConfidenceClass {
  kKnown,            // matches a known transcript (FSM); not adjudicated
  kNovelConfirmed,   // every novel junction corroborated, no artifact sign
  kNovelProbable,    // corroborated with one sign, or supported by caller agreement only
  kUnconfirmed,      // weak or partial evidence (also ISM: a truncated known transcript)
  kUnresolved,       // nothing to test, held by the circularity firewall, or a category
                     // PanIsoGuard does not adjudicate
  kReferenceBias,    // the novelty is explained by the sample's genome: a novel-vs-reference
                     // junction canonical on a personalized haplotype (variant axis) or
                     // realizable on a pangenome path (pangenome axis). The junction is real.
  kArtifact,         // evidence against
};

const char* to_string(NoveltySupport s);
const char* to_string(Mechanism m);
const char* to_string(ConfidenceClass c);

// Per-isoform evidence assembled by the adjudicator. `_evaluable` flags mark
// inputs that may be absent: an unevaluable axis is a soft-skip, never silently
// treated as evidence for or against.
struct EvidenceVector {
  std::string structural_category;
  bool is_novel = false;        // NIC or NNC
  bool chain_available = false;
  bool catalog_evaluable = false;  // chain present AND a reference catalog was supplied

  // short-read corroboration of the isoform's novel junctions; the known junctions are
  // counted too, to tell "no reads at the novel junction" from "no reads at this locus".
  bool sj_evaluable = false;
  int n_novel_junctions = 0;
  int n_novel_jx_sr_supported = 0;
  int n_known_junctions = 0;
  int n_known_jx_sr_supported = 0;

  // SQANTI3 QC priors (consumed, never recomputed). polyA_motif_found is only present
  // when SQANTI3 was run with a polyA motif list; it refines the perc_A intra-priming sign.
  bool rts_evaluable = false;    bool rts_stage = false;
  bool canon_evaluable = false;  bool noncanonical = false;
  bool percA_evaluable = false;  double perc_A_downstream_TTS = 0.0;
  bool polyA_motif_evaluable = false;  bool polyA_motif_found = false;

  // BAM read-level (mapping) axis: worst per-novel-junction fractions over reads
  // whose N op matches the junction, taken only over junctions with at least
  // bam_min_spanning reads (bam_n_junctions_gated of them).
  bool bam_evaluable = false;
  int bam_n_spanning_total = 0;
  int bam_n_junctions_gated = 0;
  double bam_max_frac_low_mapq = 0.0;
  double bam_max_frac_supplementary = 0.0;
  double bam_max_frac_softclip = 0.0;
  double bam_max_frac_indel_near = 0.0;

  // variant (reference-bias) axis: a novel junction that is non-canonical on the
  // reference but canonical on a personalized haplotype FASTA.
  bool variant_evaluable = false;
  bool variant_rescue = false;            // ALL novel junctions variant-created on ONE common haplotype
  int n_novel_jx_variant_created = 0;     // how many are created on >= 1 haplotype (per-junction view)
  bool haplotype_consistent = false;      // one haplotype carries the canonical motif at every
                                          // created junction (a chain lives on one haplotype)
  int n_novel_jx_variant_disrupted = 0;   // canonical on the reference, non-canonical on EVERY
                                          // haplotype: the canonical motif exists only on the reference
  bool variant_circular = false;  // haplotype provenance is RNA-derived/unknown (circular-risk)

  // pangenome (reference-bias) axis: every novel junction is realizable on a
  // pangenome graph haplotype path. Non-circular only when the junction set was
  // extracted from independent population assemblies (gated by provenance).
  bool pangenome_evaluable = false;
  bool pangenome_rescue = false;       // ALL novel junctions graph-supported, on one common path if known
  bool pangenome_circular = false;     // junction-set provenance is sample-derived/unknown (circular-risk)
  int n_novel_jx_pangenome = 0;
  bool pangenome_path_evaluable = false;  // every supporting junction row listed its path ids
  bool pangenome_path_consistent = false; // ...and one path carries all of them

  // multi-caller consensus (optional)
  bool consensus_evaluable = false;  int n_callers = 0;

  // whole-chain long-read support (--bam): reads whose alignment carries the isoform's
  // entire intron chain. Same reads the caller used, so methodological (not independent)
  // corroboration, like caller consensus: it can lift UNKNOWN support to PROBABLE only.
  bool chain_evaluable = false;
  int chain_n_full_reads = 0;
  int chain_n_partial_reads = 0;
};

struct Verdict {
  ConfidenceClass confidence = ConfidenceClass::kUnresolved;
  // The user's decision unit, one of: known, confirmed, reference_bias, unconfirmed,
  // artifact, not_adjudicated; and the one rule that decided it (e.g. partial_support,
  // no_short_reads, firewall_held, ism_intrapriming). Both are stable short codes.
  std::string verdict = "unconfirmed";
  std::string verdict_reason;
  NoveltySupport novelty_support = NoveltySupport::kUnknown;
  Mechanism primary_mechanism = Mechanism::kNone;
  bool circularity_flag = false;       // set when a verdict rests on RNA-derived variants
  std::vector<std::string> rule_trace; // human-readable record of which rules fired
};

}  // namespace panisoguard
