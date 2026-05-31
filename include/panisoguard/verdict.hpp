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

// Final per-isoform confidence class (the projection of the 2-axis grid).
enum class ConfidenceClass {
  kHighConfKnown,
  kHighConfNovel,
  kMediumConfNovel,
  kLowConfPartial,
  kPanRefRescuedFalseNovel,     // reference bias: a novel-vs-linear-reference junction
                                // explained by a personalized haplotype (variant axis)
                                // or realizable on a pangenome graph path (pangenome axis)
  kAmbiguous,
  kArtifact,
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

  // short-read corroboration of the isoform's novel junctions
  bool sj_evaluable = false;
  int n_novel_junctions = 0;
  int n_novel_jx_sr_supported = 0;

  // SQANTI3 QC priors (consumed, never recomputed)
  bool rts_evaluable = false;    bool rts_stage = false;
  bool canon_evaluable = false;  bool noncanonical = false;
  bool percA_evaluable = false;  double perc_A_downstream_TTS = 0.0;

  // BAM read-level (mapping) axis: worst per-novel-junction fractions over reads
  // whose N op matches the junction.
  bool bam_evaluable = false;
  int bam_n_spanning_total = 0;
  double bam_max_frac_low_mapq = 0.0;
  double bam_max_frac_supplementary = 0.0;
  double bam_max_frac_softclip = 0.0;
  double bam_max_frac_indel_near = 0.0;

  // variant (reference-bias) axis: a novel junction that is non-canonical on the
  // reference but canonical on a personalized haplotype FASTA.
  bool variant_evaluable = false;
  bool variant_rescue = false;
  bool variant_circular = false;  // haplotype provenance is RNA-derived/unknown (circular-risk)

  // pangenome (reference-bias) axis: a novel junction realizable on a pangenome
  // graph haplotype path (population-level reference data; non-circular).
  bool pangenome_evaluable = false;
  bool pangenome_rescue = false;
  int n_novel_jx_pangenome = 0;

  // multi-caller consensus (optional)
  bool consensus_evaluable = false;  int n_callers = 0;
};

struct Verdict {
  ConfidenceClass confidence = ConfidenceClass::kAmbiguous;
  NoveltySupport novelty_support = NoveltySupport::kUnknown;
  Mechanism primary_mechanism = Mechanism::kNone;
  bool circularity_flag = false;       // set when a verdict rests on RNA-derived variants
  std::vector<std::string> rule_trace; // human-readable record of which rules fired
};

}  // namespace panisoguard
