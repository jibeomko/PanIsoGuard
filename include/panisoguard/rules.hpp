#pragma once

#include <string>

#include "panisoguard/verdict.hpp"

namespace panisoguard {

// All tunable thresholds live here and are loaded from a human-readable TOML so
// the decision logic is auditable without reading C++ (the projection itself is
// documented in docs/decision_engine.md and echoed per-isoform in rule_trace).
struct RuleConfig {
  int sj_min_uniq_reads = 3;             // STAR SJ.tab n_uniq to corroborate a novel junction
  bool sj_require_canonical_motif = true;
  double perc_A_degradation_threshold = 60.0;  // SQANTI3 intra-priming default: keep <60, flag >=60
  int consensus_min_callers = 2;         // a novel chain recovered by >= this many callers
                                         // is consensus-supported (when --caller-support given)
  // BAM (mapping) axis
  int bam_min_mapq = 20;
  int bam_min_spanning = 3;                 // a novel junction's read fractions count only with >= this many
                                            // spanning reads; below it the mapping axis is not evaluated there
  std::string platform = "hifi";            // "hifi" (defaults) or "ont" (indel-near gate off; see apply_platform_preset)
  int bam_junction_window_bp = 10;
  double bam_max_low_mapq_frac = 0.5;       // > this fraction of low-MAPQ spanning reads -> mapping artifact
  double bam_max_supplementary_frac = 0.5;  // > this fraction supplementary/secondary -> mapping artifact
  double bam_max_indel_near_frac = 0.5;     // > this fraction of spanning reads with an indel adjacent to the
                                            // junction -> mapping artifact (an alignment-ambiguous indel, not a real
                                            // intron). Junction-proximal (bam_junction_window_bp); validated on chr22.
  // Pangenome (reference-bias) axis
  int pangenome_min_haplotypes = 1;         // a novel junction must be on >= this many graph haplotypes to rescue
  // Whole-chain long-read support (--bam): reads carrying the isoform's entire intron chain
  int chain_min_full_reads = 3;             // >= this many -> UNKNOWN support becomes PROBABLE (never CONFIRMED)
  std::string sqanti3_version_target = "6.0";
  // Bump (and the [meta] name in config/rules.default.toml) whenever a threshold default or
  // the projection in RuleEngine::evaluate() changes; the history is in CHANGELOG.md.
  std::string ruleset_version = "builtin-0.0.6";
};

class RuleEngine {
 public:
  RuleEngine() = default;  // built-in defaults; usable with no config file

  // Load overrides from a TOML file (missing keys keep their defaults).
  static RuleEngine from_toml(const std::string& path);

  const RuleConfig& config() const { return cfg_; }

  // The effective thresholds as "key=value ..." (one line of provenance.log), so a run stays
  // reproducible even if the config file changes later.
  std::string thresholds_summary() const;

  // Sequencing-platform preset applied after the config: "hifi" keeps the defaults, "ont"
  // disables the indel-near mapping gate (HiFi-calibrated; fires on ~90 % of ONT junctions).
  void apply_platform_preset(const std::string& platform);

  // Map an assembled EvidenceVector to a Verdict (2-axis evaluation + projection).
  Verdict evaluate(const EvidenceVector& ev) const;

 private:
  RuleConfig cfg_;
};

}  // namespace panisoguard
