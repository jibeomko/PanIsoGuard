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
  int consensus_min_callers = 999;       // disabled by default (single-caller runs)
  // BAM (mapping) axis
  int bam_min_mapq = 20;
  int bam_softclip_min_bp = 20;
  int bam_junction_window_bp = 10;
  double bam_max_low_mapq_frac = 0.5;       // > this fraction of low-MAPQ spanning reads -> mapping artifact
  double bam_max_supplementary_frac = 0.5;  // > this fraction supplementary/secondary -> mapping artifact
  // Per-axis enable switches (used by `ablate` to mask one axis at a time).
  bool use_short_read = true;
  bool use_mapping = true;
  bool use_noncanonical = true;
  bool use_rts = true;
  bool use_degradation = true;
  std::string sqanti3_version_target = "6.0";
  std::string ruleset_version = "builtin-0.0.1";
};

class RuleEngine {
 public:
  RuleEngine() = default;  // built-in defaults; usable with no config file

  // Load overrides from a TOML file (missing keys keep their defaults).
  static RuleEngine from_toml(const std::string& path);

  const RuleConfig& config() const { return cfg_; }

  // Map an assembled EvidenceVector to a Verdict (2-axis evaluation + projection).
  Verdict evaluate(const EvidenceVector& ev) const;

  // A copy of this engine with one evidence axis disabled, for ablation. Axis is
  // one of: short_read, mapping, noncanonical, rt_switch, degradation.
  RuleEngine with_axis_disabled(const std::string& axis) const;

 private:
  RuleConfig cfg_;
};

}  // namespace panisoguard
