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
  double perc_A_degradation_threshold = 59.0;  // SQANTI intra-priming heuristic
  int consensus_min_callers = 999;       // disabled by default (single-caller runs)
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

 private:
  RuleConfig cfg_;
};

}  // namespace panisoguard
