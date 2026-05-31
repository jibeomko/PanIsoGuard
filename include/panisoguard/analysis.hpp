#pragma once

#include <map>
#include <string>
#include <vector>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/rules.hpp"
#include "panisoguard/sqanti.hpp"

namespace panisoguard {

// --- R1: non-redundancy of PanIsoGuard vs SQANTI3 on novel isoforms ----------
// "flagged" = called an artifact/false-novel. PanIsoGuard flag = ARTIFACT class;
// SQANTI3 flag = filter_result == "Artifact". Restricted to novel (NIC/NNC).
struct NonRedundancy {
  long both_flag = 0;        // agree: artifact
  long both_keep = 0;        // agree: kept
  long pig_only_flag = 0;    // PanIsoGuard flags, SQANTI3 keeps  (the non-redundant catch)
  long sqanti_only_flag = 0; // SQANTI3 flags, PanIsoGuard keeps
  long n_novel = 0;
  long n_sqanti_filter_available = 0;
  double jaccard = 0.0;      // |both_flag| / |union of flagged|
  double mcnemar_chi2 = 0.0; // continuity-corrected, 1 df
  double mcnemar_p = 1.0;    // two-sided
  // isoforms in the pig_only_flag cell, with their attributed mechanism
  std::vector<std::pair<std::string, std::string>> pig_only_examples;  // (isoform_id, mechanism)
};

NonRedundancy compute_nonredundancy(const std::vector<AdjudicationResult>& results,
                                    const SqantiTable& sqanti, std::size_t max_examples = 50);

// --- R4: per-axis ablation ---------------------------------------------------
// Re-evaluate the (already-assembled) evidence with one axis disabled and count
// how many isoforms change confidence class.
struct AxisAblation {
  std::string axis;
  long n_changed = 0;
  long n_total = 0;
  std::map<std::string, long> transitions;  // "FULL->ABLATED" -> count
};

AxisAblation compute_ablation(const std::vector<AdjudicationResult>& full_results,
                              const RuleEngine& base_engine, const std::string& axis);

}  // namespace panisoguard
