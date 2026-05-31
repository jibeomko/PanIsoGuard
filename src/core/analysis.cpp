#include "panisoguard/analysis.hpp"

#include <cmath>
#include <unordered_map>

namespace panisoguard {

NonRedundancy compute_nonredundancy(const std::vector<AdjudicationResult>& results,
                                    const SqantiTable& sqanti, std::size_t max_examples) {
  // isoform_id -> SQANTI filter_result (only for rows that carry one)
  std::unordered_map<std::string, std::string> filter_of;
  filter_of.reserve(sqanti.records.size());
  for (const auto& r : sqanti.records) {
    if (!r.filter_result.empty()) filter_of.emplace(r.isoform, r.filter_result);
  }

  NonRedundancy nr;
  for (const auto& res : results) {
    // Restrict to the contested novel set.
    const bool is_novel = res.structural_category == "novel_in_catalog" ||
                          res.structural_category == "novel_not_in_catalog";
    if (!is_novel) continue;
    ++nr.n_novel;

    auto it = filter_of.find(res.isoform_id);
    if (it == filter_of.end()) continue;  // no SQANTI verdict to compare against
    ++nr.n_sqanti_filter_available;

    const bool pig_flag = res.verdict.confidence == ConfidenceClass::kArtifact;
    const bool sqanti_flag = (it->second == "Artifact");

    if (pig_flag && sqanti_flag) ++nr.both_flag;
    else if (!pig_flag && !sqanti_flag) ++nr.both_keep;
    else if (pig_flag && !sqanti_flag) {
      ++nr.pig_only_flag;
      if (nr.pig_only_examples.size() < max_examples) {
        nr.pig_only_examples.emplace_back(res.isoform_id, to_string(res.verdict.primary_mechanism));
      }
    } else {
      ++nr.sqanti_only_flag;
    }
  }

  const long union_flag = nr.both_flag + nr.pig_only_flag + nr.sqanti_only_flag;
  nr.jaccard = union_flag > 0 ? static_cast<double>(nr.both_flag) / union_flag : 0.0;

  // McNemar with continuity correction on the discordant cells (b, c), 1 df.
  const long b = nr.pig_only_flag, c = nr.sqanti_only_flag;
  if (b + c > 0) {
    const double diff = std::fabs(static_cast<double>(b - c)) - 1.0;
    nr.mcnemar_chi2 = diff > 0 ? (diff * diff) / static_cast<double>(b + c) : 0.0;
    // two-sided p for chi-square_1df: P(X > x) = erfc(sqrt(x/2))
    nr.mcnemar_p = std::erfc(std::sqrt(nr.mcnemar_chi2 / 2.0));
  }
  return nr;
}

AxisAblation compute_ablation(const std::vector<AdjudicationResult>& full_results,
                              const RuleEngine& base_engine, const std::string& axis) {
  AxisAblation a;
  a.axis = axis;
  const RuleEngine masked = base_engine.with_axis_disabled(axis);
  for (const auto& res : full_results) {
    ++a.n_total;
    const Verdict v = masked.evaluate(res.evidence);
    if (v.confidence != res.verdict.confidence) {
      ++a.n_changed;
      const std::string key =
          std::string(to_string(res.verdict.confidence)) + "->" + to_string(v.confidence);
      ++a.transitions[key];
    }
  }
  return a;
}

}  // namespace panisoguard
