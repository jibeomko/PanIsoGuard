#include "catch2/catch.hpp"

#include <string>

#include "panisoguard/analysis.hpp"

using namespace panisoguard;

namespace {
AdjudicationResult mk(const std::string& id, ConfidenceClass cls, Mechanism mech = Mechanism::kNone) {
  AdjudicationResult r;
  r.isoform_id = id;
  r.structural_category = "novel_not_in_catalog";
  r.verdict.confidence = cls;
  r.verdict.primary_mechanism = mech;
  return r;
}
SqantiRecord sq(const std::string& id, const std::string& filter) {
  SqantiRecord r;
  r.isoform = id;
  r.structural_category = "novel_not_in_catalog";
  r.filter_result = filter;
  return r;
}
}  // namespace

TEST_CASE("non-redundancy 2x2 vs SQANTI3", "[analysis]") {
  std::vector<AdjudicationResult> results = {
      mk("iso1", ConfidenceClass::kArtifact, Mechanism::kMapping),   // PIG flag
      mk("iso2", ConfidenceClass::kArtifact),                        // PIG flag
      mk("iso3", ConfidenceClass::kHighConfNovel),                   // PIG keep
      mk("iso4", ConfidenceClass::kHighConfNovel),                   // PIG keep
  };
  SqantiTable s;
  s.records = {sq("iso1", "Isoform"),   // SQANTI keep  -> pig_only_flag
               sq("iso2", "Artifact"),  // SQANTI flag  -> both_flag
               sq("iso3", "Artifact"),  // SQANTI flag  -> sqanti_only_flag
               sq("iso4", "Isoform")};  // SQANTI keep  -> both_keep

  const NonRedundancy nr = compute_nonredundancy(results, s);
  CHECK(nr.n_novel == 4);
  CHECK(nr.n_sqanti_filter_available == 4);
  CHECK(nr.both_flag == 1);
  CHECK(nr.both_keep == 1);
  CHECK(nr.pig_only_flag == 1);
  CHECK(nr.sqanti_only_flag == 1);
  CHECK(nr.jaccard == Approx(1.0 / 3.0));            // 1 / (1+1+1)
  CHECK(nr.mcnemar_chi2 == Approx(0.0));             // b==c==1 with continuity correction
  REQUIRE(nr.pig_only_examples.size() == 1);
  CHECK(nr.pig_only_examples[0].first == "iso1");
  CHECK(nr.pig_only_examples[0].second == "mapping_or_repeat");
}

TEST_CASE("ablation counts class changes when an axis is masked", "[analysis]") {
  const RuleEngine eng;

  // A novel isoform that is HIGH_CONF_NOVEL only because short-read support is on.
  EvidenceVector ev;
  ev.structural_category = "novel_not_in_catalog";
  ev.is_novel = true;
  ev.chain_available = true;
  ev.sj_evaluable = true;
  ev.n_novel_junctions = 1;
  ev.n_novel_jx_sr_supported = 1;

  AdjudicationResult r;
  r.isoform_id = "isoX";
  r.structural_category = ev.structural_category;
  r.evidence = ev;
  r.verdict = eng.evaluate(ev);
  REQUIRE(r.verdict.confidence == ConfidenceClass::kHighConfNovel);

  const AxisAblation ab = compute_ablation({r}, eng, "short_read");
  CHECK(ab.n_total == 1);
  CHECK(ab.n_changed == 1);  // masking short_read -> support UNKNOWN -> AMBIGUOUS
  CHECK(ab.transitions.at("HIGH_CONF_NOVEL->AMBIGUOUS") == 1);

  // Masking an axis the isoform does not use changes nothing.
  const AxisAblation ab2 = compute_ablation({r}, eng, "degradation");
  CHECK(ab2.n_changed == 0);
}
