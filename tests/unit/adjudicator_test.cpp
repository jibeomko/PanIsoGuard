#include "catch2/catch.hpp"

#include <map>
#include <string>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/gtf.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
SqantiRecord rec(const std::string& id, const std::string& cat, const std::string& all_canon) {
  SqantiRecord r;
  r.isoform = id;
  r.chrom = "chr1";
  r.strand = Strand::kPlus;
  r.structural_category = cat;
  r.all_canonical = all_canon;
  return r;
}
IntronChain chain(std::vector<Junction> introns) {
  IntronChain c;
  c.chrom = "chr1";
  c.strand = Strand::kPlus;
  c.introns = std::move(introns);
  return c;
}
SjRecord sj(Junction j, int n_uniq, int motif = 1) {
  return SjRecord{"chr1", j, Strand::kPlus, motif, true, n_uniq, 0, 30};
}
const AdjudicationResult* find(const std::vector<AdjudicationResult>& v, const std::string& id) {
  for (const auto& r : v)
    if (r.isoform_id == id) return &r;
  return nullptr;
}
}  // namespace

TEST_CASE("adjudicator end-to-end: SR-supported novel vs unsupported artifact", "[adjudicator]") {
  // Catalog (mini.gtf) knows introns {200,300},{400,500} (T1) and {200,500} (T2).
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoNNC", "novel_not_in_catalog", "non_canonical"));
  sqanti.records.push_back(rec("isoNNC2", "novel_not_in_catalog", "canonical"));
  sqanti.records.push_back(rec("isoFSM", "full-splice_match", "canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoNNC"] = chain({{200, 350}, {400, 500}});   // novel jx {200,350}, no SR support
  chains["isoNNC2"] = chain({{200, 360}, {400, 500}});  // novel jx {200,360}, SR-supported below

  SjTable sjt;
  sjt.add(sj({400, 500}, 8));    // known junction (ignored by novelty check)
  sjt.add(sj({200, 360}, 10));   // supports isoNNC2's novel junction
  // note: {200,350} is intentionally absent -> isoNNC's novel junction unsupported

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.sj = &sjt;
  const auto out = adjudicate(in, RuleEngine());
  REQUIRE(out.size() == 3);

  const AdjudicationResult* nnc = find(out, "isoNNC");
  REQUIRE(nnc != nullptr);
  CHECK(nnc->evidence.n_novel_junctions == 1);
  CHECK(nnc->evidence.n_novel_jx_sr_supported == 0);
  CHECK(nnc->verdict.novelty_support == NoveltySupport::kUnsupported);
  CHECK(nnc->verdict.primary_mechanism == Mechanism::kNoncanonical);
  CHECK(nnc->verdict.confidence == ConfidenceClass::kArtifact);

  const AdjudicationResult* nnc2 = find(out, "isoNNC2");
  REQUIRE(nnc2 != nullptr);
  CHECK(nnc2->evidence.n_novel_junctions == 1);
  CHECK(nnc2->evidence.n_novel_jx_sr_supported == 1);
  CHECK(nnc2->verdict.novelty_support == NoveltySupport::kSupported);
  CHECK(nnc2->verdict.confidence == ConfidenceClass::kHighConfNovel);

  const AdjudicationResult* fsm = find(out, "isoFSM");
  REQUIRE(fsm != nullptr);
  CHECK(fsm->verdict.confidence == ConfidenceClass::kHighConfKnown);
}
