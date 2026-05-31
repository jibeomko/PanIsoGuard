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

TEST_CASE("adjudicator: pangenome graph-supported novel junction -> PAN_REF_RESCUED", "[adjudicator]") {
  // Catalog knows {200,300},{400,500} (T1) and {200,500} (T2); {200,350} is novel.
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPan", "novel_not_in_catalog", "non_canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoPan"] = chain({{200, 350}, {400, 500}});  // novel jx {200,350}

  // The pangenome carries that novel junction on a graph haplotype path.
  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  const auto out = adjudicate(in, RuleEngine());
  REQUIRE(out.size() == 1);

  const AdjudicationResult* p = find(out, "isoPan");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.n_novel_junctions == 1);
  CHECK(p->evidence.pangenome_evaluable);
  CHECK(p->evidence.pangenome_rescue);
  CHECK(p->evidence.n_novel_jx_pangenome == 1);
  CHECK(p->verdict.confidence == ConfidenceClass::kPanRefRescuedFalseNovel);
  CHECK(p->verdict.primary_mechanism == Mechanism::kPopulationKnown);
  CHECK_FALSE(p->verdict.circularity_flag);  // population reference data is non-circular

  // Disabling the pangenome axis must drop the rescue (no SR/variant evidence here).
  const auto out_off = adjudicate(in, RuleEngine().with_axis_disabled("pangenome"));
  CHECK(find(out_off, "isoPan")->verdict.confidence != ConfidenceClass::kPanRefRescuedFalseNovel);
}

TEST_CASE("adjudicator: pangenome partial support does NOT rescue (all-or-nothing)", "[adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPartial", "novel_not_in_catalog", "canonical"));

  std::map<std::string, IntronChain> chains;
  // BOTH introns are novel vs the catalog ({200,350} and {600,900}).
  chains["isoPartial"] = chain({{200, 350}, {600, 900}});

  // The pangenome carries only ONE of the two novel junctions.
  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  const auto out = adjudicate(in, RuleEngine());

  const AdjudicationResult* p = find(out, "isoPartial");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.n_novel_junctions == 2);
  CHECK(p->evidence.n_novel_jx_pangenome == 1);
  CHECK_FALSE(p->evidence.pangenome_rescue);  // 1 of 2 -> the isoform is not fully reference bias
  CHECK(p->verdict.confidence != ConfidenceClass::kPanRefRescuedFalseNovel);
}

TEST_CASE("adjudicator: sample-derived pangenome provenance is held (circularity firewall)", "[adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPan", "novel_not_in_catalog", "canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoPan"] = chain({{200, 350}, {400, 500}});  // single novel jx {200,350}

  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  in.pangenome_circular = true;  // junction set may derive from the sample's own reads
  const auto out = adjudicate(in, RuleEngine());

  const AdjudicationResult* p = find(out, "isoPan");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.pangenome_rescue);  // all novel junctions are graph-supported...
  CHECK(p->evidence.pangenome_circular);
  CHECK(p->verdict.confidence == ConfidenceClass::kAmbiguous);  // ...but held, not promoted
  CHECK(p->verdict.circularity_flag);
}
