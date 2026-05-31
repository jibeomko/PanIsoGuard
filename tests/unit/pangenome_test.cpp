#include "catch2/catch.hpp"

#include "panisoguard/pangenome.hpp"
#include "panisoguard/rules.hpp"

using namespace panisoguard;

TEST_CASE("pangenome junction reader parses + exact strand-aware lookup", "[pangenome]") {
  const std::string path = std::string(PANISOGUARD_TEST_DATA_DIR) + "/mini.pangenome.tsv";
  const PangenomeJunctions pj = read_pangenome_junctions(path);
  REQUIRE(pj.size() == 3);

  // 1-based inclusive [101,200] -> 0-based half-open [100,200].
  REQUIRE(pj.find_exact("chr1", Strand::kPlus, Junction{100, 200}) != nullptr);
  // Opposite strand at the same coordinates must NOT match.
  REQUIRE(pj.find_exact("chr1", Strand::kMinus, Junction{100, 200}) == nullptr);
  // Off-by-one is not a match (exact join, no fuzz).
  REQUIRE(pj.find_exact("chr1", Strand::kPlus, Junction{100, 201}) == nullptr);
  // Absent chromosome.
  REQUIRE(pj.find_exact("chrX", Strand::kPlus, Junction{100, 200}) == nullptr);
}

TEST_CASE("pangenome supports() honours the haplotype-count threshold", "[pangenome]") {
  const std::string path = std::string(PANISOGUARD_TEST_DATA_DIR) + "/mini.pangenome.tsv";
  const PangenomeJunctions pj = read_pangenome_junctions(path);

  // chr1 + [100,200] has n_haplotypes = 3.
  REQUIRE(pj.supports("chr1", Strand::kPlus, Junction{100, 200}, 1));
  REQUIRE(pj.supports("chr1", Strand::kPlus, Junction{100, 200}, 3));
  REQUIRE_FALSE(pj.supports("chr1", Strand::kPlus, Junction{100, 200}, 4));
  // chr1 - [300,400] has default n_haplotypes = 1.
  REQUIRE(pj.supports("chr1", Strand::kMinus, Junction{300, 400}, 1));
  REQUIRE_FALSE(pj.supports("chr1", Strand::kMinus, Junction{300, 400}, 2));
  // Absent junction is never supported.
  REQUIRE_FALSE(pj.supports("chr2", Strand::kMinus, Junction{50, 150}, 1));
}

TEST_CASE("pangenome rescue: graph-supported novel junction -> PAN_REF_RESCUED", "[pangenome]") {
  const RuleEngine engine;  // built-in defaults: use_pangenome = true
  EvidenceVector ev;
  ev.structural_category = "novel_in_catalog";
  ev.is_novel = true;
  ev.chain_available = true;
  ev.n_novel_junctions = 1;
  ev.pangenome_evaluable = true;
  ev.pangenome_rescue = true;
  ev.n_novel_jx_pangenome = 1;

  SECTION("graph support promotes to PAN_REF_RESCUED with population_known mechanism") {
    const Verdict v = engine.evaluate(ev);
    REQUIRE(v.confidence == ConfidenceClass::kPanRefRescuedFalseNovel);
    REQUIRE(v.primary_mechanism == Mechanism::kPopulationKnown);
    REQUIRE_FALSE(v.circularity_flag);  // population reference data is non-circular
  }

  SECTION("pangenome rescues even when a sample variant would be circular-risk") {
    ev.variant_evaluable = true;
    ev.variant_rescue = true;
    ev.variant_circular = true;  // sample haplotype provenance rna-derived/unknown
    const Verdict v = engine.evaluate(ev);
    REQUIRE(v.confidence == ConfidenceClass::kPanRefRescuedFalseNovel);
    REQUIRE(v.primary_mechanism == Mechanism::kPopulationKnown);
    REQUIRE_FALSE(v.circularity_flag);
  }

  SECTION("disabling the pangenome axis removes the rescue") {
    const RuleEngine off = engine.with_axis_disabled("pangenome");
    const Verdict v = off.evaluate(ev);
    REQUIRE(v.confidence != ConfidenceClass::kPanRefRescuedFalseNovel);
  }
}
