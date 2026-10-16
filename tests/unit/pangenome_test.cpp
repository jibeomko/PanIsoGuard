#include "catch2/catch.hpp"

#include <map>
#include <fstream>
#include <stdexcept>
#include <string>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/pangenome.hpp"
#include "panisoguard/rules.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
// Write a tiny pangenome-junction file to a temp path and return it.
std::string write_tmp(const std::string& name, const std::string& body) {
  const std::string p = scratch(name);
  std::ofstream(p) << body;
  return p;
}
}  // namespace

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

TEST_CASE("pangenome duplicate junction rows keep strongest haplotype support", "[pangenome]") {
  PangenomeJunctions pj;
  pj.add(PangenomeJunctionRecord{"chr1", Junction{100, 200}, Strand::kPlus, 2});
  pj.add(PangenomeJunctionRecord{"chr1", Junction{100, 200}, Strand::kPlus, 5});
  pj.add(PangenomeJunctionRecord{"chr1", Junction{100, 200}, Strand::kPlus, 3});

  REQUIRE(pj.size() == 1);
  REQUIRE(pj.supports("chr1", Strand::kPlus, Junction{100, 200}, 5));
  REQUIRE_FALSE(pj.supports("chr1", Strand::kPlus, Junction{100, 200}, 6));
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

TEST_CASE("pangenome rescue: graph-supported novel junction -> REFERENCE_BIAS", "[pangenome]") {
  const RuleEngine engine;  // built-in defaults: use_pangenome = true
  EvidenceVector ev;
  ev.structural_category = "novel_in_catalog";
  ev.is_novel = true;
  ev.chain_available = true;
  ev.n_novel_junctions = 1;
  ev.pangenome_evaluable = true;
  ev.pangenome_rescue = true;
  ev.n_novel_jx_pangenome = 1;

  SECTION("graph support promotes to REFERENCE_BIAS with population_known mechanism") {
    const Verdict v = engine.evaluate(ev);
    REQUIRE(v.confidence == ConfidenceClass::kReferenceBias);
    REQUIRE(v.primary_mechanism == Mechanism::kPopulationKnown);
    REQUIRE_FALSE(v.circularity_flag);  // population reference data is non-circular
  }

  SECTION("pangenome rescues even when a sample variant would be circular-risk") {
    ev.variant_evaluable = true;
    ev.variant_rescue = true;
    ev.variant_circular = true;  // sample haplotype provenance rna-derived/unknown
    const Verdict v = engine.evaluate(ev);
    REQUIRE(v.confidence == ConfidenceClass::kReferenceBias);
    REQUIRE(v.primary_mechanism == Mechanism::kPopulationKnown);
    REQUIRE_FALSE(v.circularity_flag);
  }

  SECTION("disabling the pangenome axis removes the rescue") {
    const RuleEngine off = engine.with_axis_disabled("pangenome");
    const Verdict v = off.evaluate(ev);
    REQUIRE(v.confidence != ConfidenceClass::kReferenceBias);
  }

  SECTION("sample-derived/unknown provenance is held under the circularity firewall") {
    ev.pangenome_circular = true;  // junction set may derive from the sample's own reads
    const Verdict v = engine.evaluate(ev);
    REQUIRE(v.confidence == ConfidenceClass::kUnresolved);
    REQUIRE(v.circularity_flag);
  }
}

TEST_CASE("pangenome reader rejects malformed rows", "[pangenome]") {
  // Unrecognized strand is rejected (not silently stored as never-matching Unknown).
  REQUIRE_THROWS_AS(read_pangenome_junctions(write_tmp("pig_pg_strand.tsv", "chr1\t101\t200\t*\t3\n")),
                    std::runtime_error);
  // Non-integer coordinate throws with context instead of a bare 'stoll'.
  REQUIRE_THROWS_AS(read_pangenome_junctions(write_tmp("pig_pg_coord.tsv", "chr1\tXX\t200\t+\n")),
                    std::runtime_error);
  // n_haplotypes must be >= 1.
  REQUIRE_THROWS_AS(read_pangenome_junctions(write_tmp("pig_pg_nhap.tsv", "chr1\t101\t200\t+\t0\n")),
                    std::runtime_error);
  // A well-formed minimal (4-column) row is accepted.
  REQUIRE(read_pangenome_junctions(write_tmp("pig_pg_ok.tsv", "chr1\t101\t200\t+\n")).size() == 1);
}

TEST_CASE("pangenome reader keeps the optional haplotype/path ids and merges duplicates", "[pangenome]") {
  const PangenomeJunctions pj = read_pangenome_junctions(write_tmp(
      "pig_pg_ids.tsv",
      "chr1\t101\t200\t+\t2\tHG002#1,HG002#2\n"
      "chr1\t101\t200\t+\t1\tHG005#1\n"      // duplicate row: ids union, count follows
      "chr1\t301\t400\t+\n"));               // no ids -> path consistency not evaluable
  const PangenomeJunctionRecord* a = pj.find_exact("chr1", Strand::kPlus, Junction{100, 200});
  REQUIRE(a != nullptr);
  CHECK(a->haplotype_ids == std::vector<std::string>{"HG002#1", "HG002#2", "HG005#1"});
  CHECK(a->n_haplotypes == 3);
  const PangenomeJunctionRecord* b = pj.find_exact("chr1", Strand::kPlus, Junction{300, 400});
  REQUIRE(b != nullptr);
  CHECK(b->haplotype_ids.empty());
}

TEST_CASE("pangenome rescue requires one path carrying every novel junction", "[pangenome]") {
  // Both novel junctions are graph-supported, but on disjoint paths -> held; a shared
  // path -> rescued; rows without path ids -> rescued on counts alone (as before 0.0.6).
  const Catalog cat = build_catalog_from_gtf(std::string(PANISOGUARD_TEST_DATA_DIR) + "/mini.gtf");
  SqantiTable sqanti;
  SqantiRecord r;
  r.isoform = "iso";
  r.chrom = "chr1";
  r.strand = Strand::kPlus;
  r.structural_category = "novel_not_in_catalog";
  sqanti.records.push_back(r);
  std::map<std::string, IntronChain> chains;
  IntronChain c;
  c.chrom = "chr1";
  c.strand = Strand::kPlus;
  c.introns = {{200, 350}, {400, 550}};  // neither in mini.gtf
  chains["iso"] = c;
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome_circular = false;
  auto run = [&](const std::string& body) {
    PangenomeJunctions pan = read_pangenome_junctions(write_tmp("pig_pg_path.tsv", body));
    in.pangenome = &pan;
    return adjudicate(in, RuleEngine()).front();
  };
  const auto disjoint = run("chr1\t201\t350\t+\t1\tP1\nchr1\t401\t550\t+\t1\tP2\n");
  CHECK(disjoint.evidence.n_novel_jx_pangenome == 2);
  CHECK(disjoint.evidence.pangenome_path_evaluable);
  CHECK_FALSE(disjoint.evidence.pangenome_path_consistent);
  CHECK(disjoint.verdict.confidence == ConfidenceClass::kUnresolved);
  CHECK(disjoint.verdict.verdict_reason == "path_inconsistent");

  const auto shared = run("chr1\t201\t350\t+\t2\tP1,P2\nchr1\t401\t550\t+\t1\tP2\n");
  CHECK(shared.evidence.pangenome_path_consistent);
  CHECK(shared.verdict.confidence == ConfidenceClass::kReferenceBias);
  CHECK(shared.verdict.verdict_reason == "pangenome_path");

  const auto no_ids = run("chr1\t201\t350\t+\t1\nchr1\t401\t550\t+\t1\n");
  CHECK_FALSE(no_ids.evidence.pangenome_path_evaluable);
  CHECK(no_ids.verdict.confidence == ConfidenceClass::kReferenceBias);
}
