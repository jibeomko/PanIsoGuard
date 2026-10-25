#include "catch2/catch.hpp"

#include "panisoguard/sqanti.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("SQANTI3 classification parses with name-indexed columns", "[sqanti]") {
  const SqantiTable t = read_sqanti_classification(tiny("mini_classification.txt"));

  REQUIRE(t.records.size() == 4);
  CHECK(t.n_columns == 14);
  CHECK(t.records[0].structural_category == "full-splice_match");
  CHECK(t.records[3].structural_category == "incomplete-splice_match");
}

TEST_CASE("SQANTI3 novelty and QC priors", "[sqanti]") {
  const SqantiTable t = read_sqanti_classification(tiny("mini_classification.txt"));

  std::size_t n_novel = 0;
  for (const auto& r : t.records)
    if (r.is_novel()) ++n_novel;
  CHECK(n_novel == 2);

  const SqantiRecord& iso1 = t.records[0];
  CHECK(iso1.isoform == "iso1");
  CHECK(iso1.perc_A_downstream_TTS == Approx(12.5));
  CHECK_FALSE(sqanti_is_na(iso1.perc_A_downstream_TTS));

  const SqantiRecord& iso2 = t.records[1];  // NIC
  CHECK(iso2.is_novel());

  const SqantiRecord& iso3 = t.records[2];  // NNC artifact
  CHECK(iso3.structural_category == "novel_not_in_catalog");
  CHECK(iso3.all_canonical == "non_canonical");

  CHECK(t.records[3].strand == Strand::kMinus);
}
