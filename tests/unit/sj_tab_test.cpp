#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "panisoguard/sj_tab.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("SJ.tab reader reports a non-integer field with line context", "[sj_tab]") {
  const std::string p = scratch("sj_bad.SJ.tab");
  std::ofstream(p) << "chr1\tNOPE\t300\t1\t1\t1\t10\t0\t30\n";
  REQUIRE_THROWS_AS(read_sj_tab(p), std::runtime_error);
  std::remove(p.c_str());
}

TEST_CASE("SJ.tab parses and supports exact coordinate lookup", "[sj_tab]") {
  const SjTable t = read_sj_tab(tiny("mini.SJ.tab"));
  REQUIRE(t.size() == 3);

  // 1-based 201..300 normalizes to 0-based half-open {200, 300} on + strand.
  const SjRecord* r1 = t.find_exact("chr1", Strand::kPlus, Junction{200, 300});
  REQUIRE(r1 != nullptr);
  CHECK(r1->n_uniq == 10);
  CHECK(r1->canonical());
  CHECK(r1->strand == Strand::kPlus);

  const SjRecord* r2 = t.find_exact("chr1", Strand::kPlus, Junction{400, 500});
  REQUIRE(r2 != nullptr);
  CHECK(r2->n_uniq == 8);

  // A noncanonical junk junction (minus strand in the fixture) is present but
  // flagged noncanonical.
  const SjRecord* junk = t.find_exact("chr1", Strand::kMinus, Junction{200, 450});
  REQUIRE(junk != nullptr);
  CHECK_FALSE(junk->canonical());
  CHECK(junk->n_uniq == 2);

  // A junction absent from SJ.tab returns nullptr.
  CHECK(t.find_exact("chr1", Strand::kPlus, Junction{200, 500}) == nullptr);
  CHECK(t.find_exact("chr2", Strand::kPlus, Junction{200, 300}) == nullptr);
}

TEST_CASE("SJ.tab duplicate junction rows are aggregated", "[sj_tab]") {
  SjTable t;
  t.add(SjRecord{"chr1", Junction{200, 360}, Strand::kPlus, 0, false, 2, 1, 18});
  t.add(SjRecord{"chr1", Junction{200, 360}, Strand::kPlus, 1, true, 4, 3, 31});

  REQUIRE(t.size() == 1);
  const SjRecord* hit = t.find_exact("chr1", Strand::kPlus, Junction{200, 360});
  REQUIRE(hit != nullptr);
  CHECK(hit->canonical());
  CHECK(hit->annotated);
  CHECK(hit->n_uniq == 6);
  CHECK(hit->n_multi == 4);
  CHECK(hit->max_overhang == 31);
}

TEST_CASE("SJ lookup is strand-aware (no cross-strand corroboration)", "[sj_tab]") {
  SjTable t;
  t.add(SjRecord{"chr1", Junction{200, 360}, Strand::kPlus, 1, true, 10, 0, 30});
  t.add(SjRecord{"chr1", Junction{700, 800}, Strand::kMinus, 1, true, 7, 0, 25});

  // same coordinates, opposite strand -> not found
  CHECK(t.find_exact("chr1", Strand::kPlus, Junction{200, 360}) != nullptr);
  CHECK(t.find_exact("chr1", Strand::kMinus, Junction{200, 360}) == nullptr);
  CHECK(t.find_exact("chr1", Strand::kMinus, Junction{700, 800}) != nullptr);
  CHECK(t.find_exact("chr1", Strand::kPlus, Junction{700, 800}) == nullptr);
}

TEST_CASE("SjTable: any_in_span sees a spliced read anywhere in a span, on the strand or unstranded", "[sj_tab]") {
  SjTable t;
  t.add(SjRecord{"chr1", {100, 200}, Strand::kPlus, 1, false, 3, 0, 20});
  t.add(SjRecord{"chr1", {1000, 5000}, Strand::kMinus, 1, false, 3, 0, 20});     // long intron on the other strand
  t.add(SjRecord{"chr1", {7000, 7100}, Strand::kUnknown, 0, false, 1, 0, 20});   // STAR strand 0
  CHECK(t.any_in_span("chr1", Strand::kPlus, 150, 160));
  CHECK(t.any_in_span("chr1", Strand::kPlus, 0, 101));
  CHECK_FALSE(t.any_in_span("chr1", Strand::kPlus, 200, 300));        // half-open: the intron ends at 200
  CHECK_FALSE(t.any_in_span("chr1", Strand::kPlus, 2000, 2100));      // only the minus-strand intron spans it
  CHECK(t.any_in_span("chr1", Strand::kMinus, 2000, 2100));
  CHECK(t.any_in_span("chr1", Strand::kUnknown, 2000, 2100));
  CHECK(t.any_in_span("chr1", Strand::kMinus, 7050, 7060));           // an unstranded junction counts for either strand
  CHECK_FALSE(t.any_in_span("chr2", Strand::kPlus, 150, 160));
  t.add(SjRecord{"chr2", {150, 160}, Strand::kPlus, 1, false, 3, 0, 20});   // the index is rebuilt after a later add
  CHECK(t.any_in_span("chr2", Strand::kPlus, 100, 200));
}
