#include "catch2/catch.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "panisoguard/sj_tab.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("SJ.tab reader reports a non-integer field with line context", "[sj_tab]") {
  const auto p = std::filesystem::temp_directory_path() / "pig_sj_bad.SJ.tab";
  std::ofstream(p) << "chr1\tNOPE\t300\t1\t1\t1\t10\t0\t30\n";
  REQUIRE_THROWS_AS(read_sj_tab(p.string()), std::runtime_error);
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
