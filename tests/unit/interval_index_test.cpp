#include "catch2/catch.hpp"

#include <algorithm>
#include <stdexcept>

#include "panisoguard/interval_index.hpp"

using namespace panisoguard;

namespace {
std::vector<uint32_t> sorted(std::vector<uint32_t> v) {
  std::sort(v.begin(), v.end());
  return v;
}
}  // namespace

TEST_CASE("IntervalIndex returns overlapping payloads after finalize", "[interval]") {
  IntervalIndex idx;
  idx.add(10, 20, 1);
  idx.add(15, 25, 2);
  idx.add(100, 110, 3);
  idx.finalize();

  REQUIRE(idx.is_indexed());
  REQUIRE(idx.size() == 3);
  REQUIRE(sorted(idx.overlap(18, 19)) == std::vector<uint32_t>{1, 2});
  REQUIRE(sorted(idx.overlap(12, 13)) == std::vector<uint32_t>{1});
  REQUIRE(idx.overlap(50, 60).empty());
  REQUIRE(sorted(idx.overlap(105, 106)) == std::vector<uint32_t>{3});
}

TEST_CASE("IntervalIndex enforces the build/query phase split", "[interval]") {
  IntervalIndex idx;
  idx.add(0, 5, 7);
  REQUIRE_THROWS_AS(idx.overlap(1, 2), std::logic_error);  // query before finalize

  idx.finalize();
  REQUIRE_THROWS_AS(idx.add(6, 7, 8), std::logic_error);   // add after finalize
}
