#include "catch2/catch.hpp"

#include "panisoguard/fingerprint.hpp"

using namespace panisoguard;

namespace {
IntronChain make_chain(const char* chrom, Strand s, std::vector<Junction> introns) {
  IntronChain c;
  c.chrom = chrom;
  c.strand = s;
  c.introns = std::move(introns);
  return c;
}
}  // namespace

TEST_CASE("fingerprint is order-independent over introns", "[fingerprint]") {
  const auto a = make_chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});
  const auto b = make_chain("chr1", Strand::kPlus, {{400, 500}, {200, 300}});
  REQUIRE(fingerprint_intron_chain(a) == fingerprint_intron_chain(b));
}

TEST_CASE("fingerprint distinguishes a differing junction", "[fingerprint]") {
  const auto a = make_chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});
  const auto c = make_chain("chr1", Strand::kPlus, {{200, 300}, {400, 501}});
  REQUIRE(fingerprint_intron_chain(a) != fingerprint_intron_chain(c));
}

TEST_CASE("fingerprint depends on chrom and strand", "[fingerprint]") {
  const auto a = make_chain("chr1", Strand::kPlus, {{200, 300}});
  REQUIRE(fingerprint_intron_chain(a) !=
          fingerprint_intron_chain(make_chain("chr2", Strand::kPlus, {{200, 300}})));
  REQUIRE(fingerprint_intron_chain(a) !=
          fingerprint_intron_chain(make_chain("chr1", Strand::kMinus, {{200, 300}})));
}
