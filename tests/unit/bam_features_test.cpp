#include "catch2/catch.hpp"

#include "panisoguard/bam_features.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("BAM features at a junction count spanning reads and artifact signatures", "[bam]") {
  BamReader r(tiny("mini.bam"));
  BamFeatureParams p;  // min_mapq=20, softclip_min_bp=20, window=10

  // The fixture has 5 reads with an N op exactly at 0-based [2000,3000):
  // r1 clean, r2 MAPQ 5, r3 supplementary, r4 30bp soft-clip, r5 indel near 3' boundary.
  const JunctionBamFeatures f = r.features_at_junction("chr1", Junction{2000, 3000}, p);
  REQUIRE(f.evaluated);
  CHECK(f.n_spanning == 5);
  CHECK(f.n_low_mapq == 1);
  CHECK(f.n_supplementary == 1);
  CHECK(f.n_softclip == 1);
  CHECK(f.n_indel_near == 1);
  CHECK(f.frac_low_mapq() == Approx(0.2));

  // A junction with no spanning reads.
  const JunctionBamFeatures none = r.features_at_junction("chr1", Junction{7000, 8000}, p);
  CHECK(none.evaluated);
  CHECK(none.n_spanning == 0);

  // A contig absent from the header is not_evaluable.
  const JunctionBamFeatures absent = r.features_at_junction("chrZ", Junction{2000, 3000}, p);
  CHECK_FALSE(absent.evaluated);
}

TEST_CASE("BAM reader throws on a missing file", "[bam]") {
  CHECK_THROWS(BamReader("/nonexistent/path/to.bam"));
}
