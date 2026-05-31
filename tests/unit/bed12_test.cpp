#include "catch2/catch.hpp"

#include "panisoguard/bed12.hpp"
#include "panisoguard/fingerprint.hpp"
#include "panisoguard/gtf.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("BED12 decodes intron chain and FLAIR sample prefix", "[bed12]") {
  const auto recs = read_bed12(tiny("mini.bed"));
  REQUIRE(recs.size() == 1);

  const Bed12Record& r = recs[0];
  CHECK(r.id_prefix == "1-1");  // FLAIR running id, not a biological sample
  CHECK(r.chain.strand == Strand::kPlus);
  CHECK(r.chain.chrom == "chr1");
  REQUIRE(r.chain.introns == std::vector<Junction>{{200, 300}, {400, 500}});
}

TEST_CASE("BED12 isoform fingerprint matches the catalog chain (cross-reader)", "[bed12][fingerprint]") {
  // The BED isoform is structurally T1; its fingerprint must be a known chain in
  // a catalog built independently from the GTF. This is the core cross-reader
  // coordinate-convergence guarantee.
  const auto recs = read_bed12(tiny("mini.bed"));
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  REQUIRE(cat.has_chain(fingerprint_intron_chain(recs[0].chain)));
}
