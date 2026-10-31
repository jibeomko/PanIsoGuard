#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>

#include "panisoguard/bam_features.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("BAM features at a junction count spanning reads and artifact signatures", "[bam]") {
  BamReader r(tiny("mini.bam"));
  BamFeatureParams p;  // min_mapq=20, window=10

  // The fixture has 5 reads with an N op exactly at 0-based [2000,3000):
  // r1 clean, r2 MAPQ 5, r3 supplementary, r4 30bp soft-clip, r5 indel near 3' boundary.
  const JunctionBamFeatures f = r.features_at_junction("chr1", Junction{2000, 3000}, Strand::kPlus, p);
  REQUIRE(f.evaluated);
  CHECK(f.n_spanning == 5);
  CHECK(f.n_low_mapq == 1);
  CHECK(f.n_supplementary == 1);
  CHECK(f.n_indel_near == 1);
  CHECK(f.frac_low_mapq() == Approx(0.2));

  // A junction with no spanning reads.
  const JunctionBamFeatures none = r.features_at_junction("chr1", Junction{7000, 8000}, Strand::kPlus, p);
  CHECK(none.evaluated);
  CHECK(none.n_spanning == 0);

  // A contig absent from the header is not_evaluable.
  const JunctionBamFeatures absent = r.features_at_junction("chrZ", Junction{2000, 3000}, Strand::kPlus, p);
  CHECK_FALSE(absent.evaluated);
}

TEST_CASE("BAM reader throws on a missing file", "[bam]") {
  CHECK_THROWS(BamReader("/nonexistent/path/to.bam"));
}

TEST_CASE("BAM chain support counts primary reads carrying the whole intron chain", "[bam][chain]") {
  BamReader r(tiny("mini.bam"));
  IntronChain one;
  one.chrom = "chr1";
  one.strand = Strand::kPlus;
  one.introns = {{2000, 3000}};
  const ChainBamSupport s1 = r.chain_support(one);
  REQUIRE(s1.evaluated);
  CHECK(s1.n_full == 4);      // r1, r2 (low MAPQ still counts), r4, r5; r3 is supplementary
  CHECK(s1.n_partial == 0);

  IntronChain two = one;      // a second intron no read has: the reads are consistent but too short
  two.introns.push_back({3500, 4000});
  const ChainBamSupport s2 = r.chain_support(two);
  CHECK(s2.n_full == 0);
  CHECK(s2.n_partial == 4);

  IntronChain mono = one;
  mono.introns.clear();
  CHECK_FALSE(r.chain_support(mono).evaluated);
  IntronChain absent = one;
  absent.chrom = "chrZ";
  CHECK_FALSE(r.chain_support(absent).evaluated);
}

TEST_CASE("BAM reads on the other strand do not support a junction or a chain", "[bam][strand]") {
  // strand.sam: five reads spanning chr1 [1000,1100): two XS:A:- (STAR style), one ts:A:+ on a
  // forward alignment (= plus), one ts:A:+ on a reverse alignment (minimap2: = minus), one untagged.
  BamReader r(tiny("strand.bam"));
  BamFeatureParams p;
  const Junction jx{1000, 1100};
  CHECK(r.features_at_junction("chr1", jx, Strand::kPlus, p).n_spanning == 2);     // ts+fwd, untagged
  CHECK(r.features_at_junction("chr1", jx, Strand::kMinus, p).n_spanning == 4);    // XS- x2, ts+rev, untagged
  CHECK(r.features_at_junction("chr1", jx, Strand::kUnknown, p).n_spanning == 5);  // unknown strand: every read
  IntronChain plus;
  plus.chrom = "chr1";
  plus.strand = Strand::kPlus;
  plus.introns = {jx};
  CHECK(r.chain_support(plus).n_full == 2);
  IntronChain minus = plus;
  minus.strand = Strand::kMinus;
  CHECK(r.chain_support(minus).n_full == 4);
  IntronChain unknown = plus;
  unknown.strand = Strand::kUnknown;
  CHECK(r.chain_support(unknown).n_full == 5);
}

TEST_CASE("BAM read errors are errors, not regions with no reads", "[bam]") {
  // The first 300 bytes of mini.bam: the header block is whole, the record block is cut.
  const std::string trunc = scratch("truncated.bam");
  {
    std::ifstream in(tiny("mini.bam"), std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::ofstream(trunc, std::ios::binary) << bytes.substr(0, 300);
    std::ifstream idx(tiny("mini.bam.bai"), std::ios::binary);
    std::ofstream(trunc + ".bai", std::ios::binary) << idx.rdbuf();
  }
  BamReader r(trunc);
  BamFeatureParams p;
  CHECK_THROWS(r.features_at_junction("chr1", Junction{2000, 3000}, Strand::kPlus, p));
  IntronChain one;
  one.chrom = "chr1";
  one.strand = Strand::kPlus;
  one.introns = {{2000, 3000}};
  CHECK_THROWS(r.chain_support(one));
  std::remove(trunc.c_str());
  std::remove((trunc + ".bai").c_str());
}
