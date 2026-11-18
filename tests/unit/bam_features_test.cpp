#include "catch2/catch.hpp"

#include <cmath>
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

TEST_CASE("indel-near is counted in excess of the reads' own indel rate", "[bam]") {
  BamFeatureParams p;
  // noisy.bam: three reads with an indel every ~7 bases (6 ops over 40 aligned bases) and one
  // clean read, all spanning [2000,3000). The noisy reads each have an indel next to the
  // junction, but at their own rate the chance of one inside the two 20-bp windows is ~1, so
  // the excess over background is ~0 and the gate does not fire.
  BamReader noisy(tiny("noisy.bam"));
  const JunctionBamFeatures f = noisy.features_at_junction("chr1", Junction{2000, 3000}, Strand::kPlus, p);
  REQUIRE(f.evaluated);
  CHECK(f.n_spanning == 4);
  CHECK(f.n_indel_near == 3);
  CHECK(f.frac_indel_near() == Approx(0.75));
  CHECK(f.indel_near_background == Approx(0.75 * (1.0 - std::exp(-6.0))).margin(1e-6));   // 3 reads at rate 6/40, 1 at 0
  CHECK(f.frac_indel_near_excess() < 0.01);
  // mini.bam: clean reads with one indel in one read; the background stays small and the
  // excess is the raw fraction minus it.
  BamReader hifi(tiny("mini.bam"));
  const JunctionBamFeatures g = hifi.features_at_junction("chr1", Junction{2000, 3000}, Strand::kPlus, p);
  CHECK(g.indel_near_background < 0.1);
  CHECK(g.frac_indel_near_excess() == Approx(g.frac_indel_near() - g.indel_near_background));
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
  CHECK(s1.n_span == 4);      // the same four reach from before the intron to past it: the chain is all of the depth

  IntronChain two = one;      // a second intron no read has: the reads are consistent but too short
  two.introns.push_back({3500, 4000});
  const ChainBamSupport s2 = r.chain_support(two);
  CHECK(s2.n_full == 0);
  CHECK(s2.n_partial == 4);
  CHECK(s2.n_span == 0);      // no read reaches past the second intron

  IntronChain mono = one;
  mono.introns.clear();
  CHECK_FALSE(r.chain_support(mono).evaluated);
  IntronChain absent = one;
  absent.chrom = "chrZ";
  CHECK_FALSE(r.chain_support(absent).evaluated);
}

TEST_CASE("junction_reads returns the primary alignments that splice exactly at a junction", "[bam][phase]") {
  BamReader r(tiny("mini.bam"));
  const auto reads = r.junction_reads("chr1", {2000, 3000}, Strand::kPlus, 200);
  REQUIRE(reads.size() == 4);   // r1, r2, r4, r5 in file order; r3 is supplementary
  for (const JunctionRead& x : reads) {
    CHECK(x.start == 1900);
    CHECK(x.blocks.front().ref == 1900);
    CHECK(x.blocks.back().ref + x.blocks.back().len == x.end);
    CHECK(x.base_at(2500) == 0);   // inside the intron (the fixture's reads carry no sequence; bases are tested on indel.bam)
  }
  CHECK(reads[0].end == 3100);
  CHECK(reads[3].end == 3095);   // r5: five of its last hundred bases are an insertion
  CHECK(reads[2].blocks.front().qry == 30);   // r4: thirty soft-clipped bases precede the first aligned one
  REQUIRE(reads[3].insertions.size() == 1);   // r5: an insertion right after the intron
  CHECK(reads[3].insertions[0].first == 3000);
  CHECK(reads[3].insertions[0].second == 5);
  CHECK(r.junction_reads("chr1", {2000, 3000}, Strand::kPlus, 2).size() == 2);   // the cap
  CHECK(r.junction_reads("chr1", {2000, 3001}, Strand::kPlus, 200).empty());     // another junction
  CHECK(r.junction_reads("chrZ", {2000, 3000}, Strand::kPlus, 200).empty());     // contig absent
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
