#include "catch2/catch.hpp"

#include "panisoguard/decoys.hpp"

using namespace panisoguard;

namespace {
Transcript tx(const char* id, std::initializer_list<std::pair<int64_t, int64_t>> introns) {
  Transcript t;
  t.id = id;
  t.gene_id = "G1";
  t.chain.chrom = "chr1";
  t.chain.strand = Strand::kPlus;
  for (const auto& p : introns) t.chain.introns.push_back(Junction{p.first, p.second});
  return t;
}
}  // namespace

TEST_CASE("decoys: prefix + suffix at a shared intron, never an annotated chain or a piece of one", "[decoys]") {
  const Transcript A = tx("A", {{100, 200}, {300, 400}, {500, 600}});
  const Transcript B = tx("B", {{120, 220}, {300, 400}, {520, 620}});
  const Transcript C = tx("C", {{100, 200}, {300, 400}});          // a 5' piece of A
  std::vector<Decoy> d = make_decoys({A, B, C}, 0, 1);
  REQUIRE(d.size() == 2);                                           // A|B and B|A; C adds nothing new
  for (const auto& x : d) {
    CHECK(x.gene_id == "G1");
    CHECK(x.chain.introns.size() == 3);
    CHECK(x.chain.introns[1] == Junction{300, 400});                // the shared intron
    CHECK(x.chain.introns != A.chain.introns);
    CHECK(x.chain.introns != B.chain.introns);
  }
  CHECK(d[0].chain.introns.front() != d[1].chain.introns.front());  // one starts like A, one like B
  // the cap samples deterministically
  CHECK(make_decoys({A, B, C}, 1, 1).size() == 1);
  CHECK(make_decoys({A, B, C}, 1, 1)[0].chain.introns == make_decoys({A, B, C}, 1, 1)[0].chain.introns);
  // a different gene never recombines with this one
  Transcript D = tx("D", {{130, 230}, {300, 400}, {540, 640}});
  D.gene_id = "G2";
  CHECK(make_decoys({A, D}, 0, 1).empty());
  // monoexonic / single-intron transcripts and BED-derived (no gene_id) ones are skipped
  Transcript E = tx("E", {{300, 400}});
  CHECK(make_decoys({A, E}, 0, 1).empty());
}
