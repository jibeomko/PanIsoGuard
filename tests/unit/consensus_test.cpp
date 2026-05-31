#include "catch2/catch.hpp"

#include "panisoguard/consensus.hpp"
#include "panisoguard/gtf.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
IntronChain chain(const char* chrom, Strand s, std::vector<Junction> introns) {
  IntronChain c;
  c.chrom = chrom;
  c.strand = s;
  c.introns = std::move(introns);
  return c;
}
const ConsensusIsoform* by_introns(const std::vector<ConsensusIsoform>& v,
                                   const std::vector<Junction>& introns) {
  for (const auto& iso : v)
    if (iso.chain.introns == introns) return &iso;
  return nullptr;
}
}  // namespace

TEST_CASE("consensus merges callers by intron chain, preserving native IDs", "[consensus]") {
  // chainA (= GTF T1) called by flair and isoquant under different running IDs;
  // chainB called only by flair.
  const auto chainA = chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});
  const auto chainB = chain("chr1", Strand::kPlus, {{200, 500}});  // = GTF T2

  ConsensusBuilder b;
  b.add("flair", "1-1_x", chainA);
  b.add("isoquant", "transcript.9", chainA);  // same chain, different running ID
  b.add("flair", "2-1_y", chainB);

  const auto iso = b.build();
  REQUIRE(iso.size() == 2);

  const ConsensusIsoform* a = by_introns(iso, chainA.introns);
  REQUIRE(a != nullptr);
  CHECK(a->n_callers() == 2);
  CHECK(a->n_isoforms() == 2);
  REQUIRE(a->support.at("flair") == std::vector<std::string>{"1-1_x"});
  REQUIRE(a->support.at("isoquant") == std::vector<std::string>{"transcript.9"});

  const ConsensusIsoform* bb = by_introns(iso, chainB.introns);
  REQUIRE(bb != nullptr);
  CHECK(bb->n_callers() == 1);

  // PIG ids are assigned and unique.
  CHECK(iso[0].pig_id != iso[1].pig_id);
  CHECK(iso[0].pig_id.rfind("PIG.", 0) == 0);
}

TEST_CASE("consensus annotates known/novel against a catalog", "[consensus][catalog]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  const auto known = chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});  // T1: known
  const auto novel = chain("chr1", Strand::kPlus, {{200, 350}, {400, 500}});  // novel site

  ConsensusBuilder b;
  b.add("flair", "k1", known);
  b.add("bambu", "n1", novel);
  const auto iso = b.build(&cat);

  const ConsensusIsoform* k = by_introns(iso, known.introns);
  const ConsensusIsoform* n = by_introns(iso, novel.introns);
  REQUIRE((k != nullptr && n != nullptr));
  CHECK(k->catalog_checked);
  CHECK(k->known_in_catalog);
  CHECK(n->catalog_checked);
  CHECK_FALSE(n->known_in_catalog);
}
