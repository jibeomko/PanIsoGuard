#include "catch2/catch.hpp"

#include <cstdio>
#include <map>
#include <utility>

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

TEST_CASE("read_caller_support round-trips the combine matrix to native_id -> n_callers",
          "[consensus]") {
  // chainA: 2 callers (flair, isoquant); chainB: 1 caller; flair has two native IDs on
  // chainA's caller to exercise the '|' id-join split.
  const auto chainA = chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});
  const auto chainB = chain("chr1", Strand::kPlus, {{200, 500}});
  ConsensusBuilder b;
  b.add("flair", "1-1_x", chainA);
  b.add("flair", "1-2_x", chainA);  // second flair id on the same chain
  b.add("isoquant", "transcript.9", chainA);
  b.add("bambu", "n1", chainB);

  const std::string path = scratch("caller_support_roundtrip.tsv");
  write_caller_support_matrix(path, b.build());

  const CallerSupport support = read_caller_support(path);
  // Every native running id on chainA maps to n_callers=2; chainB's id maps to 1.
  CHECK(support.by_native_id.at("1-1_x") == 2);
  CHECK(support.by_native_id.at("1-2_x") == 2);
  CHECK(support.by_native_id.at("transcript.9") == 2);
  CHECK(support.by_native_id.at("n1") == 1);
  CHECK(support.by_native_id.count("absent") == 0);
  // ...and the chain fingerprint (the join key adjudicate uses) carries the same counts.
  REQUIRE(support.has_fingerprints());
  CHECK(support.by_chain.at(fingerprint_intron_chain(chainA)) == 2);
  CHECK(support.by_chain.at(fingerprint_intron_chain(chainB)) == 1);
  CHECK(support.n_native_id_collisions == 0);
  std::remove(path.c_str());
}

TEST_CASE("consensus: duplicate native ids across callers stay separate per chain", "[consensus]") {
  // Two callers both name an isoform "transcript_1" but mean different chains. A
  // native-id join would hand chainX's 2-caller support to caller B's transcript_1
  // (chainY, 1 caller); the fingerprint join does not.
  const auto chainX = chain("chr1", Strand::kPlus, {{200, 300}, {400, 500}});
  const auto chainY = chain("chr1", Strand::kPlus, {{200, 500}});
  ConsensusBuilder b;
  b.add("callerA", "transcript_1", chainX);
  b.add("callerB", "transcript_2", chainX);
  b.add("callerB", "transcript_1", chainY);

  const std::string path = scratch("caller_support_collision.tsv");
  write_caller_support_matrix(path, b.build());
  const CallerSupport s = read_caller_support(path);
  CHECK(s.n_native_id_collisions == 1);                 // reported, not silently max'ed
  CHECK(s.by_chain.at(fingerprint_intron_chain(chainX)) == 2);
  CHECK(s.by_chain.at(fingerprint_intron_chain(chainY)) == 1);
  std::remove(path.c_str());
}

TEST_CASE("consensus: monoexonic PIG ids do not depend on input order", "[consensus]") {
  const auto mono = chain("chr1", Strand::kPlus, {});
  auto ids = [&](std::vector<std::pair<const char*, const char*>> order) {
    ConsensusBuilder b;
    for (const auto& p : order) b.add(p.first, p.second, mono);
    std::map<std::string, std::string> id_of;  // native id -> pig id
    for (const auto& iso : b.build()) id_of[iso.support.begin()->second.front()] = iso.pig_id;
    return id_of;
  };
  const auto forward = ids({{"flair", "m1"}, {"flair", "m2"}, {"bambu", "m3"}});
  const auto reversed = ids({{"bambu", "m3"}, {"flair", "m2"}, {"flair", "m1"}});
  REQUIRE(forward.size() == 3);
  CHECK(forward == reversed);
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
