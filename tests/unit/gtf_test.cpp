#include "catch2/catch.hpp"

#include <algorithm>

#include "panisoguard/gtf.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
const Transcript* find_tx(const std::vector<Transcript>& v, const std::string& id) {
  auto it = std::find_if(v.begin(), v.end(), [&](const Transcript& t) { return t.id == id; });
  return it == v.end() ? nullptr : &*it;
}
}  // namespace

TEST_CASE("GTF attribute extraction", "[gtf]") {
  const std::string attrs = "gene_id \"G1\"; transcript_id \"T1\"; exon_number 2;";
  CHECK(extract_gtf_attr(attrs, "transcript_id") == "T1");
  CHECK(extract_gtf_attr(attrs, "gene_id") == "G1");
  CHECK(extract_gtf_attr(attrs, "missing") == "");
}

TEST_CASE("GTF transcripts derive correct intron chains", "[gtf]") {
  const auto tx = read_gtf_transcripts(tiny("mini.gtf"));
  REQUIRE(tx.size() == 2);

  const Transcript* t1 = find_tx(tx, "T1");
  REQUIRE(t1 != nullptr);
  REQUIRE(t1->chain.introns == std::vector<Junction>{{200, 300}, {400, 500}});

  const Transcript* t2 = find_tx(tx, "T2");  // exon-skip: single intron
  REQUIRE(t2 != nullptr);
  REQUIRE(t2->chain.introns == std::vector<Junction>{{200, 500}});
}

TEST_CASE("Catalog records known chains and introns", "[gtf][catalog]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  for (const auto& t : read_gtf_transcripts(tiny("mini.gtf"))) CHECK(cat.has_chain(fingerprint_intron_chain(t.chain)));

  CHECK(cat.has_intron("chr1", Strand::kPlus, Junction{200, 300}));
  CHECK(cat.has_intron("chr1", Strand::kPlus, Junction{200, 500}));
  CHECK_FALSE(cat.has_intron("chr1", Strand::kPlus, Junction{200, 400}));  // novel combination
  CHECK_FALSE(cat.has_intron("chr1", Strand::kMinus, Junction{200, 300}));  // strand-aware
}

TEST_CASE("A catalog without transcripts answers chains and introns only", "[gtf][catalog]") {
  const Catalog slim = build_catalog_from_gtf(tiny("mini.gtf"), /*with_transcripts=*/false);
  const Catalog full = build_catalog_from_gtf(tiny("mini.gtf"));
  CHECK_FALSE(slim.has_transcripts());
  CHECK(full.has_transcripts());
  for (const auto& t : read_gtf_transcripts(tiny("mini.gtf"))) CHECK(slim.has_chain(fingerprint_intron_chain(t.chain)));
  CHECK(slim.has_intron("chr1", Strand::kPlus, Junction{200, 300}));
  CHECK_FALSE(slim.has_intron("chr1", Strand::kPlus, Junction{200, 400}));
  CHECK_FALSE(slim.has_intron_start("chr1", Strand::kPlus, 200));   // the site sets belong to the transcript index
  CHECK(full.has_intron_start("chr1", Strand::kPlus, 200));
}

TEST_CASE("GTF attribute extraction reads an unquoted value", "[gtf]") {
  CHECK(extract_gtf_attr("gene_id \"G1\"; exon_number 2; transcript_id \"T1\";", "exon_number") == "2");
  CHECK(extract_gtf_attr("gene_id \"G1\"; exon_number 7", "exon_number") == "7");
}
