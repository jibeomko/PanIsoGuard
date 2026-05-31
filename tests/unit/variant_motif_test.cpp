#include "catch2/catch.hpp"

#include <memory>

#include "panisoguard/variant_motif.hpp"
#include "panisoguard/rules.hpp"
#include "test_util.hpp"

using namespace panisoguard;

TEST_CASE("splice_motif reads donor/acceptor and judges canonicality", "[variant]") {
  FastaFetcher ref(tiny("mini_ref.fa"));
  FastaFetcher hap(tiny("mini_hap.fa"));
  const Junction intron{20, 40};

  const SpliceMotif rm = splice_motif(ref, "chrT", intron, Strand::kPlus);
  REQUIRE(rm.evaluable);
  CHECK(rm.motif == "GG-AG");
  CHECK_FALSE(rm.canonical);

  const SpliceMotif hm = splice_motif(hap, "chrT", intron, Strand::kPlus);
  REQUIRE(hm.evaluable);
  CHECK(hm.motif == "GT-AG");
  CHECK(hm.canonical);

  // Minus strand reads donor/acceptor reverse-complemented.
  const SpliceMotif rmm = splice_motif(ref, "chrT", intron, Strand::kMinus);
  CHECK(rmm.motif == "CT-CC");  // revcomp(AG)-revcomp(GG)
  CHECK_FALSE(rmm.canonical);
}

TEST_CASE("HaplotypeProvider classifies variant-created vs disrupted", "[variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("mini_ref.fa"));
  auto hap = std::make_shared<FastaFetcher>(tiny("mini_hap.fa"));
  const Junction intron{20, 40};

  // reference non-canonical, haplotype canonical -> created (reference bias)
  HaplotypeProvider created(ref, {hap});
  CHECK(created.classify("chrT", intron, Strand::kPlus) == VariantVerdict::kCreated);

  // swap roles: "reference" canonical, "haplotype" non-canonical -> disrupted
  HaplotypeProvider disrupted(hap, {ref});
  CHECK(disrupted.classify("chrT", intron, Strand::kPlus) == VariantVerdict::kDisrupted);

  // both canonical -> none
  HaplotypeProvider none(hap, {hap});
  CHECK(none.classify("chrT", intron, Strand::kPlus) == VariantVerdict::kNone);

  // unknown contig -> not evaluable
  CHECK(created.classify("chrZ", intron, Strand::kPlus) == VariantVerdict::kNotEvaluable);
}

TEST_CASE("rules: variant rescue -> PAN_REF_RESCUED, circular-risk held", "[variant][rules]") {
  const RuleEngine eng;
  EvidenceVector ev;
  ev.structural_category = "novel_not_in_catalog";
  ev.is_novel = true;
  ev.chain_available = true;
  ev.variant_evaluable = true;
  ev.variant_rescue = true;

  ev.variant_circular = false;  // independent (wgs/external) provenance
  Verdict v = eng.evaluate(ev);
  CHECK(v.primary_mechanism == Mechanism::kVariant);
  CHECK(v.confidence == ConfidenceClass::kPanRefRescuedFalseNovel);
  CHECK_FALSE(v.circularity_flag);

  ev.variant_circular = true;   // RNA-derived/unknown -> never promoted
  Verdict vc = eng.evaluate(ev);
  CHECK(vc.confidence == ConfidenceClass::kAmbiguous);
  CHECK(vc.circularity_flag);
}
