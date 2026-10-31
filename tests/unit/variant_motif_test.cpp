#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>
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
  CHECK(created.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kCreated);

  // swap roles: "reference" canonical, "haplotype" non-canonical -> disrupted
  HaplotypeProvider disrupted(hap, {ref});
  CHECK(disrupted.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kDisrupted);

  // both canonical -> none
  HaplotypeProvider none(hap, {hap});
  CHECK(none.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kNone);

  // unknown contig -> not evaluable
  CHECK(created.motifs("chrZ", intron, Strand::kPlus).verdict == VariantVerdict::kNotEvaluable);
}

TEST_CASE("rules: variant rescue -> REFERENCE_BIAS, circular-risk held", "[variant][rules]") {
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
  CHECK(v.confidence == ConfidenceClass::kReferenceBias);
  CHECK_FALSE(v.circularity_flag);

  ev.variant_circular = true;   // RNA-derived/unknown -> never promoted
  Verdict vc = eng.evaluate(ev);
  CHECK(vc.confidence == ConfidenceClass::kUnresolved);
  CHECK(vc.circularity_flag);
}

namespace {
// A one-contig FASTA on disk (faidx builds the .fai); removed with the returned guard.
struct TempFasta {
  std::string path;
  TempFasta(const std::string& name, const std::string& contig, const std::string& seq) : path(scratch(name)) {
    std::ofstream(path) << ">" << contig << "\n" << seq << "\n";
  }
  ~TempFasta() {
    std::remove(path.c_str());
    std::remove((path + ".fai").c_str());
  }
};
}  // namespace

TEST_CASE("an N at the splice site is not evaluable, so it cannot create or disrupt a motif", "[variant]") {
  // intron {20,40}: donor = bases 20-21, acceptor = bases 38-39
  const std::string a20(20, 'A'), c16(16, 'C');
  TempFasta n_ref("n_ref.fa", "chrT", a20 + "NT" + c16 + "AG" + a20);
  TempFasta gt_hap("gt_hap.fa", "chrT", a20 + "GT" + c16 + "AG" + a20);
  const Junction intron{20, 40};
  CHECK_FALSE(splice_motif(FastaFetcher(n_ref.path), "chrT", intron, Strand::kPlus).evaluable);
  HaplotypeProvider hp(std::make_shared<FastaFetcher>(n_ref.path), {std::make_shared<FastaFetcher>(gt_hap.path)});
  CHECK(hp.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kNotEvaluable);
  // an unknown strand cannot orient donor and acceptor
  CHECK_FALSE(splice_motif(FastaFetcher(gt_hap.path), "chrT", intron, Strand::kUnknown).evaluable);
}

TEST_CASE("'disrupted on every haplotype' needs every haplotype readable", "[variant]") {
  auto canonical = std::make_shared<FastaFetcher>(tiny("mini_hap.fa"));      // GT-AG at {20,40}
  auto noncanonical = std::make_shared<FastaFetcher>(tiny("mini_ref.fa"));   // GG-AG
  TempFasta other("other_contig.fa", "chrX", std::string(60, 'A'));
  auto missing = std::make_shared<FastaFetcher>(other.path);                 // no chrT at all
  const Junction intron{20, 40};
  // reference canonical, one haplotype non-canonical, the other unreadable: undecided
  HaplotypeProvider partial(canonical, {noncanonical, missing});
  CHECK(partial.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kNotEvaluable);
  // both haplotypes non-canonical and readable: disrupted
  HaplotypeProvider both(canonical, {noncanonical, noncanonical});
  CHECK(both.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kDisrupted);
  // a canonical haplotype is positive evidence even when another is unreadable
  HaplotypeProvider created(noncanonical, {missing, canonical});
  CHECK(created.motifs("chrT", intron, Strand::kPlus).verdict == VariantVerdict::kCreated);
}
