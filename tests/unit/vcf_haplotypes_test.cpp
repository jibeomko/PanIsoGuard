#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <memory>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/variant_motif.hpp"
#include "test_util.hpp"

using namespace panisoguard;

// two_jx_ref.fa: chrT2 with two non-canonical introns {20,40} and {60,80} (GG donors).
// two_jx_split.vcf.gz: S1 carries the G>T at 1-based 22 on haplotype 1 and the one at 62 on
// haplotype 2 (= two_jx_hapA.fa and two_jx_hapB.fa); two_jx_joint.vcf.gz puts both on
// haplotype 1 (= two_jx_hapC.fa + reference); two_jx_unphased.vcf.gz has both as 0/1 plus an
// insertion inside the first donor window; two_jx_twosample.vcf.gz has samples S1 and S2.

TEST_CASE("VcfHaplotypes applies the sample's alleles by haplotype, like bcftools consensus -H 1/2", "[vcf]") {
  VcfHaplotypes vcf(tiny("two_jx_split.vcf.gz"));
  CHECK(vcf.n_haplotypes() == 2);
  CHECK(vcf.has_chrom("chrT2"));
  CHECK_FALSE(vcf.has_chrom("chrZ"));
  bool unphased = false;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", unphased) == std::vector<std::string>{"GT", "GG"});
  CHECK(vcf.overlay("chrT2", 60, 62, "GG", unphased) == std::vector<std::string>{"GG", "GT"});
  CHECK(vcf.overlay("chrT2", 38, 40, "AG", unphased) == std::vector<std::string>{"AG", "AG"});   // no variant: reference
  CHECK_FALSE(unphased);
  CHECK(vcf.overlay("chrZ", 20, 22, "GG", unphased).empty());                                       // contig not in the header
  CHECK(vcf.n_unphased_het_applied() == 0);
  CHECK(vcf.n_non_snv_skipped() == 0);
}

TEST_CASE("a VCF-backed HaplotypeProvider gives the FASTA provider's motifs and verdicts", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  HaplotypeProvider from_fasta(ref, {std::make_shared<FastaFetcher>(tiny("two_jx_hapA.fa")),
                                     std::make_shared<FastaFetcher>(tiny("two_jx_hapB.fa"))});
  HaplotypeProvider from_vcf(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_split.vcf.gz"))});
  CHECK(from_vcf.n_haplotypes() == 2);
  for (const Junction intron : {Junction{20, 40}, Junction{60, 80}}) {
    const JunctionMotifs a = from_fasta.motifs("chrT2", intron, Strand::kPlus);
    const JunctionMotifs b = from_vcf.motifs("chrT2", intron, Strand::kPlus);
    REQUIRE(b.haps.size() == a.haps.size());
    CHECK(b.ref.motif == a.ref.motif);
    for (std::size_t h = 0; h < a.haps.size(); ++h) {
      CHECK(b.haps[h].motif == a.haps[h].motif);
      CHECK(b.haps[h].canonical == a.haps[h].canonical);
    }
    CHECK(b.verdict == a.verdict);
    CHECK(b.verdict == VariantVerdict::kCreated);
    CHECK(b.phase == MotifPhase::kPhased);
    CHECK(a.phase == MotifPhase::kUnknown);
  }
  // a contig the VCF does not know: every haplotype unreadable, so nothing is decided
  const JunctionMotifs z = from_vcf.motifs("chrT2", {20, 40}, Strand::kUnknown);
  CHECK(z.verdict == VariantVerdict::kNotEvaluable);
}

TEST_CASE("adjudicator with --haplotype-vcf: one haplotype must carry every novel junction", "[vcf][adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  SqantiTable sqanti;
  SqantiRecord r;
  r.isoform = "isoBoth"; r.chrom = "chrT2"; r.strand = Strand::kPlus;
  r.structural_category = "novel_not_in_catalog"; r.all_canonical = "non_canonical";
  sqanti.records.push_back(r);
  std::map<std::string, IntronChain> chains;
  IntronChain c; c.chrom = "chrT2"; c.strand = Strand::kPlus; c.introns = {{20, 40}, {60, 80}};
  chains["isoBoth"] = c;
  AdjudicateInputs in;
  in.sqanti = &sqanti; in.chains = &chains; in.catalog = &cat; in.haplotype_circular = false;
  std::vector<JunctionRow> jx;

  SECTION("junctions on different haplotypes: held") {
    HaplotypeProvider split(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_split.vcf.gz"))});
    in.haplotype = &split;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    REQUIRE(out.size() == 1);
    CHECK(out[0].evidence.n_novel_jx_variant_created == 2);
    CHECK_FALSE(out[0].evidence.haplotype_consistent);
    CHECK(out[0].verdict.confidence == ConfidenceClass::kUnresolved);
    CHECK(out[0].verdict.verdict_reason == "haplotype_inconsistent");
    REQUIRE(jx.size() == 2);
    CHECK(jx[0].variant_phase == "phased");
    CHECK(jx[0].canonical_on == std::vector<int>{1});
    CHECK(jx[1].canonical_on == std::vector<int>{2});
  }
  SECTION("both on haplotype 1: rescued") {
    HaplotypeProvider joint(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_joint.vcf.gz"))});
    in.haplotype = &joint;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    CHECK(out[0].evidence.haplotype_consistent);
    CHECK(out[0].verdict.confidence == ConfidenceClass::kReferenceBias);
    CHECK(out[0].evidence.n_novel_jx_phase_uncertain == 0);
  }
  SECTION("unphased genotypes: applied by allele order, flagged, indel skipped") {
    auto vcf = std::make_shared<VcfHaplotypes>(tiny("two_jx_unphased.vcf.gz"));
    HaplotypeProvider unph(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});
    in.haplotype = &unph;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    // 0/1 at both sites puts both T alleles on haplotype 2, so the chain is "consistent" only by allele order
    CHECK(out[0].verdict.confidence == ConfidenceClass::kReferenceBias);
    CHECK(out[0].evidence.n_novel_jx_phase_uncertain == 2);
    REQUIRE(jx.size() == 2);
    CHECK(jx[0].variant_phase == "unphased");
    CHECK(jx[0].canonical_on == std::vector<int>{2});
    CHECK(vcf->n_unphased_het_applied() == 2);
    CHECK(vcf->n_non_snv_skipped() == 1);   // the G>GA insertion inside the first donor window
    bool found = false;
    for (const auto& t : out[0].verdict.rule_trace) if (t.find("unphased allele") != std::string::npos) found = true;
    CHECK(found);
  }
}

TEST_CASE("VcfHaplotypes refuses what it cannot read safely", "[vcf]") {
  REQUIRE_THROWS_WITH(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz")), Catch::Contains("--vcf-sample"));
  CHECK_NOTHROW(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz"), "S2"));
  REQUIRE_THROWS_WITH(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz"), "S9"), Catch::Contains("S9"));
  // S2 of the two-sample file carries both variants on haplotype 1
  VcfHaplotypes s2(tiny("two_jx_twosample.vcf.gz"), "S2");
  bool unphased = false;
  CHECK(s2.overlay("chrT2", 20, 22, "GG", unphased) == std::vector<std::string>{"GT", "GG"});
  // an uncompressed, unindexed VCF is refused
  const std::string plain = scratch("plain.vcf");
  { std::ofstream(plain) << "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"; }
  REQUIRE_THROWS_WITH(VcfHaplotypes(plain), Catch::Contains("index"));
  std::remove(plain.c_str());
  REQUIRE_THROWS(VcfHaplotypes(tiny("does_not_exist.vcf.gz")));
}
