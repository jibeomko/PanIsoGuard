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
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"GT", "GG"});
  CHECK(vcf.overlay("chrT2", 60, 62, "GG", wp) == std::vector<std::string>{"GG", "GT"});
  CHECK(vcf.overlay("chrT2", 38, 40, "AG", wp) == std::vector<std::string>{"AG", "AG"});   // no variant: reference
  CHECK_FALSE(wp.unphased);
  CHECK(vcf.overlay("chrZ", 20, 22, "GG", wp).empty());                                       // contig not in the header
  CHECK(vcf.n_unphased_het_applied() == 0);
  CHECK(vcf.n_unapplied_alleles() == 0);
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
  SECTION("unphased genotypes: applied by allele order, flagged") {
    auto vcf = std::make_shared<VcfHaplotypes>(tiny("two_jx_unphased_snv.vcf.gz"));
    HaplotypeProvider unph(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});
    in.haplotype = &unph;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    // 0/1 at both sites puts both T alleles on haplotype 2 by allele order only: each junction is created
    // on a haplotype, but the file does not say it is the same one, so the chain is held
    CHECK(out[0].verdict.confidence == ConfidenceClass::kUnresolved);
    CHECK(out[0].verdict.verdict_reason == "haplotype_phase_unknown");
    CHECK(out[0].evidence.n_novel_jx_variant_created == 2);
    CHECK(out[0].evidence.haplotype_phase_unknown);
    CHECK(out[0].evidence.n_novel_jx_phase_uncertain == 2);
    REQUIRE(jx.size() == 2);
    CHECK(jx[0].variant_phase == "unphased");
    CHECK(jx[0].canonical_on == std::vector<int>{2});
    CHECK(vcf->n_unphased_het_applied() == 2);
    bool found = false;
    for (const auto& t : out[0].verdict.rule_trace) if (t.find("one common haplotype is not established") != std::string::npos) found = true;
    CHECK(found);
  }
  SECTION("an unphased insertion inside the donor next to an unphased substitution: the donor is not known") {
    // two_jx_unphased.vcf.gz adds G>GA (0/1) between the two donor bases of the first intron. By
    // allele order both alleles sit on haplotype 2, whose donor then reads GA; with the insertion
    // on the other haplotype the substitution alone would make GT. The file does not say which, so
    // no haplotype's donor is known there; the second junction is still created.
    auto vcf = std::make_shared<VcfHaplotypes>(tiny("two_jx_unphased.vcf.gz"));
    HaplotypeProvider unph(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});
    in.haplotype = &unph;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    CHECK(out[0].evidence.n_novel_jx_variant_created == 1);
    CHECK(out[0].verdict.confidence == ConfidenceClass::kUnresolved);
    REQUIRE(jx.size() == 2);
    CHECK(jx[0].variant_phase == "unlinked");
    CHECK(jx[0].canonical_on.empty());
    CHECK(jx[1].canonical_on == std::vector<int>{2});
    CHECK(vcf->n_indels_applied() == 1);
    CHECK(vcf->n_unapplied_alleles() == 0);
  }
}

TEST_CASE("VcfHaplotypes refuses what it cannot read safely", "[vcf]") {
  REQUIRE_THROWS_WITH(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz")), Catch::Contains("--vcf-sample"));
  CHECK_NOTHROW(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz"), "S2"));
  REQUIRE_THROWS_WITH(VcfHaplotypes(tiny("two_jx_twosample.vcf.gz"), "S9"), Catch::Contains("S9"));
  // S2 of the two-sample file carries both variants on haplotype 1
  VcfHaplotypes s2(tiny("two_jx_twosample.vcf.gz"), "S2");
  WindowPhase wp;
  CHECK(s2.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"GT", "GG"});
  // an uncompressed, unindexed VCF is refused
  const std::string plain = scratch("plain.vcf");
  { std::ofstream(plain) << "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"; }
  REQUIRE_THROWS_WITH(VcfHaplotypes(plain), Catch::Contains("index"));
  std::remove(plain.c_str());
  REQUIRE_THROWS(VcfHaplotypes(tiny("does_not_exist.vcf.gz")));
}

// two_jx_edges.vcf.gz: a haploid first record (pos 5), an MNP over the first donor (21-22, 1|0), a `*`
// allele (40) and a <DEL> (42) at motif positions, a multiallelic site whose second ALT is carried
// (62, 2|0), a missing genotype (80, ./.) and two records at one position (100: C on haplotype 1, G on 2).
TEST_CASE("VcfHaplotypes edge cases: ploidy from several records, MNPs, alleles without a sequence, missing and stacked genotypes", "[vcf]") {
  VcfHaplotypes vcf(tiny("two_jx_edges.vcf.gz"));
  CHECK(vcf.n_haplotypes() == 2);   // not the haploid first record's 1
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"GT", "GG"});   // MNP applied base by base
  CHECK(vcf.overlay("chrT2", 60, 62, "GG", wp) == std::vector<std::string>{"GC", "GG"});   // allele 2 of a multiallelic site
  // `*`: the base is removed by a deletion the file does not carry, so it is unknown, not reference
  CHECK(vcf.overlay("chrT2", 38, 40, "AG", wp) == std::vector<std::string>{"AN", "AG"});
  CHECK(vcf.n_unapplied_alleles() == 1);
  // <DEL> without END or SVLEN: how far it reaches is not known, so the bases from it on are not
  CHECK(vcf.overlay("chrT2", 40, 42, "AA", wp) == std::vector<std::string>{"AN", "AA"});
  CHECK(vcf.n_unapplied_alleles() == 2);
  CHECK(vcf.overlay("chrT2", 78, 80, "AG", wp) == std::vector<std::string>{"AN", "AN"});   // ./. is unknown, not reference
  CHECK(wp.missing);
  CHECK(vcf.overlay("chrT2", 99, 101, "AA", wp) == std::vector<std::string>{"CA", "GA"});  // two records, one position
  CHECK_FALSE(wp.unphased);
  CHECK(vcf.n_unphased_het_applied() == 0);
}

TEST_CASE("a VCF without the junction's contig leaves the variant axis not evaluable", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("mini_ref.fa"));   // chr1; the VCF header names chrT2 only
  HaplotypeProvider from_vcf(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_split.vcf.gz"))});
  const JunctionMotifs m = from_vcf.motifs("chr1", {200, 300}, Strand::kPlus);
  REQUIRE(m.haps.size() == 2);
  CHECK_FALSE(m.haps[0].evaluable);
  CHECK_FALSE(m.haps[1].evaluable);
  CHECK(m.verdict == VariantVerdict::kNotEvaluable);
  CHECK(m.phase == MotifPhase::kUnknown);
}

// two_jx_missing.vcf.gz: at the first donor the genotype is 0|. (allele missing on haplotype 2) and
// 0|1 one base on; haplotype 2 could carry GT-AG (the T) but the missing base is unknown, so
// nothing may be called canonical there and the junction is not rescued.
TEST_CASE("a missing allele at a motif base leaves that haplotype's motif unknown, not reference", "[vcf][variant]") {
  VcfHaplotypes vcf(tiny("two_jx_missing.vcf.gz"));
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"GG", "NT"});
  CHECK(wp.missing);
  CHECK(vcf.n_missing_alleles() == 1);
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  HaplotypeProvider prov(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_missing.vcf.gz"))});
  const JunctionMotifs m = prov.motifs("chrT2", {20, 40}, Strand::kPlus);
  CHECK(m.haps[0].evaluable); CHECK_FALSE(m.haps[0].canonical);
  CHECK_FALSE(m.haps[1].evaluable);
  CHECK(m.verdict == VariantVerdict::kNotEvaluable);   // no haplotype is known to be canonical
}

// two_jx_haploid.vcf.gz: a diploid record elsewhere (ploidy 2) and a haploid genotype `1` that breaks the
// canonical acceptor of the first intron; the absent second allele is unknown, not a surviving reference
// copy, so the junction is not "none" but not evaluable.
TEST_CASE("a haploid genotype in a diploid file does not invent a reference haplotype", "[vcf][variant]") {
  VcfHaplotypes vcf(tiny("two_jx_haploid.vcf.gz"));
  CHECK(vcf.n_haplotypes() == 2);
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 38, 40, "AG", wp) == std::vector<std::string>{"CG", "NG"});
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  HaplotypeProvider prov(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_haploid.vcf.gz"))});
  // the reference intron {20,40} is GG-AG; with the acceptor broken on the only present allele
  const JunctionMotifs m = prov.motifs("chrT2", {20, 40}, Strand::kPlus);
  CHECK(m.haps[0].evaluable); CHECK_FALSE(m.haps[0].canonical);
  CHECK_FALSE(m.haps[1].evaluable);
}

// two_jx_phasesets*.vcf.gz: both novel junctions are created by heterozygous alleles written first
// (1|0). With PS 1 and PS 2 the file does not say the two alleles share a haplotype, so the rescue is
// held; flipping the second block's order changes nothing; one phase set (PS 7 twice) rescues.
TEST_CASE("adjudicator with --haplotype-vcf: a common haplotype needs one phase set", "[vcf][adjudicator]") {
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
  auto run = [&](const char* vcf) {
    HaplotypeProvider hp(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny(vcf))});
    in.haplotype = &hp;
    return adjudicate(in, RuleEngine(), nullptr).at(0);
  };
  const AdjudicationResult two_sets = run("two_jx_phasesets.vcf.gz");
  CHECK(two_sets.evidence.n_novel_jx_variant_created == 2);
  CHECK(two_sets.evidence.haplotype_phase_unknown);
  CHECK(two_sets.verdict.confidence == ConfidenceClass::kUnresolved);
  CHECK(two_sets.verdict.verdict_reason == "haplotype_phase_unknown");
  const AdjudicationResult flipped = run("two_jx_phasesets_flip.vcf.gz");
  CHECK(flipped.verdict.confidence == ConfidenceClass::kUnresolved);          // the same information, the same call
  CHECK(flipped.verdict.verdict_reason == "haplotype_phase_unknown");
  const AdjudicationResult one_set = run("two_jx_sameset.vcf.gz");
  CHECK_FALSE(one_set.evidence.haplotype_phase_unknown);
  CHECK(one_set.verdict.confidence == ConfidenceClass::kReferenceBias);
  const AdjudicationResult no_ps = run("two_jx_joint.vcf.gz");                 // phased, no PS: one default set
  CHECK(no_ps.verdict.confidence == ConfidenceClass::kReferenceBias);
}

// two_jx_twohet_*.vcf.gz: the second intron {60,80} (reference GG-AG) gets a donor G>T (62) and an acceptor
// A>C (79), both heterozygous. In phase sets 1 and 2 the file does not say which haplotype carries both,
// so no haplotype's motif is known, whichever way the alleles are written; in one set it does.
TEST_CASE("two heterozygous records inside one motif need one phase set", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  auto motifs = [&](const char* vcf) {
    HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny(vcf))});
    return p.motifs("chrT2", {60, 80}, Strand::kPlus);
  };
  for (const char* f : {"two_jx_twohet_sets.vcf.gz", "two_jx_twohet_flip.vcf.gz", "two_jx_twohet_unphased.vcf.gz"}) {
    const JunctionMotifs m = motifs(f);
    CHECK(m.phase_unlinked);
    CHECK_FALSE(m.haps[0].evaluable); CHECK_FALSE(m.haps[1].evaluable);
    CHECK(m.verdict == VariantVerdict::kNotEvaluable);
  }
  const JunctionMotifs one_set = motifs("two_jx_twohet_sameps.vcf.gz");   // 1|0 and 0|1 in set 1: haplotype 1 is GT-AG
  CHECK_FALSE(one_set.phase_unlinked);
  CHECK(one_set.haps[0].motif == "GT-AG"); CHECK(one_set.haps[1].motif == "GG-CG");
  CHECK(one_set.verdict == VariantVerdict::kCreated);
  // a single heterozygous record is certain on its own, whichever haplotype carries it
  const JunctionMotifs single = motifs("two_jx_unphased.vcf.gz");
  CHECK_FALSE(single.phase_unlinked);
  CHECK(single.verdict == VariantVerdict::kCreated);
}

// two_jx_nogt.vcf.gz: the first donor base has a record without GT (FORMAT DP only), the second a 0|1.
// A record that names a variant but no genotype leaves the base unknown on every haplotype.
TEST_CASE("a record without GT at a motif base is unknown, not reference", "[vcf][variant]") {
  VcfHaplotypes vcf(tiny("two_jx_nogt.vcf.gz"));
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"NG", "NT"});
  CHECK(wp.missing);
  CHECK(vcf.n_missing_alleles() == 2);
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny("two_jx_nogt.vcf.gz"))});
  CHECK(p.motifs("chrT2", {20, 40}, Strand::kPlus).verdict == VariantVerdict::kNotEvaluable);
}

// two_jx_partial_*.vcf.gz: partly missing genotypes. 1|. places its known ALT on haplotype 1 only by the
// writing, so two such records (or one with a fully heterozygous one) need one phase set like any two
// heterozygous records; alone, a 1|. is certain (the ALT exists on a haplotype); 1/. is unphased.
TEST_CASE("partly missing genotypes place alleles by index and need one phase set too", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  auto motifs = [&](const char* vcf, Junction j) {
    HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny(vcf))});
    return p.motifs("chrT2", j, Strand::kPlus);
  };
  // donor 1|.:1 (G>T) and acceptor 0|.:2 (A>C, haplotype 1 keeps A): GT-AG on haplotype 1 by writing only
  for (const char* f : {"two_jx_partial_sets.vcf.gz", "two_jx_partial_unphased.vcf.gz"}) {
    const JunctionMotifs m = motifs(f, {60, 80});
    CHECK(m.phase_unlinked);
    CHECK(m.verdict == VariantVerdict::kNotEvaluable);
  }
  const JunctionMotifs same = motifs("two_jx_partial_sameps.vcf.gz", {60, 80});   // one phase set: the ALT and the REF are linked
  CHECK_FALSE(same.phase_unlinked);
  CHECK(same.haps[0].motif == "GT-AG");
  CHECK(same.verdict == VariantVerdict::kCreated);
  // a reference allele placed by a partly missing record (0|.:1) next to a heterozygous one (1|0:2)
  const JunctionMotifs ref_sets = motifs("two_jx_partial_ref.vcf.gz", {20, 40});
  CHECK(ref_sets.phase_unlinked);
  CHECK(ref_sets.verdict == VariantVerdict::kNotEvaluable);
  const JunctionMotifs ref_same = motifs("two_jx_partial_ref_sameps.vcf.gz", {20, 40});
  CHECK_FALSE(ref_same.phase_unlinked);
  CHECK(ref_same.haps[0].motif == "GT-AG");
  CHECK(ref_same.verdict == VariantVerdict::kCreated);
  // a single 1|. is certain on its own
  const JunctionMotifs single = motifs("two_jx_partial_single.vcf.gz", {60, 80});
  CHECK_FALSE(single.phase_unlinked);
  CHECK(single.verdict == VariantVerdict::kCreated);
}

// two_jx_mnp_span*.vcf.gz: one 20-base MNP over the whole second intron {60,80}, changing the donor GG to GT
// and nothing else; it is read once per window but is one record, so an unphased 0/1 is a single variant
// (certain on its own), not two unlinked ones.
TEST_CASE("an MNP spanning both motif windows is one record, phased or not", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  for (const char* f : {"two_jx_mnp_span.vcf.gz", "two_jx_mnp_span_phased.vcf.gz"}) {
    HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny(f))});
    const JunctionMotifs m = p.motifs("chrT2", {60, 80}, Strand::kPlus);
    CHECK_FALSE(m.phase_unlinked);
    CHECK(m.haps[1].motif == "GT-AG");
    CHECK(m.verdict == VariantVerdict::kCreated);
  }
}

// two_jx_samepos*.vcf.gz: two records at one position of the second donor (G>T 1|0 and G>C 0|1); in different
// phase sets they are two unlinked records whichever order the rows are in, in one set they are linked.
TEST_CASE("two records at one position are two records, in either row order", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  auto motifs = [&](const char* vcf) {
    HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{std::make_shared<VcfHaplotypes>(tiny(vcf))});
    return p.motifs("chrT2", {60, 80}, Strand::kPlus);
  };
  for (const char* f : {"two_jx_samepos_sets.vcf.gz", "two_jx_samepos_sets_swapped.vcf.gz"}) {
    const JunctionMotifs m = motifs(f);
    CHECK(m.phase_unlinked);
    CHECK(m.verdict == VariantVerdict::kNotEvaluable);
  }
  const JunctionMotifs one_set = motifs("two_jx_samepos_sameps.vcf.gz");
  CHECK_FALSE(one_set.phase_unlinked);
  CHECK(one_set.haps[0].motif == "GT-AG");
  CHECK(one_set.verdict == VariantVerdict::kCreated);
}

// two_jx_dup_*.vcf.gz: the second donor site written twice with the same alleles. With another genotype
// (1|0:1 then 0|1:2) the file contradicts itself and the bases are unknown in either row order; an identical
// repeat is one record.
TEST_CASE("a repeated site with another genotype is a contradiction, not two haplotypes", "[vcf][variant]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  // both row orders; the first genotype written again after the contradiction; another ALT over the same base
  for (const char* f : {"two_jx_dup_conflict.vcf.gz", "two_jx_dup_conflict_swapped.vcf.gz", "two_jx_dup_conflict_aba.vcf.gz",
                        "two_jx_dup_conflict_otheralt.vcf.gz"}) {
    auto vcf = std::make_shared<VcfHaplotypes>(tiny(f));
    HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});
    const JunctionMotifs m = p.motifs("chrT2", {60, 80}, Strand::kPlus);
    CHECK_FALSE(m.haps[0].evaluable); CHECK_FALSE(m.haps[1].evaluable);
    CHECK(m.verdict == VariantVerdict::kNotEvaluable);
    CHECK(vcf->n_conflicting_records() == 1);
  }
  auto same = std::make_shared<VcfHaplotypes>(tiny("two_jx_dup_identical.vcf.gz"));
  HaplotypeProvider p(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{same});
  const JunctionMotifs m = p.motifs("chrT2", {60, 80}, Strand::kPlus);
  CHECK(m.haps[0].motif == "GT-AG");
  CHECK(m.verdict == VariantVerdict::kCreated);
  CHECK(same->n_conflicting_records() == 0);
}

TEST_CASE("VcfHaplotypes leaves unknown what a record cannot say", "[vcf]") {
  // two_jx_defects.vcf.gz: an allele index the record does not define (22 G>T, 2|0); an MNP and
  // a substitution that give haplotype 1 different bases at one position (39 AG>CT, 40 G>A); a
  // REF that is not the reference (61 A>T where the reference has G); a `*` allele under the
  // deletion that removes its base (86 AAAA>A, 88 A>*).
  VcfHaplotypes vcf(tiny("two_jx_defects.vcf.gz"));
  WindowPhase wp;
  CHECK(vcf.overlay("chrT2", 20, 22, "GG", wp) == std::vector<std::string>{"GN", "GG"});
  CHECK(vcf.n_missing_alleles() == 1);
  CHECK(vcf.overlay("chrT2", 38, 40, "AG", wp) == std::vector<std::string>{"NN", "AG"});
  CHECK(vcf.n_conflicting_records() == 1);
  CHECK(vcf.overlay("chrT2", 60, 62, "GG", wp) == std::vector<std::string>{"NG", "NG"});
  CHECK(vcf.n_ref_mismatches() == 1);
  CHECK(wp.missing);
  CHECK(vcf.overlay("chrT2", 85, 90, "AAAAA", wp) == std::vector<std::string>{"ANNNA", "AAAAA"});
  CHECK(vcf.n_unapplied_alleles() == 0);   // the deletion speaks for the `*`
}

// indel_ref.fa, indel.vcf.gz and indel.bam: seven introns on chrI, one case each
// (tests/data/tiny/make_indel_fixtures.py lists them).
TEST_CASE("insertions and deletions are applied at splice motifs", "[vcf][indel]") {
  auto ref = std::make_shared<FastaFetcher>(tiny("indel_ref.fa"));
  auto vcf = std::make_shared<VcfHaplotypes>(tiny("indel.vcf.gz"));
  HaplotypeProvider hp(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});

  SECTION("a deletion or an insertion inside the motif") {
    const JunctionMotifs a = hp.motifs("chrI", {40, 70}, Strand::kPlus);
    CHECK(a.ref.motif == "GA-AG");
    CHECK(a.verdict == VariantVerdict::kCreated);
    CHECK(a.haps[0].motif == "GT-AG");
    CHECK(a.indel);
    const JunctionMotifs b = hp.motifs("chrI", {110, 140}, Strand::kPlus);
    CHECK(b.ref.motif == "GT-TG");
    CHECK(b.verdict == VariantVerdict::kCreated);
    CHECK(b.haps[1].motif == "GT-AG");
    CHECK(vcf->n_indels_applied() == 2);
  }
  SECTION("an insertion exactly at the exon boundary is placed by the reads") {
    JunctionMotifs c = hp.motifs("chrI", {180, 210}, Strand::kPlus);
    CHECK(c.verdict == VariantVerdict::kNotEvaluable);   // intronic or exonic: the file cannot tell
    CHECK(c.boundary_insertion_right == 4);
    CHECK(c.boundary_insertion_left == 0);
    c = hp.motifs("chrI", {180, 210}, Strand::kPlus, [](bool, int) { return BoundaryReads::kWithout; });
    CHECK(c.verdict == VariantVerdict::kCreated);        // the reads lack it: it is intronic and ends in AG
    CHECK(c.haps[0].motif == "GT-AG");
    CHECK(c.boundary_insertion_right == 0);
    c = hp.motifs("chrI", {180, 210}, Strand::kPlus, [](bool left, int len) {
      CHECK_FALSE(left);
      CHECK(len == 4);
      return BoundaryReads::kWith;
    });
    CHECK(c.verdict == VariantVerdict::kNone);           // the reads carry it: the intron ends as on the reference
    CHECK(c.haps[0].motif == "GT-TC");
  }
  SECTION("a deletion equal to the intron, wherever the reference lets it be written") {
    const JunctionMotifs d = hp.motifs("chrI", {320, 350}, Strand::kPlus);   // the file writes [321, 351)
    CHECK(d.verdict == VariantVerdict::kDeleted);
    CHECK(d.deleted == std::vector<bool>{false, true});
    CHECK(d.phase == MotifPhase::kPhased);
    CHECK(vcf->n_deletions_equal_to_junctions() == 1);
    CHECK(hp.motifs("chrI", {321, 351}, Strand::kPlus).verdict == VariantVerdict::kDeleted);
    CHECK(hp.motifs("chrI", {322, 352}, Strand::kPlus).verdict != VariantVerdict::kDeleted);
    CHECK(hp.motifs("chrI", {320, 349}, Strand::kPlus).verdict != VariantVerdict::kDeleted);
  }
  SECTION("a haplotype FASTA of another contig length is not read") {
    HaplotypeProvider fa(ref, std::vector<std::shared_ptr<FastaFetcher>>{std::make_shared<FastaFetcher>(tiny("indel_hap_short.fa"))});
    const JunctionMotifs e = fa.motifs("chrI", {390, 420}, Strand::kPlus);
    CHECK(e.verdict == VariantVerdict::kNotEvaluable);
    CHECK_FALSE(e.haps[0].evaluable);
    fa.motifs("chrI", {460, 490}, Strand::kPlus);
    CHECK(fa.n_length_mismatches() == 1);   // once per file and contig
  }
}

TEST_CASE("adjudicator: indels, a deletion equal to a junction, and the haplotype of the junction reads", "[vcf][indel][adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  auto ref = std::make_shared<FastaFetcher>(tiny("indel_ref.fa"));
  auto vcf = std::make_shared<VcfHaplotypes>(tiny("indel.vcf.gz"));
  HaplotypeProvider hp(ref, std::vector<std::shared_ptr<VcfHaplotypes>>{vcf});
  BamReader bam(tiny("indel.bam"));
  SqantiTable sqanti;
  std::map<std::string, IntronChain> chains;
  const std::vector<std::pair<std::string, int64_t>> introns = {{"A", 40}, {"B", 110}, {"C", 180}, {"G", 250}, {"D", 320}, {"E", 390}, {"F", 460}};
  for (const auto& kv : introns) {
    SqantiRecord r; r.isoform = kv.first; r.chrom = "chrI"; r.strand = Strand::kPlus;
    r.structural_category = "novel_not_in_catalog"; r.all_canonical = "non_canonical";
    sqanti.records.push_back(r);
    IntronChain c; c.chrom = "chrI"; c.strand = Strand::kPlus; c.introns = {{kv.second, kv.second + 30}};
    chains[kv.first] = c;
  }
  AdjudicateInputs in;
  in.sqanti = &sqanti; in.chains = &chains; in.catalog = &cat; in.haplotype = &hp;
  std::vector<JunctionRow> jx;   // in genomic order: A B C G D E F
  auto traced = [](const AdjudicationResult& r, const char* what) {
    for (const auto& t : r.verdict.rule_trace) if (t.find(what) != std::string::npos) return true;
    return false;
  };

  SECTION("with the reads") {
    in.bam = &bam;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    REQUIRE(out.size() == 7);
    REQUIRE(jx.size() == 7);
    for (const std::size_t i : {0u, 1u, 2u}) {   // A, B: the indel makes the motif; C: the reads lack the insertion
      CHECK(out[i].verdict.confidence == ConfidenceClass::kReferenceBias);
      CHECK(out[i].verdict.verdict_reason.rfind("haplotype_variant", 0) == 0);
      CHECK(jx[i].variant == "created");
      CHECK(jx[i].read_phase == "NA");             // on both haplotypes: nothing to phase
    }
    CHECK(out[3].verdict.confidence != ConfidenceClass::kReferenceBias);   // G: the reads carry the insertion
    CHECK(jx[3].variant == "none");
    CHECK(out[4].verdict.confidence == ConfidenceClass::kReferenceBias);   // D
    CHECK(out[4].verdict.verdict_reason.rfind("haplotype_deletion", 0) == 0);
    CHECK(out[4].evidence.n_novel_jx_variant_deleted == 1);
    CHECK(traced(out[4], "deletion equal to every novel junction"));
    CHECK(jx[4].variant == "deleted");
    CHECK(jx[4].hap_motifs == std::vector<std::string>{"CC-TT", "DEL"});
    CHECK(jx[4].read_phase == "consistent");
    CHECK(jx[4].read_phase_n_motif == 4);
    CHECK(out[5].verdict.confidence == ConfidenceClass::kReferenceBias);   // E: reads of the haplotype with the motif
    CHECK(jx[5].read_phase == "consistent");
    CHECK(jx[5].read_phase_n_other == 0);
    CHECK(out[6].verdict.confidence != ConfidenceClass::kReferenceBias);   // F: reads of the other haplotype
    CHECK(out[6].evidence.n_novel_jx_read_phase_contradicted == 1);
    CHECK(out[6].evidence.n_novel_jx_variant_created == 0);
    CHECK(traced(out[6], "not counted as explained"));
    CHECK(jx[6].variant == "created");
    CHECK(jx[6].read_phase == "contradicted");
    CHECK(jx[6].read_phase_n_motif == 0);
    CHECK(jx[6].read_phase_n_other == 4);
  }
  SECTION("without reads the boundary insertion stays open and nothing contradicts") {
    const auto out = adjudicate(in, RuleEngine(), &jx);
    CHECK(out[2].verdict.confidence != ConfidenceClass::kReferenceBias);   // C
    CHECK(jx[2].variant == "NA");
    CHECK(out[6].verdict.confidence == ConfidenceClass::kReferenceBias);   // F
    CHECK(jx[6].read_phase == "NA");
  }
}

TEST_CASE("read_phase assigns a read only by a margin of sites", "[adjudicator]") {
  JunctionRead r;
  r.start = 0; r.end = 60;
  r.blocks = {{0, 0, 10}, {50, 10, 10}};
  r.seq = "ACGTACGTACACGTACGTAC";
  CHECK(r.base_at(3) == 'T');
  CHECK(r.base_at(10) == 0);    // inside the intron
  CHECK(r.base_at(52) == 'G');
  const std::vector<bool> on = {false, true};
  // the read carries haplotype 2 at both sites: assigned
  std::vector<PhaseSite> sites = {{1, "AC"}, {52, "TG"}};
  ReadPhase p = read_phase({r}, sites, on, 2);
  CHECK(p.n_motif == 1);
  CHECK(p.n_other == 0);
  // one site each way, or one site only: not assigned
  sites = {{1, "CA"}, {52, "TG"}};
  p = read_phase({r}, sites, on, 2);
  CHECK(p.n_motif + p.n_other == 0);
  sites = {{52, "TG"}, {30, "AC"}};
  CHECK(read_phase({r}, sites, on, 2).n_motif == 0);
}
