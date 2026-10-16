#include "catch2/catch.hpp"

#include <cstdio>
#include <map>
#include <memory>
#include <string>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/bam_features.hpp"
#include "panisoguard/variant_motif.hpp"
#include "panisoguard/consensus.hpp"
#include "panisoguard/gtf.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
SqantiRecord rec(const std::string& id, const std::string& cat, const std::string& all_canon) {
  SqantiRecord r;
  r.isoform = id;
  r.chrom = "chr1";
  r.strand = Strand::kPlus;
  r.structural_category = cat;
  r.all_canonical = all_canon;
  return r;
}
IntronChain chain(std::vector<Junction> introns) {
  IntronChain c;
  c.chrom = "chr1";
  c.strand = Strand::kPlus;
  c.introns = std::move(introns);
  return c;
}
SjRecord sj(Junction j, int n_uniq, int motif = 1) {
  return SjRecord{"chr1", j, Strand::kPlus, motif, true, n_uniq, 0, 30};
}
const AdjudicationResult* find(const std::vector<AdjudicationResult>& v, const std::string& id) {
  for (const auto& r : v)
    if (r.isoform_id == id) return &r;
  return nullptr;
}
}  // namespace

TEST_CASE("adjudicator end-to-end: SR-supported novel vs unsupported artifact", "[adjudicator]") {
  // Catalog (mini.gtf) knows introns {200,300},{400,500} (T1) and {200,500} (T2).
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoNNC", "novel_not_in_catalog", "non_canonical"));
  sqanti.records.push_back(rec("isoNNC2", "novel_not_in_catalog", "canonical"));
  sqanti.records.push_back(rec("isoFSM", "full-splice_match", "canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoNNC"] = chain({{200, 350}, {400, 500}});   // novel jx {200,350}, no SR support
  chains["isoNNC2"] = chain({{200, 360}, {400, 500}});  // novel jx {200,360}, SR-supported below

  SjTable sjt;
  sjt.add(sj({400, 500}, 8));    // known junction (ignored by novelty check)
  sjt.add(sj({200, 360}, 10));   // supports isoNNC2's novel junction
  // note: {200,350} is intentionally absent -> isoNNC's novel junction unsupported

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.sj = &sjt;
  const auto out = adjudicate(in, RuleEngine());
  REQUIRE(out.size() == 3);

  const AdjudicationResult* nnc = find(out, "isoNNC");
  REQUIRE(nnc != nullptr);
  CHECK(nnc->evidence.n_novel_junctions == 1);
  CHECK(nnc->evidence.n_novel_jx_sr_supported == 0);
  CHECK(nnc->verdict.novelty_support == NoveltySupport::kUnsupported);
  CHECK(nnc->verdict.primary_mechanism == Mechanism::kNoncanonical);
  CHECK(nnc->verdict.confidence == ConfidenceClass::kArtifact);

  const AdjudicationResult* nnc2 = find(out, "isoNNC2");
  REQUIRE(nnc2 != nullptr);
  CHECK(nnc2->evidence.n_novel_junctions == 1);
  CHECK(nnc2->evidence.n_novel_jx_sr_supported == 1);
  CHECK(nnc2->verdict.novelty_support == NoveltySupport::kSupported);
  CHECK(nnc2->verdict.confidence == ConfidenceClass::kNovelConfirmed);

  const AdjudicationResult* fsm = find(out, "isoFSM");
  REQUIRE(fsm != nullptr);
  CHECK(fsm->verdict.confidence == ConfidenceClass::kKnown);
}

TEST_CASE("adjudicator: pangenome graph-supported novel junction -> REFERENCE_BIAS", "[adjudicator]") {
  // Catalog knows {200,300},{400,500} (T1) and {200,500} (T2); {200,350} is novel.
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPan", "novel_not_in_catalog", "non_canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoPan"] = chain({{200, 350}, {400, 500}});  // novel jx {200,350}

  // The pangenome carries that novel junction on a graph haplotype path.
  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  const auto out = adjudicate(in, RuleEngine());
  REQUIRE(out.size() == 1);

  const AdjudicationResult* p = find(out, "isoPan");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.n_novel_junctions == 1);
  CHECK(p->evidence.pangenome_evaluable);
  CHECK(p->evidence.pangenome_rescue);
  CHECK(p->evidence.n_novel_jx_pangenome == 1);
  CHECK(p->verdict.confidence == ConfidenceClass::kReferenceBias);
  CHECK(p->verdict.primary_mechanism == Mechanism::kPopulationKnown);
  CHECK_FALSE(p->verdict.circularity_flag);  // population reference data is non-circular

  // Disabling the pangenome axis must drop the rescue (no SR/variant evidence here).
  const auto out_off = adjudicate(in, RuleEngine().with_axis_disabled("pangenome"));
  CHECK(find(out_off, "isoPan")->verdict.confidence != ConfidenceClass::kReferenceBias);
}

TEST_CASE("adjudicator: pangenome partial support does NOT rescue (all-or-nothing)", "[adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPartial", "novel_not_in_catalog", "canonical"));

  std::map<std::string, IntronChain> chains;
  // BOTH introns are novel vs the catalog ({200,350} and {600,900}).
  chains["isoPartial"] = chain({{200, 350}, {600, 900}});

  // The pangenome carries only ONE of the two novel junctions.
  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  const auto out = adjudicate(in, RuleEngine());

  const AdjudicationResult* p = find(out, "isoPartial");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.n_novel_junctions == 2);
  CHECK(p->evidence.n_novel_jx_pangenome == 1);
  CHECK_FALSE(p->evidence.pangenome_rescue);  // 1 of 2 -> the isoform is not fully reference bias
  CHECK(p->verdict.confidence != ConfidenceClass::kReferenceBias);
}

TEST_CASE("adjudicator: sample-derived pangenome provenance is held (circularity firewall)", "[adjudicator]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoPan", "novel_not_in_catalog", "canonical"));

  std::map<std::string, IntronChain> chains;
  chains["isoPan"] = chain({{200, 350}, {400, 500}});  // single novel jx {200,350}

  PangenomeJunctions pan;
  pan.add(PangenomeJunctionRecord{"chr1", Junction{200, 350}, Strand::kPlus, 2});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.pangenome = &pan;
  in.pangenome_circular = true;  // junction set may derive from the sample's own reads
  const auto out = adjudicate(in, RuleEngine());

  const AdjudicationResult* p = find(out, "isoPan");
  REQUIRE(p != nullptr);
  CHECK(p->evidence.pangenome_rescue);  // all novel junctions are graph-supported...
  CHECK(p->evidence.pangenome_circular);
  CHECK(p->verdict.confidence == ConfidenceClass::kUnresolved);  // ...but held, not promoted
  CHECK(p->verdict.circularity_flag);
}

namespace {
IntronChain chain_on(const char* chrom, std::vector<Junction> introns) {
  IntronChain c;
  c.chrom = chrom;
  c.strand = Strand::kPlus;
  c.introns = std::move(introns);
  return c;
}
}  // namespace

TEST_CASE("adjudicator: variant rescue requires one haplotype carrying every novel junction",
          "[adjudicator][variant]") {
  // chrT2: two non-canonical junctions on the reference. hapA fixes only the first, hapB
  // only the second, hapC both. Neither junction is in the (chr1) catalog -> both novel.
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  auto hapA = std::make_shared<FastaFetcher>(tiny("two_jx_hapA.fa"));
  auto hapB = std::make_shared<FastaFetcher>(tiny("two_jx_hapB.fa"));
  auto hapC = std::make_shared<FastaFetcher>(tiny("two_jx_hapC.fa"));

  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoBoth", "novel_not_in_catalog", "non_canonical"));
  sqanti.records.push_back(rec("isoFirst", "novel_not_in_catalog", "non_canonical"));
  std::map<std::string, IntronChain> chains;
  chains["isoBoth"] = chain_on("chrT2", {{20, 40}, {60, 80}});
  chains["isoFirst"] = chain_on("chrT2", {{20, 40}});

  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.haplotype_circular = false;  // wgs provenance

  SECTION("junction 1 on hapA, junction 2 on hapB: every junction created, no rescue") {
    HaplotypeProvider split(ref, {hapA, hapB});
    in.haplotype = &split;
    const auto out = adjudicate(in, RuleEngine());
    const AdjudicationResult* both = find(out, "isoBoth");
    REQUIRE(both != nullptr);
    CHECK(both->evidence.n_novel_jx_variant_created == 2);  // the per-junction view still says 2/2
    CHECK_FALSE(both->evidence.haplotype_consistent);
    CHECK_FALSE(both->evidence.variant_rescue);
    CHECK(both->verdict.confidence == ConfidenceClass::kUnresolved);
    CHECK(both->verdict.verdict_reason == "haplotype_inconsistent");
    // A single created junction is trivially on one haplotype -> rescued as before.
    const AdjudicationResult* first = find(out, "isoFirst");
    REQUIRE(first != nullptr);
    CHECK(first->evidence.haplotype_consistent);
    CHECK(first->verdict.confidence == ConfidenceClass::kReferenceBias);
  }
  SECTION("one haplotype (hapC) carries both: rescued") {
    HaplotypeProvider joint(ref, {hapA, hapC});
    in.haplotype = &joint;
    const auto out = adjudicate(in, RuleEngine());
    const AdjudicationResult* both = find(out, "isoBoth");
    REQUIRE(both != nullptr);
    CHECK(both->evidence.haplotype_consistent);
    CHECK(both->verdict.confidence == ConfidenceClass::kReferenceBias);
    CHECK(both->verdict.verdict_reason == "haplotype_variant");
  }
}

TEST_CASE("adjudicator: caller support is joined by chain fingerprint, not native id",
          "[adjudicator][consensus]") {
  // Caller B's "transcript_1" is chainY (1 caller); caller A's "transcript_1" is chainX
  // (2 callers). Adjudicating caller B must see 1, whatever the name says.
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  const auto chainX = chain({{200, 300}, {400, 550}});
  const auto chainY = chain({{200, 550}});
  ConsensusBuilder b;
  b.add("callerA", "transcript_1", chainX);
  b.add("callerB", "transcript_2", chainX);
  b.add("callerB", "transcript_1", chainY);
  const std::string path = scratch("adj_caller_support_collision.tsv");
  write_caller_support_matrix(path, b.build());
  const CallerSupport support = read_caller_support(path);
  std::remove(path.c_str());

  SqantiTable sqanti;
  sqanti.records.push_back(rec("transcript_1", "novel_not_in_catalog", "canonical"));
  sqanti.records.push_back(rec("transcript_2", "novel_not_in_catalog", "canonical"));
  std::map<std::string, IntronChain> chains;  // caller B's isoforms
  chains["transcript_1"] = chainY;
  chains["transcript_2"] = chainX;
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.caller_support = &support;
  const auto out = adjudicate(in, RuleEngine());
  const AdjudicationResult* t1 = find(out, "transcript_1");
  const AdjudicationResult* t2 = find(out, "transcript_2");
  REQUIRE(t1 != nullptr);
  REQUIRE(t2 != nullptr);
  CHECK(t1->evidence.n_callers == 1);
  CHECK(t1->verdict.confidence == ConfidenceClass::kUnresolved);     // no promotion from a name clash
  CHECK(t2->evidence.n_callers == 2);
  CHECK(t2->verdict.confidence == ConfidenceClass::kNovelProbable);
}

TEST_CASE("adjudicator: SQANTI3 NA QC columns are not evaluable (not false)", "[adjudicator][sqanti]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  SqantiTable sqanti;
  SqantiRecord r = rec("isoNA", "novel_not_in_catalog", "NA");
  r.RTS_stage = "NA";
  r.polyA_motif_found = "NA";
  sqanti.records.push_back(r);
  std::map<std::string, IntronChain> chains;
  chains["isoNA"] = chain({{200, 350}, {400, 500}});
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  const auto out = adjudicate(in, RuleEngine());
  const AdjudicationResult* na = find(out, "isoNA");
  REQUIRE(na != nullptr);
  CHECK_FALSE(na->evidence.rts_evaluable);
  CHECK_FALSE(na->evidence.canon_evaluable);
  CHECK_FALSE(na->evidence.polyA_motif_evaluable);
  CHECK(na->verdict.primary_mechanism == Mechanism::kNone);
}

TEST_CASE("adjudicator: the junction table carries the per-junction evidence", "[adjudicator][junctions]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoNNC", "novel_not_in_catalog", "non_canonical"));
  sqanti.records.push_back(rec("isoNNC2", "novel_not_in_catalog", "canonical"));
  std::map<std::string, IntronChain> chains;
  chains["isoNNC"] = chain({{200, 350}, {400, 500}});
  chains["isoNNC2"] = chain({{200, 360}, {400, 500}});
  SjTable sjt;
  sjt.add(sj({400, 500}, 8));
  sjt.add(sj({200, 360}, 10));
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.sj = &sjt;
  std::vector<JunctionRow> jx;
  adjudicate(in, RuleEngine(), &jx);
  REQUIRE(jx.size() == 3);  // {200,350} {200,360} {400,500}, genomic order
  CHECK(jx[0].intron == Junction{200, 350});
  CHECK_FALSE(jx[0].known);
  CHECK(jx[0].sj_evaluable);
  CHECK(jx[0].sj_n_uniq == 0);
  CHECK_FALSE(jx[0].sr_supported);
  CHECK(jx[0].isoforms == std::vector<std::string>{"isoNNC"});
  CHECK(jx[1].intron == Junction{200, 360});
  CHECK(jx[1].sj_n_uniq == 10);
  CHECK(jx[1].sr_supported);
  CHECK(jx[2].intron == Junction{400, 500});
  CHECK(jx[2].known);
  CHECK(jx[2].sj_n_uniq == 8);
  CHECK(jx[2].isoforms.size() == 2);   // shared by both isoforms
  CHECK(jx[2].variant == "NA");        // the variant axis is not evaluated on known junctions
}

TEST_CASE("adjudicator: a motif canonical on the reference only is a non-canonical sign", "[adjudicator][variant]") {
  // "reference" = hapC (both junctions canonical); the sample's haplotype = the GG/GG sequence.
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  auto ref = std::make_shared<FastaFetcher>(tiny("two_jx_hapC.fa"));
  auto hap = std::make_shared<FastaFetcher>(tiny("two_jx_ref.fa"));
  HaplotypeProvider provider(ref, {hap});
  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoD", "novel_not_in_catalog", "canonical"));  // SQANTI3 saw the reference
  std::map<std::string, IntronChain> chains;
  chains["isoD"] = chain_on("chrT2", {{20, 40}, {60, 80}});
  SjTable sjt;  // evaluable, no reads anywhere -> UNSUPPORTED (no known junctions to gauge coverage)
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.sj = &sjt;
  in.haplotype = &provider;

  SECTION("independent provenance: the sign fires, with its own reason suffix") {
    in.haplotype_circular = false;
    std::vector<JunctionRow> jx;
    const auto out = adjudicate(in, RuleEngine(), &jx);
    REQUIRE(out.size() == 1);
    CHECK(out[0].evidence.n_novel_jx_variant_disrupted == 2);
    CHECK(out[0].verdict.primary_mechanism == Mechanism::kNoncanonical);
    CHECK(out[0].verdict.confidence == ConfidenceClass::kArtifact);
    CHECK(out[0].verdict.verdict_reason == "unsupported_noncanonical;haplotype_disrupted");
    REQUIRE(jx.size() == 2);
    CHECK(jx[0].variant == "disrupted");
    CHECK(jx[0].canonical_on.empty());
    CHECK(jx[0].ref_motif == "GT-AG");
    CHECK(jx[0].hap_motifs == std::vector<std::string>{"GG-AG"});
    // Ablating the variant axis removes the sign.
    const auto off = adjudicate(in, RuleEngine().with_axis_disabled("variant"));
    CHECK(off[0].verdict.confidence == ConfidenceClass::kUnconfirmed);
    CHECK(off[0].verdict.verdict_reason == "unsupported_no_sign");
  }
  SECTION("circular-risk provenance: the sign is ignored") {
    in.haplotype_circular = true;
    const auto out = adjudicate(in, RuleEngine());
    CHECK(out[0].evidence.n_novel_jx_variant_disrupted == 2);
    CHECK(out[0].verdict.confidence == ConfidenceClass::kUnconfirmed);
    CHECK(out[0].verdict.verdict_reason == "unsupported_no_sign");
  }
}

TEST_CASE("adjudicator: chain support from the BAM fills the evidence and lifts a held isoform", "[adjudicator][chain]") {
  const Catalog cat = build_catalog_from_gtf(tiny("mini.gtf"));
  BamReader bam(tiny("mini.bam"));
  SqantiTable sqanti;
  sqanti.records.push_back(rec("isoChain", "novel_not_in_catalog", "canonical"));
  std::map<std::string, IntronChain> chains;
  chains["isoChain"] = chain({{2000, 3000}});   // 4 primary reads carry it; no short reads given
  AdjudicateInputs in;
  in.sqanti = &sqanti;
  in.chains = &chains;
  in.catalog = &cat;
  in.bam = &bam;
  const auto out = adjudicate(in, RuleEngine());
  REQUIRE(out.size() == 1);
  CHECK(out[0].evidence.chain_evaluable);
  CHECK(out[0].evidence.chain_n_full_reads == 4);
  CHECK(out[0].verdict.novelty_support == NoveltySupport::kUnknown);
  CHECK(out[0].verdict.confidence == ConfidenceClass::kNovelProbable);
  CHECK(out[0].verdict.verdict_reason == "chain_support");
}
