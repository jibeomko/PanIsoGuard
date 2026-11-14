#include "catch2/catch.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <string>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/classify.hpp"
#include "test_util.hpp"

using namespace panisoguard;

namespace {
// chrC, plus strand: A has introns {100,200},{300,400},{500,600}; B has {100,700},{800,900}.
// M on the minus strand has {2100,2200}. Nothing is annotated past 2300.
std::string write_gtf() {
  const std::string path = scratch("classify_ref.gtf");
  std::ofstream g(path);
  auto tx = [&](const char* id, const char* strand, std::initializer_list<std::pair<int, int>> exons) {
    for (const auto& e : exons)
      g << "chrC\tt\texon\t" << e.first << '\t' << e.second << "\t.\t" << strand << "\t.\tgene_id \"g_" << id << "\"; transcript_id \"" << id << "\";\n";
  };
  tx("A", "+", {{1, 100}, {201, 300}, {401, 500}, {601, 700}});
  tx("B", "+", {{1, 100}, {701, 800}, {901, 1000}});
  tx("M", "-", {{2001, 2100}, {2201, 2300}});
  tx("S", "+", {{3001, 3400}});                                   // a mono-exonic transcript
  return path;
}
// chrC, 6000 bases: GT at every annotated intron start and AG before every annotated intron
// end; the 20 bases after position 700 are A and the 20 after 1000 are C.
std::string write_fasta() {
  const std::string path = scratch("classify_ref.fa");
  std::string seq(6000, 'C');
  for (int s : {100, 300, 500, 700, 800}) { seq[s] = 'G'; seq[s + 1] = 'T'; }
  for (int e : {200, 400, 600, 700, 900}) { seq[e - 2] = 'A'; seq[e - 1] = 'G'; }
  for (int i = 700; i < 720; ++i) seq[i] = 'A';
  std::ofstream f(path);
  f << ">chrC\n";
  for (std::size_t i = 0; i < seq.size(); i += 60) f << seq.substr(i, 60) << '\n';
  f.close();
  std::remove((path + ".fai").c_str());
  return path;
}
CallerIsoform iso(const char* id, Strand strand, std::vector<Junction> introns, int64_t start, int64_t end) {
  CallerIsoform c; c.id = id; c.chain.chrom = "chrC"; c.chain.strand = strand; c.chain.introns = std::move(introns); c.start = start; c.end = end;
  return c;
}
}  // namespace

TEST_CASE("classify_isoforms gives SQANTI3's structural categories from the chains", "[classify]") {
  const std::string gtf = write_gtf();
  const Catalog cat = build_catalog_from_gtf(gtf);
  std::vector<CallerIsoform> isos = {
      iso("fsm", Strand::kPlus, {{100, 200}, {300, 400}, {500, 600}}, 0, 700),
      iso("ism", Strand::kPlus, {{300, 400}, {500, 600}}, 250, 700),
      iso("nic", Strand::kPlus, {{300, 400}, {800, 900}}, 250, 1000),
      iso("nic_sites", Strand::kPlus, {{100, 400}}, 0, 500),           // known donor and acceptor, new junction
      iso("nic_retained", Strand::kPlus, {{100, 200}, {300, 400}}, 0, 700),   // A's first two introns, span retains the third
      iso("nnc", Strand::kPlus, {{100, 200}, {350, 400}}, 0, 500),
      iso("anti", Strand::kPlus, {{2100, 2200}}, 2000, 2300),
      iso("inter", Strand::kPlus, {{5000, 5100}}, 4900, 5200),
      iso("mono_fsm", Strand::kPlus, {}, 3050, 3300),                  // inside S's single exon
      iso("mono_ism", Strand::kPlus, {}, 210, 290),                    // inside A's second exon
      iso("mono_genic", Strand::kPlus, {}, 250, 350),                  // crosses A's exon end
      iso("mono_intron", Strand::kPlus, {}, 110, 190),                 // inside A's first intron
      iso("mono_inter", Strand::kPlus, {}, 5000, 5100),
  };
  const SqantiTable t = classify_isoforms(isos, cat, nullptr);
  REQUIRE(t.records.size() == 13);
  std::map<std::string, std::string> cat_of;
  for (const auto& r : t.records) cat_of[r.isoform] = r.structural_category;
  CHECK(cat_of["fsm"] == "full-splice_match");
  CHECK(cat_of["ism"] == "incomplete-splice_match");
  CHECK(cat_of["nic"] == "novel_in_catalog");
  CHECK(cat_of["nic_sites"] == "novel_in_catalog");
  CHECK(cat_of["nic_retained"] == "novel_in_catalog");
  CHECK(cat_of["nnc"] == "novel_not_in_catalog");
  CHECK(cat_of["anti"] == "antisense");
  CHECK(cat_of["inter"] == "intergenic");
  CHECK(cat_of["mono_fsm"] == "full-splice_match");
  CHECK(cat_of["mono_ism"] == "incomplete-splice_match");
  CHECK(cat_of["mono_genic"] == "genic");
  CHECK(cat_of["mono_intron"] == "genic_intron");
  CHECK(cat_of["mono_inter"] == "intergenic");
  // without the genome the QC priors are not evaluable, never false
  CHECK(sqanti_is_na_token(t.records[0].all_canonical));
  CHECK(sqanti_is_na(t.records[0].perc_A_downstream_TTS));
  CHECK(sqanti_is_na_token(t.records[0].RTS_stage));
  CHECK(t.records[0].chrom == "chrC");
  CHECK(t.records[0].strand == Strand::kPlus);
  std::remove(gtf.c_str());
}

TEST_CASE("classify_isoforms reads splice motifs and the 3' A content from the genome", "[classify]") {
  const std::string gtf = write_gtf(), fa = write_fasta();
  const Catalog cat = build_catalog_from_gtf(gtf);
  FastaFetcher ref(fa);
  std::vector<CallerIsoform> isos = {
      iso("fsm", Strand::kPlus, {{100, 200}, {300, 400}, {500, 600}}, 0, 700),      // GT-AG everywhere, A-rich 3' end
      iso("nnc", Strand::kPlus, {{100, 200}, {350, 400}}, 0, 500),                  // {350,400}: CC-AG
      iso("bend", Strand::kPlus, {{100, 700}, {800, 900}}, 0, 1000),                // ends before 20 C bases
      iso("minus", Strand::kMinus, {{2100, 2200}}, 1000, 2300),                     // 3' end at 1000: bases 980..1000 are C, no T
      iso("mono", Strand::kPlus, {}, 5000, 5100),
  };
  const SqantiTable t = classify_isoforms(isos, cat, &ref);
  std::map<std::string, SqantiRecord> r;
  for (const auto& rec : t.records) r[rec.isoform] = rec;
  CHECK(r["fsm"].all_canonical == "canonical");
  CHECK(r["fsm"].perc_A_downstream_TTS == Approx(100.0));
  CHECK(r["nnc"].all_canonical == "non_canonical");
  CHECK(r["bend"].all_canonical == "canonical");
  CHECK(r["bend"].perc_A_downstream_TTS == Approx(0.0));
  CHECK(r["minus"].perc_A_downstream_TTS == Approx(0.0));
  CHECK(sqanti_is_na_token(r["mono"].all_canonical));
  // the rows feed the adjudicator like SQANTI3's: FSM is KNOWN, mono-exon is not adjudicated
  std::map<std::string, IntronChain> chains;
  for (const auto& i : isos) chains[i.id] = i.chain;
  AdjudicateInputs in; in.sqanti = &t; in.chains = &chains; in.catalog = &cat;
  const auto out = adjudicate(in, RuleEngine());
  for (const auto& o : out) {
    if (o.isoform_id == "fsm") CHECK(o.verdict.confidence == ConfidenceClass::kKnown);
    if (o.isoform_id == "mono") CHECK(o.verdict.verdict == "not_adjudicated");   // intergenic
    if (o.isoform_id == "nnc") CHECK(o.evidence.n_novel_junctions == 1);
  }
  std::remove(gtf.c_str()); std::remove(fa.c_str()); std::remove((fa + ".fai").c_str());
}
