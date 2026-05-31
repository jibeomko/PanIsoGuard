#include "run_common.hpp"

#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

#include "panisoguard/bed12.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/sj_tab.hpp"
#include "panisoguard/variant_motif.hpp"

namespace panisoguard {

bool consume_common_arg(const std::string& a, int& i, int argc, char** argv, CommonArgs& c) {
  auto next = [&](const char* name) -> std::string {
    if (i + 1 >= argc) {
      std::fprintf(stderr, "panisoguard: %s requires an argument\n", name);
      std::exit(1);
    }
    return argv[++i];
  };
  if (a == "--classification")     { c.classification = next("--classification"); return true; }
  if (a == "--isoforms-bed")       { c.isoforms_bed = next("--isoforms-bed"); return true; }
  if (a == "--isoforms-gtf")       { c.isoforms_gtf = next("--isoforms-gtf"); return true; }
  if (a == "--ref-gtf")            { c.ref_gtf = next("--ref-gtf"); return true; }
  if (a == "--sj-tab")             { c.sj_tab = next("--sj-tab"); return true; }
  if (a == "--bam")                { c.bam_path = next("--bam"); return true; }
  if (a == "--reference")          { c.reference = next("--reference"); return true; }
  if (a == "--reference-haplotype"){ c.reference_haplotypes.push_back(next("--reference-haplotype")); return true; }
  if (a == "--haplotype-provenance"){ c.haplotype_provenance = next("--haplotype-provenance"); return true; }
  if (a == "--config")             { c.config = next("--config"); return true; }
  return false;
}

bool common_args_ok(const CommonArgs& c, std::string& err) {
  if (c.classification.empty()) { err = "--classification is required"; return false; }
  if (c.isoforms_bed.empty() && c.isoforms_gtf.empty()) {
    err = "one of --isoforms-bed / --isoforms-gtf is required";
    return false;
  }
  return true;
}

LoadedRun load_and_adjudicate(const CommonArgs& c) {
  LoadedRun run;
  run.sqanti = read_sqanti_classification(c.classification);
  std::fprintf(stderr, "read %zu SQANTI records (%zu cols)\n", run.sqanti.records.size(),
               run.sqanti.n_columns);

  std::map<std::string, IntronChain> chains;
  if (!c.isoforms_bed.empty()) {
    for (auto& r : read_bed12(c.isoforms_bed)) chains.emplace(r.name, std::move(r.chain));
  } else {
    for (auto& t : read_gtf_transcripts(c.isoforms_gtf)) chains.emplace(t.id, std::move(t.chain));
  }
  std::fprintf(stderr, "read %zu caller isoform chains\n", chains.size());

  Catalog catalog;
  const Catalog* catalog_ptr = nullptr;
  if (!c.ref_gtf.empty()) {
    std::fprintf(stderr, "building reference catalog...\n");
    catalog = build_catalog_from_gtf(c.ref_gtf);
    catalog_ptr = &catalog;
  }

  SjTable sj;
  const SjTable* sj_ptr = nullptr;
  if (!c.sj_tab.empty()) {
    sj = read_sj_tab(c.sj_tab);
    sj_ptr = &sj;
    std::fprintf(stderr, "read %zu short-read junctions\n", sj.size());
  }

  std::unique_ptr<BamReader> bam;
  if (!c.bam_path.empty()) {
    bam = std::make_unique<BamReader>(c.bam_path, c.reference);
    std::fprintf(stderr, "opened BAM for mapping axis: %s\n", c.bam_path.c_str());
  }

  std::shared_ptr<FastaFetcher> ref_fa;
  std::vector<std::shared_ptr<FastaFetcher>> hap_fas;
  std::unique_ptr<HaplotypeProvider> haplo;
  if (!c.reference_haplotypes.empty()) {
    if (c.reference.empty()) {
      throw std::runtime_error("--reference (genome FASTA) is required with --reference-haplotype");
    }
    ref_fa = std::make_shared<FastaFetcher>(c.reference);
    for (const auto& h : c.reference_haplotypes) hap_fas.push_back(std::make_shared<FastaFetcher>(h));
    haplo = std::make_unique<HaplotypeProvider>(ref_fa, hap_fas);
    std::fprintf(stderr, "variant axis: reference + %zu haplotype FASTA(s), provenance=%s\n",
                 hap_fas.size(), c.haplotype_provenance.c_str());
  }

  run.engine = c.config.empty() ? RuleEngine() : RuleEngine::from_toml(c.config);

  AdjudicateInputs in;
  in.sqanti = &run.sqanti;
  in.chains = &chains;
  in.catalog = catalog_ptr;
  in.sj = sj_ptr;
  in.bam = bam.get();
  in.haplotype = haplo.get();
  in.haplotype_circular = haplotype_provenance_is_circular(c.haplotype_provenance);
  run.results = adjudicate(in, run.engine);
  return run;
}

}  // namespace panisoguard
