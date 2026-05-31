// `panisoguard adjudicate` -- Tier-0 end-to-end: consume a SQANTI3 classification
// + caller isoform chains (+ optional reference GTF and short-read SJ.tab), assign
// a confidence class with mechanistic attribution, and write three outputs.

#include <cstdio>
#include <map>
#include <memory>
#include <string>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/bed12.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/result_writer.hpp"
#include "panisoguard/rules.hpp"
#include "panisoguard/sj_tab.hpp"
#include "panisoguard/sqanti.hpp"
#include "panisoguard/version.hpp"

namespace {

void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard adjudicate [options] --out-prefix PREFIX\n"
      "\n"
      "  --classification PATH   SQANTI3 *_classification.txt (required)\n"
      "  --isoforms-bed PATH     caller isoform BED12 (provide this or --isoforms-gtf)\n"
      "  --isoforms-gtf PATH     caller isoform GTF\n"
      "  --ref-gtf PATH          reference GTF; enables novel-junction classification (recommended)\n"
      "  --sj-tab PATH           STAR SJ.tab for short-read corroboration (optional)\n"
      "  --bam PATH              indexed BAM/CRAM for the read-level mapping axis (optional)\n"
      "  --reference PATH        reference FASTA (required only to decode a CRAM --bam)\n"
      "  --config PATH           rules TOML (optional; built-in defaults otherwise)\n"
      "  --out-prefix PREFIX     output prefix (writes .adjudicated.tsv/.attribution.jsonl/.provenance.log)\n");
}

}  // namespace

int cmd_adjudicate(int argc, char** argv) {
  std::string classification, isoforms_bed, isoforms_gtf, ref_gtf, sj_tab, bam_path, reference, config, out_prefix;

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "panisoguard adjudicate: %s requires an argument\n", name);
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "--classification")     classification = next("--classification");
    else if (a == "--isoforms-bed")  isoforms_bed = next("--isoforms-bed");
    else if (a == "--isoforms-gtf")  isoforms_gtf = next("--isoforms-gtf");
    else if (a == "--ref-gtf")       ref_gtf = next("--ref-gtf");
    else if (a == "--sj-tab")        sj_tab = next("--sj-tab");
    else if (a == "--bam")           bam_path = next("--bam");
    else if (a == "--reference")     reference = next("--reference");
    else if (a == "--config")        config = next("--config");
    else if (a == "--out-prefix")    out_prefix = next("--out-prefix");
    else if (a == "-h" || a == "--help") { usage(); return 0; }
    else { std::fprintf(stderr, "panisoguard adjudicate: unknown option '%s'\n", a.c_str()); usage(); return 1; }
  }

  if (classification.empty() || out_prefix.empty() ||
      (isoforms_bed.empty() && isoforms_gtf.empty())) {
    std::fprintf(stderr, "panisoguard adjudicate: --classification, one of --isoforms-bed/--isoforms-gtf, and --out-prefix are required\n");
    usage();
    return 1;
  }

  using namespace panisoguard;
  try {
    const SqantiTable sqanti = read_sqanti_classification(classification);
    std::fprintf(stderr, "read %zu SQANTI records (%zu cols)\n", sqanti.records.size(), sqanti.n_columns);

    // Caller chains keyed by isoform id (== SQANTI isoform id).
    std::map<std::string, IntronChain> chains;
    if (!isoforms_bed.empty()) {
      for (auto& r : read_bed12(isoforms_bed)) chains.emplace(r.name, std::move(r.chain));
    } else {
      for (auto& t : read_gtf_transcripts(isoforms_gtf)) chains.emplace(t.id, std::move(t.chain));
    }
    std::fprintf(stderr, "read %zu caller isoform chains\n", chains.size());

    Catalog catalog;
    const Catalog* catalog_ptr = nullptr;
    if (!ref_gtf.empty()) {
      std::fprintf(stderr, "building reference catalog...\n");
      catalog = build_catalog_from_gtf(ref_gtf);
      catalog_ptr = &catalog;
    }

    SjTable sj;
    const SjTable* sj_ptr = nullptr;
    if (!sj_tab.empty()) {
      sj = read_sj_tab(sj_tab);
      sj_ptr = &sj;
      std::fprintf(stderr, "read %zu short-read junctions\n", sj.size());
    }

    std::unique_ptr<BamReader> bam;
    if (!bam_path.empty()) {
      bam = std::make_unique<BamReader>(bam_path, reference);
      std::fprintf(stderr, "opened BAM for mapping axis: %s\n", bam_path.c_str());
    }

    const RuleEngine engine = config.empty() ? RuleEngine() : RuleEngine::from_toml(config);

    AdjudicateInputs in;
    in.sqanti = &sqanti;
    in.chains = &chains;
    in.catalog = catalog_ptr;
    in.sj = sj_ptr;
    in.bam = bam.get();
    const auto results = adjudicate(in, engine);

    RunProvenance prov;
    prov.tool_version = kVersion;
    prov.ruleset_version = engine.config().ruleset_version;
    prov.sqanti3_version_target = engine.config().sqanti3_version_target;
    prov.classification_path = classification;
    prov.isoforms_path = isoforms_bed.empty() ? isoforms_gtf : isoforms_bed;
    prov.ref_gtf_path = ref_gtf;
    prov.sj_tab_path = sj_tab;
    prov.bam_path = bam_path;
    prov.config_path = config;
    write_adjudication_outputs(out_prefix, results, prov);

    // Summary.
    std::map<std::string, std::size_t> counts;
    for (const auto& r : results) ++counts[to_string(r.verdict.confidence)];
    std::fprintf(stderr, "adjudicated %zu isoforms:\n", results.size());
    for (const auto& kv : counts) std::fprintf(stderr, "  %-28s %zu\n", kv.first.c_str(), kv.second);
    std::fprintf(stderr, "wrote %s.{adjudicated.tsv,attribution.jsonl,provenance.log}\n", out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "panisoguard adjudicate: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
