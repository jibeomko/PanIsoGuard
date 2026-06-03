// `panisoguard adjudicate` -- Tier-0+ end-to-end: consume a SQANTI3 classification
// + caller isoform chains (+ optional reference GTF, STAR SJ.tab, and BAM), assign
// a confidence class with mechanistic attribution, and write three outputs.

#include <cstdio>
#include <map>
#include <string>

#include "panisoguard/result_writer.hpp"
#include "panisoguard/version.hpp"
#include "run_common.hpp"

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
      "  --reference PATH        genome FASTA (for the variant axis; also decodes a CRAM --bam)\n"
      "  --reference-haplotype PATH   personalized haplotype FASTA for the variant axis (repeatable)\n"
      "  --haplotype-provenance X     rna_derived|wgs|external|unknown (default unknown=circular-risk)\n"
      "  --pangenome-junctions PATH   graph-supported splice junctions (pangenome reference-bias axis; optional)\n"
      "  --pangenome-provenance X     population|external|sample_derived|unknown (default unknown=circular-risk)\n"
      "  --caller-support PATH        `panisoguard combine` matrix for the multi-caller consensus axis (optional)\n"
      "  --config PATH           rules TOML (optional; built-in defaults otherwise)\n"
      "  --out-prefix PREFIX     writes .adjudicated.tsv / .attribution.jsonl / .provenance.log\n"
      "\n"
      "For a visual PDF summary of the output:\n"
      "  pip install ./python && panisoguard-report --prefix PREFIX   (-> PREFIX.report.pdf)\n");
}

}  // namespace

int cmd_adjudicate(int argc, char** argv) {
  using namespace panisoguard;
  CommonArgs common;
  std::string out_prefix;

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (consume_common_arg(a, i, argc, argv, common)) continue;
    if (a == "--out-prefix") {
      if (i + 1 >= argc) { std::fprintf(stderr, "adjudicate: --out-prefix requires an argument\n"); return 1; }
      out_prefix = argv[++i];
    } else if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else {
      std::fprintf(stderr, "adjudicate: unknown option '%s'\n", a.c_str());
      usage();
      return 1;
    }
  }

  std::string err;
  if (!common_args_ok(common, err) || out_prefix.empty()) {
    std::fprintf(stderr, "adjudicate: %s\n", out_prefix.empty() ? "--out-prefix is required" : err.c_str());
    usage();
    return 1;
  }

  try {
    const LoadedRun run = load_and_adjudicate(common);

    RunProvenance prov;
    prov.tool_version = kVersion;
    prov.ruleset_version = run.engine.config().ruleset_version;
    prov.sqanti3_version_target = run.engine.config().sqanti3_version_target;
    prov.classification_path = common.classification;
    prov.isoforms_path = common.isoforms_bed.empty() ? common.isoforms_gtf : common.isoforms_bed;
    prov.ref_gtf_path = common.ref_gtf;
    prov.sj_tab_path = common.sj_tab;
    prov.bam_path = common.bam_path;
    prov.variant_axis_on = !common.reference_haplotypes.empty();
    prov.variant_circular = haplotype_provenance_is_circular(common.haplotype_provenance);
    prov.pangenome_axis_on = !common.pangenome_junctions.empty();
    prov.pangenome_circular = pangenome_provenance_is_circular(common.pangenome_provenance);
    prov.config_path = common.config;
    write_adjudication_outputs(out_prefix, run.results, prov);

    std::map<std::string, std::size_t> counts;
    for (const auto& r : run.results) ++counts[to_string(r.verdict.confidence)];
    std::fprintf(stderr, "adjudicated %zu isoforms:\n", run.results.size());
    for (const auto& kv : counts) std::fprintf(stderr, "  %-28s %zu\n", kv.first.c_str(), kv.second);
    std::fprintf(stderr, "wrote %s.{adjudicated.tsv,attribution.jsonl,provenance.log}\n", out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "adjudicate: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
