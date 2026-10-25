// `panisoguard adjudicate` -- Tier-0+ end-to-end: consume a SQANTI3 classification
// + caller isoform chains (+ optional reference GTF, STAR SJ.tab, and BAM), assign
// a confidence class with mechanistic attribution, and write three outputs.

#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

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
      "  --validation-bam PATH   a second, independent library: whole-chain read support is counted\n"
      "                          there instead of in --bam (optional; evidence is then independent)\n"
      "  --reference PATH        genome FASTA (for the variant axis; also decodes a CRAM --bam)\n"
      "  --reference-haplotype PATH   personalized haplotype FASTA for the variant axis (repeatable)\n"
      "  --haplotype-provenance X     rna_derived|wgs|external|unknown (default unknown=circular-risk)\n"
      "  --platform hifi|ont          BAM-axis threshold preset (default hifi; ont turns the\n"
      "                               HiFi-calibrated indel-near gate off)\n"
      "  --pangenome-junctions PATH   graph-supported splice junctions (pangenome reference-bias axis; optional)\n"
      "  --pangenome-provenance X     population|external|sample_derived|unknown (default unknown=circular-risk)\n"
      "  --caller-support PATH        `panisoguard combine` matrix for the multi-caller consensus axis (optional)\n"
      "  --write-kept PATH            also write the input isoform GTF/BED12 filtered to the kept verdicts\n"
      "  --keep-verdicts LIST         verdicts --write-kept keeps (comma-separated; default\n"
      "                               known,known_fragment,confirmed,reference_bias)\n"
      "  --config PATH           rules TOML (optional; built-in defaults otherwise)\n"
      "  --out-prefix PREFIX     writes .adjudicated.tsv / .junctions.tsv / .attribution.jsonl / .provenance.log\n");
}

// Comma-separated verdict names -> set; returns the first unknown name, or "" when all are valid.
std::string parse_verdict_list(const std::string& list, std::set<std::string>& out) {
  static const std::set<std::string> known = {"known", "known_fragment", "confirmed", "reference_bias",
                                              "unconfirmed", "artifact", "not_adjudicated"};
  std::stringstream ss(list);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) continue;
    if (known.count(tok) == 0) return tok;
    out.insert(tok);
  }
  return "";
}

}  // namespace

int cmd_adjudicate(int argc, char** argv) {
  using namespace panisoguard;
  CommonArgs common;
  std::string out_prefix, write_kept;
  std::string keep_verdicts = "known,known_fragment,confirmed,reference_bias";

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (consume_common_arg(a, i, argc, argv, common)) continue;
    if (a == "--out-prefix") {
      if (i + 1 >= argc) { std::fprintf(stderr, "adjudicate: --out-prefix requires an argument\n"); return 1; }
      out_prefix = argv[++i];
    } else if (a == "--write-kept") {
      if (i + 1 >= argc) { std::fprintf(stderr, "adjudicate: --write-kept requires an argument\n"); return 1; }
      write_kept = argv[++i];
    } else if (a == "--keep-verdicts") {
      if (i + 1 >= argc) { std::fprintf(stderr, "adjudicate: --keep-verdicts requires an argument\n"); return 1; }
      keep_verdicts = argv[++i];
    } else if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else {
      std::fprintf(stderr, "adjudicate: unknown option '%s'\n", a.c_str());
      usage();
      return 1;
    }
  }

  std::set<std::string> keep;
  if (const std::string bad = parse_verdict_list(keep_verdicts, keep); !bad.empty()) {
    std::fprintf(stderr, "adjudicate: --keep-verdicts: unknown verdict '%s'\n", bad.c_str());
    return 1;
  }

  std::string err;
  if (!common_args_ok(common, err) || out_prefix.empty()) {
    std::fprintf(stderr, "adjudicate: %s\n", out_prefix.empty() ? "--out-prefix is required" : err.c_str());
    usage();
    return 1;
  }

  try {
    const LoadedRun run = load_and_adjudicate(common);
    RunProvenance prov = provenance_for(common, run);
    prov.kept_path = write_kept;
    prov.kept_verdicts = keep_verdicts;
    write_adjudication_outputs(out_prefix, run.results, prov);
    write_junctions_tsv(out_prefix + ".junctions.tsv", run.junctions);
    if (!write_kept.empty()) {
      std::unordered_set<std::string> ids;
      for (const auto& r : run.results) if (keep.count(r.verdict.verdict)) ids.insert(r.isoform_id);
      const std::size_t n = write_kept_isoforms(prov.isoforms_path, !common.isoforms_bed.empty(), write_kept, ids);
      std::fprintf(stderr, "kept %zu of %zu isoforms (%s) -> %s\n", n, run.results.size(), keep_verdicts.c_str(),
                   write_kept.c_str());
    }

    std::map<std::string, std::size_t> counts;
    for (const auto& r : run.results) ++counts[to_string(r.verdict.confidence)];
    std::fprintf(stderr, "adjudicated %zu isoforms:\n", run.results.size());
    for (const auto& kv : counts) std::fprintf(stderr, "  %-28s %zu\n", kv.first.c_str(), kv.second);
    std::fprintf(stderr, "wrote %s.{adjudicated.tsv,attribution.jsonl,junctions.tsv,provenance.log}\n", out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "adjudicate: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
