// `panisoguard benchmark` -- non-redundancy of PanIsoGuard vs SQANTI3 on novel
// isoforms (the R1 deliverable): 2x2 keep/flag contingency, McNemar, Jaccard, and
// the per-isoform mechanistic reason for calls PanIsoGuard flags that SQANTI3 keeps.

#include <cstdio>
#include <fstream>
#include <string>

#include "panisoguard/analysis.hpp"
#include "run_common.hpp"

namespace {
void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard benchmark [adjudicate options] --out DIR_OR_PREFIX\n"
      "\n"
      "Runs adjudication, then compares PanIsoGuard (ARTIFACT) vs SQANTI3\n"
      "(filter_result==Artifact) on novel isoforms. Requires --classification with a\n"
      "SQANTI3 filter column, caller isoforms, and (recommended) --ref-gtf/--sj-tab/--bam.\n");
}
}  // namespace

int cmd_benchmark(int argc, char** argv) {
  using namespace panisoguard;
  CommonArgs common;
  std::string out_prefix;

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (consume_common_arg(a, i, argc, argv, common)) continue;
    if (a == "--out") {
      if (i + 1 >= argc) { std::fprintf(stderr, "benchmark: --out requires an argument\n"); return 1; }
      out_prefix = argv[++i];
    } else if (a == "-h" || a == "--help") { usage(); return 0; }
    else { std::fprintf(stderr, "benchmark: unknown option '%s'\n", a.c_str()); usage(); return 1; }
  }

  std::string err;
  if (!common_args_ok(common, err) || out_prefix.empty()) {
    std::fprintf(stderr, "benchmark: %s\n", out_prefix.empty() ? "--out is required" : err.c_str());
    usage();
    return 1;
  }

  try {
    const LoadedRun run = load_and_adjudicate(common);
    const NonRedundancy nr = compute_nonredundancy(run.results, run.sqanti);

    if (nr.n_sqanti_filter_available == 0) {
      std::fprintf(stderr, "benchmark: no SQANTI3 filter_result column found; cannot compare.\n");
      return 1;
    }

    // Summary TSV.
    {
      std::ofstream out(out_prefix + ".nonredundancy.tsv");
      if (!out) throw std::runtime_error("cannot write " + out_prefix + ".nonredundancy.tsv");
      out << "metric\tvalue\n";
      out << "n_novel\t" << nr.n_novel << '\n';
      out << "n_with_sqanti_filter\t" << nr.n_sqanti_filter_available << '\n';
      out << "both_flag(agree_artifact)\t" << nr.both_flag << '\n';
      out << "both_keep(agree_kept)\t" << nr.both_keep << '\n';
      out << "panisoguard_only_flag\t" << nr.pig_only_flag << '\n';
      out << "sqanti_only_flag\t" << nr.sqanti_only_flag << '\n';
      out << "jaccard_flagged\t" << nr.jaccard << '\n';
      out << "mcnemar_chi2\t" << nr.mcnemar_chi2 << '\n';
      out << "mcnemar_p\t" << nr.mcnemar_p << '\n';
    }
    // The non-redundant catch set, with mechanistic reason.
    {
      std::ofstream out(out_prefix + ".panisoguard_only_flagged.tsv");
      if (!out) throw std::runtime_error("cannot write " + out_prefix + ".panisoguard_only_flagged.tsv");
      out << "isoform_id\tprimary_mechanism\n";
      for (const auto& e : nr.pig_only_examples) out << e.first << '\t' << e.second << '\n';
    }

    std::fprintf(stderr,
        "non-redundancy on %ld novel isoforms (%ld with SQANTI filter):\n"
        "  both flag artifact : %ld\n  both keep          : %ld\n"
        "  PanIsoGuard-only   : %ld   <- non-redundant catch (SQANTI keeps)\n"
        "  SQANTI-only        : %ld\n  Jaccard(flagged)   : %.3f\n"
        "  McNemar chi2=%.3f  p=%.3g\n  wrote %s.{nonredundancy,panisoguard_only_flagged}.tsv\n",
        nr.n_novel, nr.n_sqanti_filter_available, nr.both_flag, nr.both_keep, nr.pig_only_flag,
        nr.sqanti_only_flag, nr.jaccard, nr.mcnemar_chi2, nr.mcnemar_p, out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "benchmark: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
