// `panisoguard ablate` -- per-axis ablation (the R4 deliverable): re-evaluate the
// assembled evidence with one axis masked at a time and report how many isoforms
// change confidence class, with the class-transition breakdown.

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "panisoguard/analysis.hpp"
#include "run_common.hpp"

namespace {

std::vector<std::string> split_csv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) if (!tok.empty()) out.push_back(tok);
  return out;
}

void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard ablate [adjudicate options] [--axes a,b,..] --out PREFIX\n"
      "\n"
      "Re-runs the rule engine with one evidence axis disabled at a time and reports\n"
      "per-axis confidence-class changes vs the full model.\n"
      "Axes: short_read, mapping, noncanonical, rt_switch, degradation (default: all).\n");
}

}  // namespace

int cmd_ablate(int argc, char** argv) {
  using namespace panisoguard;
  CommonArgs common;
  std::string out_prefix;
  std::vector<std::string> axes;

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (consume_common_arg(a, i, argc, argv, common)) continue;
    if (a == "--out") {
      if (i + 1 >= argc) { std::fprintf(stderr, "ablate: --out requires an argument\n"); return 1; }
      out_prefix = argv[++i];
    } else if (a == "--axes") {
      if (i + 1 >= argc) { std::fprintf(stderr, "ablate: --axes requires an argument\n"); return 1; }
      axes = split_csv(argv[++i]);
    } else if (a == "-h" || a == "--help") { usage(); return 0; }
    else { std::fprintf(stderr, "ablate: unknown option '%s'\n", a.c_str()); usage(); return 1; }
  }

  std::string err;
  if (!common_args_ok(common, err) || out_prefix.empty()) {
    std::fprintf(stderr, "ablate: %s\n", out_prefix.empty() ? "--out is required" : err.c_str());
    usage();
    return 1;
  }
  if (axes.empty()) axes = {"short_read", "mapping", "noncanonical", "rt_switch", "degradation"};

  try {
    const LoadedRun run = load_and_adjudicate(common);

    std::ofstream summary(out_prefix + ".ablation.tsv");
    std::ofstream trans(out_prefix + ".ablation_transitions.tsv");
    if (!summary || !trans) throw std::runtime_error("cannot write ablation outputs under " + out_prefix);
    summary << "axis\tn_changed\tn_total\tfrac_changed\n";
    trans << "axis\ttransition\tcount\n";

    std::fprintf(stderr, "ablation over %zu isoforms:\n", run.results.size());
    for (const auto& axis : axes) {
      const AxisAblation ab = compute_ablation(run.results, run.engine, axis);
      const double frac = ab.n_total ? static_cast<double>(ab.n_changed) / ab.n_total : 0.0;
      summary << ab.axis << '\t' << ab.n_changed << '\t' << ab.n_total << '\t' << frac << '\n';
      for (const auto& kv : ab.transitions) trans << ab.axis << '\t' << kv.first << '\t' << kv.second << '\n';
      std::fprintf(stderr, "  %-14s %8ld / %ld changed (%.3f%%)\n", ab.axis.c_str(), ab.n_changed,
                   ab.n_total, frac * 100.0);
    }
    std::fprintf(stderr, "wrote %s.{ablation,ablation_transitions}.tsv\n", out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ablate: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
