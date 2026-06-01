// PanIsoGuard CLI entry point (samtools-style subcommand dispatch).
//
// Each subcommand (adjudicate / benchmark / ablate / combine) is implemented in its
// own translation unit and dispatched here; `version` reports the build's linked
// htslib and compiled-in capabilities.

#include <clocale>
#include <cstdio>
#include <cstring>
#include <string>

#include <htslib/hts.h>

#include "panisoguard/version.hpp"

// Real subcommands implemented in their own translation units.
int cmd_combine(int argc, char** argv);
int cmd_adjudicate(int argc, char** argv);
int cmd_benchmark(int argc, char** argv);
int cmd_ablate(int argc, char** argv);

namespace {

void print_version() {
  std::printf("panisoguard %s\n", panisoguard::kVersion);
  std::printf("  linked htslib : %s\n", hts_version());
#ifdef HAVE_GBWTGRAPH
  std::printf("  capabilities  : pangenome-gbz=on (in-process, experimental)\n");
#else
  std::printf("  capabilities  : pangenome-gbz=off (GFA / rpvg-text subprocess tier only)\n");
#endif
}

void print_usage() {
  std::fprintf(stderr,
      "panisoguard %s -- caller-agnostic adjudication of long-read novel isoforms\n"
      "\n"
      "Usage: panisoguard <subcommand> [options]\n"
      "\n"
      "Subcommands:\n"
      "  adjudicate        Classify novel isoforms into confidence classes + mechanistic attribution\n"
      "  benchmark         Non-redundancy vs SQANTI3 (McNemar / 2x2 / Jaccard)\n"
      "  ablate            Per-evidence-axis ablation (confidence-class changes)\n"
      "  combine           Integrate multiple callers' isoforms by intron-chain fingerprint\n"
      "  version           Print version, linked htslib, and compiled-in capabilities\n"
      "\n"
      "Options:\n"
      "  -h, --help        Show this help\n"
      "  -v, --version     Show version information\n",
      panisoguard::kVersion);
}

}  // namespace

int main(int argc, char** argv) {
  // Force C numeric locale so float parsing (std::stod) is decimal-point stable
  // regardless of the user's environment locale.
  std::setlocale(LC_NUMERIC, "C");

  if (argc < 2) {
    print_usage();
    return 1;
  }

  const std::string sub = argv[1];

  if (sub == "version" || sub == "--version" || sub == "-v") {
    print_version();
    return 0;
  }
  if (sub == "help" || sub == "--help" || sub == "-h") {
    print_usage();
    return 0;
  }
  if (sub == "adjudicate")      return cmd_adjudicate(argc, argv);
  if (sub == "benchmark")       return cmd_benchmark(argc, argv);
  if (sub == "ablate")          return cmd_ablate(argc, argv);
  if (sub == "combine")         return cmd_combine(argc, argv);

  std::fprintf(stderr, "panisoguard: unknown subcommand '%s'\n\n", sub.c_str());
  print_usage();
  return 1;
}
