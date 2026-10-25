// PanIsoGuard CLI entry point (samtools-style subcommand dispatch).
//
// Each subcommand (adjudicate / combine / decoys / example) is implemented in its own
// translation unit and dispatched here; `version` reports the build's linked htslib.

#include <clocale>
#include <cstdio>
#include <cstring>
#include <string>

#include <htslib/hts.h>

#include "panisoguard/version.hpp"

// Real subcommands implemented in their own translation units.
int cmd_combine(int argc, char** argv);
int cmd_adjudicate(int argc, char** argv);
int cmd_example(int argc, char** argv);
int cmd_decoys(int argc, char** argv);

namespace {

void print_version() {
  std::printf("panisoguard %s\n", panisoguard::kVersion);
  std::printf("  linked htslib : %s\n", hts_version());
}

void print_usage() {
  std::fprintf(stderr,
      "panisoguard %s -- caller-agnostic adjudication of long-read novel isoforms\n"
      "\n"
      "Usage: panisoguard <subcommand> [options]\n"
      "\n"
      "Subcommands:\n"
      "  adjudicate        Classify novel isoforms into confidence classes + mechanistic attribution\n"
      "  combine           Integrate multiple callers' isoforms by intron-chain fingerprint\n"
      "  decoys            Per-dataset false-lift check of the whole-chain read support (recombined known chains)\n"
      "  example           Copy a bundled example (tiny, multi_caller, worked_locus) into a folder\n"
      "  version           Print version and linked htslib\n"
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
  if (sub == "combine")         return cmd_combine(argc, argv);
  if (sub == "decoys")          return cmd_decoys(argc, argv);
  if (sub == "example")         return cmd_example(argc, argv);

  std::fprintf(stderr, "panisoguard: unknown subcommand '%s'\n\n", sub.c_str());
  print_usage();
  return 1;
}
