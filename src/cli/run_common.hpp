#pragma once

#include <string>
#include <vector>

#include "panisoguard/adjudicator.hpp"
#include "panisoguard/rules.hpp"
#include "panisoguard/sqanti.hpp"

namespace panisoguard {

// Inputs shared by the adjudicate / benchmark / ablate subcommands.
struct CommonArgs {
  std::string classification;
  std::string isoforms_bed;
  std::string isoforms_gtf;
  std::string ref_gtf;
  std::string sj_tab;
  std::string bam_path;
  std::string reference;
  std::string config;
};

// If argv[i] is a shared flag, store it (advancing i past its value) and return
// true; otherwise return false so the caller can handle subcommand-specific flags.
bool consume_common_arg(const std::string& a, int& i, int argc, char** argv, CommonArgs& c);

// Validate the minimum required inputs; on failure sets `err` and returns false.
bool common_args_ok(const CommonArgs& c, std::string& err);

// One adjudication run plus the SQANTI table (needed for benchmark) and the
// engine (needed for ablation).
struct LoadedRun {
  SqantiTable sqanti;
  std::vector<AdjudicationResult> results;
  RuleEngine engine;
};

// Read every supplied input and run adjudication. Throws std::exception on I/O
// errors. Prints brief progress to stderr.
LoadedRun load_and_adjudicate(const CommonArgs& c);

}  // namespace panisoguard
