#pragma once

#include <string>
#include <vector>

#include "panisoguard/adjudicator.hpp"

namespace panisoguard {

// Provenance metadata recorded in <prefix>.provenance.log so a run is auditable.
struct RunProvenance {
  std::string tool_version;
  std::string ruleset_version;
  std::string sqanti3_version_target;
  std::string classification_path;
  std::string isoforms_path;
  std::string ref_gtf_path;   // "" if none
  std::string sj_tab_path;    // "" if none
  std::string bam_path;       // "" if none
  std::string caller_support_path;  // "" if none
  bool variant_axis_on = false;     // a haplotype FASTA was supplied
  bool variant_circular = false;    // haplotype provenance is circular-risk
  bool pangenome_axis_on = false;   // a pangenome junction file was supplied
  bool pangenome_circular = false;  // pangenome junction-set provenance is circular-risk
  bool consensus_axis_on = false;   // a caller-support matrix was supplied
  std::string config_path;    // "" if built-in defaults
};

// Write <prefix>.adjudicated.tsv, <prefix>.attribution.jsonl, and
// <prefix>.provenance.log. Throws std::runtime_error on write failure.
void write_adjudication_outputs(const std::string& out_prefix,
                                const std::vector<AdjudicationResult>& results,
                                const RunProvenance& prov);

}  // namespace panisoguard
