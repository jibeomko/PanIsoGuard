#pragma once

#include <string>
#include <unordered_set>
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
  std::string validation_bam_path;  // "" if none (independent library for whole-chain reads)
  std::string caller_support_path;  // "" if none
  std::string reference_path;       // "" if none (genome FASTA)
  std::vector<std::string> reference_haplotype_paths;  // personalized FASTAs (variant axis)
  std::string pangenome_path;       // "" if none
  std::string kept_path;            // "" if none (--write-kept: the isoform file filtered to kept verdicts)
  std::string kept_verdicts;        // the --keep-verdicts list that selected them
  bool variant_axis_on = false;     // a haplotype FASTA was supplied
  bool variant_circular = false;    // haplotype provenance is circular-risk
  bool pangenome_axis_on = false;   // a pangenome junction file was supplied
  bool pangenome_circular = false;  // pangenome junction-set provenance is circular-risk
  bool consensus_axis_on = false;   // a caller-support matrix was supplied
  std::string config_path;    // "" if built-in defaults
  std::string thresholds;     // effective rule thresholds, "key=value ..." (so a run is
                              // reproducible even if the config file changes later)
};

// Copy the caller isoform file (GTF, or BED12 when is_bed) to out_path keeping only the
// records whose transcript_id / name is in `keep`; GTF comment lines are kept, gene lines
// (no transcript_id) are dropped. Returns the number of distinct ids written.
std::size_t write_kept_isoforms(const std::string& in_path, bool is_bed, const std::string& out_path,
                                const std::unordered_set<std::string>& keep);

// Write <prefix>.junctions.tsv: one row per junction of the adjudicated isoforms with the
// per-junction evidence (intron coordinates 1-based inclusive, as in SJ.tab and GTF).
void write_junctions_tsv(const std::string& path, const std::vector<JunctionRow>& rows);

// Write <prefix>.adjudicated.tsv, <prefix>.attribution.jsonl, and
// <prefix>.provenance.log. Throws std::runtime_error on write failure.
void write_adjudication_outputs(const std::string& out_prefix,
                                const std::vector<AdjudicationResult>& results,
                                const RunProvenance& prov);

}  // namespace panisoguard
