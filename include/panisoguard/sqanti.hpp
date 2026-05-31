#pragma once

#include <map>
#include <string>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// A SQANTI3 classification row. QC columns are consumed as PRIORS and never
// recomputed. Numeric priors are NaN / INT-sentinel when the source value is
// "NA" or the column is absent (evaluable flags expose which).
struct SqantiRecord {
  std::string isoform;
  std::string chrom;
  Strand strand = Strand::kUnknown;
  std::string structural_category;
  std::string subcategory;
  std::string associated_gene;
  std::string associated_transcript;

  // QC priors
  std::string RTS_stage;       // "TRUE"/"FALSE"/...
  std::string all_canonical;   // "canonical"/"non_canonical"
  double perc_A_downstream_TTS = kNaN;
  int n_indels_junc = kIntNA;
  double dist_to_CAGE_peak = kNaN;
  double dist_to_polyA_site = kNaN;
  std::string filter_result;   // appended SQANTI3 filter column ("Isoform"/"Artifact"), if present

  static constexpr double kNaN = -1e300;  // sentinel; use is_na()
  static constexpr int kIntNA = -2147483647;

  bool is_novel() const {
    return structural_category == "novel_in_catalog" ||
           structural_category == "novel_not_in_catalog";
  }
};

inline bool sqanti_is_na(double v) { return v <= SqantiRecord::kNaN; }

// Result of reading a classification file, including schema diagnostics.
struct SqantiTable {
  std::vector<SqantiRecord> records;
  std::map<std::string, std::size_t> category_counts;  // structural_category -> n
  std::vector<std::string> missing_columns;            // expected-but-absent (not_evaluable)
  std::size_t n_columns = 0;

  std::size_t count(const std::string& category) const {
    auto it = category_counts.find(category);
    return it == category_counts.end() ? 0 : it->second;
  }
};

// Parse a SQANTI3 *_classification.txt (header + tab-separated rows), tolerant of
// column reordering and appended columns. Throws std::runtime_error on open
// failure or if the mandatory `isoform` column is absent.
SqantiTable read_sqanti_classification(const std::string& path);

}  // namespace panisoguard
