#pragma once

#include <string>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// A SQANTI3 classification row: the columns the verdict reads. The QC columns are
// consumed as PRIORS and never recomputed; an absent column or an "NA" value is the
// sentinel below, which the adjudicator turns into a not_evaluable flag.
struct SqantiRecord {
  std::string isoform;
  std::string chrom;
  Strand strand = Strand::kUnknown;
  std::string structural_category;

  // QC priors
  std::string RTS_stage;       // "TRUE"/"FALSE"/...
  std::string all_canonical;   // "canonical"/"non_canonical"
  double perc_A_downstream_TTS = kNaN;
  std::string polyA_motif_found;  // SQANTI3 intra-priming polyA-motif flag, if present

  static constexpr double kNaN = -1e300;  // sentinel; use sqanti_is_na()

  bool is_novel() const {
    return structural_category == "novel_in_catalog" ||
           structural_category == "novel_not_in_catalog";
  }
};

inline bool sqanti_is_na(double v) { return v <= SqantiRecord::kNaN; }
// SQANTI3's textual NA spellings (an absent/NA QC column is not_evaluable, never false).
inline bool sqanti_is_na_token(const std::string& s) {
  return s.empty() || s == "NA" || s == "NaN" || s == "nan" || s == ".";
}

// Result of reading a classification file.
struct SqantiTable {
  std::vector<SqantiRecord> records;
  std::size_t n_columns = 0;
};

// Parse a SQANTI3 *_classification.txt (header + tab-separated rows), tolerant of
// column reordering and appended columns. Throws std::runtime_error on open
// failure or if the mandatory `isoform` column is absent.
SqantiTable read_sqanti_classification(const std::string& path);

}  // namespace panisoguard
