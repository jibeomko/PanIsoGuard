#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// One STAR SJ.out.tab / combined_SJ.tab row.
//   col1 chrom, col2 intron_start(1-based), col3 intron_end(1-based, inclusive),
//   col4 strand{0=undef,1=+,2=-}, col5 motif{0=noncanonical,1..6 canonical},
//   col6 annotated{0,1}, col7 n_uniq, col8 n_multi, col9 max_overhang.
struct SjRecord {
  std::string chrom;
  Junction intron;  // normalized 0-based half-open
  Strand strand = Strand::kUnknown;
  int motif = 0;
  bool annotated = false;
  int n_uniq = 0;
  int n_multi = 0;
  int max_overhang = 0;

  bool canonical() const { return motif >= 1; }
};

// Short-read junction evidence with O(1) exact lookup by intron coordinates,
// matching the design's "exact donor-acceptor coordinate join".
class SjTable {
 public:
  void add(SjRecord rec);

  // Exact match on (chrom, intron.start, intron.end). Returns nullptr if absent.
  const SjRecord* find_exact(const std::string& chrom, const Junction& intron) const;

  std::size_t size() const { return records_.size(); }
  const std::vector<SjRecord>& records() const { return records_; }

 private:
  static std::string key(const std::string& chrom, const Junction& j);

  std::vector<SjRecord> records_;
  std::unordered_map<std::string, std::size_t> by_coord_;
};

// Parse a STAR SJ.tab (9 columns, no header). Throws std::runtime_error on
// open failure or malformed rows.
SjTable read_sj_tab(const std::string& path);

}  // namespace panisoguard
