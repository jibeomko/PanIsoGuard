#pragma once

#include <string>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// One BED12 record decoded into an intron chain (e.g. a FLAIR isoform BED). `name` is
// the isoform id, the same string as SQANTI3's `isoform` column, so it is the join key
// between structure and QC.
struct Bed12Record {
  std::string name;
  IntronChain chain;
  int64_t start = 0;   // chromStart .. chromEnd, 0-based half-open
  int64_t end = 0;
};

// Parse a BED12 file. Throws std::runtime_error on open failure / malformed rows.
std::vector<Bed12Record> read_bed12(const std::string& path);

}  // namespace panisoguard
