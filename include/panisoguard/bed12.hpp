#pragma once

#include <string>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// One BED12 record decoded into an intron chain (used for FLAIR3 combined.bed
// and read-level *_inconsistent / *_cannot_verify BEDs).
struct Bed12Record {
  std::string name;
  std::string id_prefix;  // token before the first '_' in the FLAIR id (see below)
  int score = 0;          // BED col5; FLAIR read-level BEDs carry MAPQ here
  IntronChain chain;
};

// FLAIR combined.bed isoform IDs look like "1-1_ENST..._ENSG..." -- the leading
// token is a FLAIR running id ("{n}-1"), NOT a biological sample (verified on the
// real cohort: combined.bed renumbers merged isoforms, and the BAM carries no
// @RG SM). The same string is the SQANTI3 `isoform` id, so it is the join key
// between structure and QC. Biological-sample provenance must come from
// per-sample inputs or a count matrix (deferred milestone), not this prefix.
inline std::string flair_id_prefix(const std::string& name) {
  const std::size_t us = name.find('_');
  return us == std::string::npos ? name : name.substr(0, us);
}

// Parse a BED12 file. Throws std::runtime_error on open failure / malformed rows.
std::vector<Bed12Record> read_bed12(const std::string& path);

}  // namespace panisoguard
