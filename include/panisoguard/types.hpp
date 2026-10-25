#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace panisoguard {

// ---------------------------------------------------------------------------
// Coordinate convention (used by ALL readers and the catalog):
//   Every internal coordinate is 0-based, half-open.
//   A Junction is an INTRON [start, end): start = first intronic base,
//   end = one-past-the-last intronic base.
//
// Reader conversions that converge on this convention:
//   * STAR SJ.tab (1-based inclusive intron s..e)      -> {s-1, e}
//   * GTF exons (1-based inclusive s..e); intron between
//     exon_i and exon_{i+1}                            -> {exon_i.e, exon_{i+1}.s - 1}
//   * BED12 blocks (0-based, relative to chromStart);
//     intron between block_i and block_{i+1}           -> {block_i.end, block_{i+1}.start}
// so the SAME intron yields the SAME {start,end} regardless of source, which is
// what makes fingerprint matching and SJ.tab corroboration exact.
// ---------------------------------------------------------------------------

enum class Strand : char { kPlus = '+', kMinus = '-', kUnknown = '.' };

inline Strand parse_strand(char c) {
  switch (c) {
    case '+': return Strand::kPlus;
    case '-': return Strand::kMinus;
    default:  return Strand::kUnknown;
  }
}
inline char strand_char(Strand s) { return static_cast<char>(s); }

// An intron, 0-based half-open [start, end).
struct Junction {
  int64_t start = 0;
  int64_t end = 0;

  bool operator==(const Junction& o) const { return start == o.start && end == o.end; }
  bool operator!=(const Junction& o) const { return !(*this == o); }
  bool operator<(const Junction& o) const {
    return start != o.start ? start < o.start : end < o.end;
  }
};

// Ordered intron chain for one transcript / isoform.
struct IntronChain {
  std::string chrom;
  Strand strand = Strand::kUnknown;
  std::vector<Junction> introns;  // kept sorted ascending after sort_introns()

  bool monoexonic() const { return introns.empty(); }
  void sort_introns() { std::sort(introns.begin(), introns.end()); }
};

// "chrom|strand:start-end": the lookup key of every exact-coordinate index (reference
// catalog, SJ.tab, pangenome junctions, the BAM-feature cache). Strand-aware, so two
// opposite-strand junctions at the same coordinates never collide.
inline std::string junction_key(const std::string& chrom, Strand strand, const Junction& j) {
  return chrom + '|' + strand_char(strand) + ':' + std::to_string(j.start) + '-' + std::to_string(j.end);
}

// true if `sub` occurs as a contiguous run inside `full`: the test behind whole-chain read
// support (a read's introns contain the isoform's chain) and the decoy generator.
inline bool contiguous_subchain(const std::vector<Junction>& sub, const std::vector<Junction>& full) {
  if (sub.empty() || sub.size() > full.size()) return false;
  for (std::size_t i = 0; i + sub.size() <= full.size(); ++i) {
    bool ok = true;
    for (std::size_t k = 0; k < sub.size() && ok; ++k) ok = full[i + k] == sub[k];
    if (ok) return true;
  }
  return false;
}

// 64-bit intron-chain fingerprint (see fingerprint.hpp).
using Fingerprint = uint64_t;

}  // namespace panisoguard
