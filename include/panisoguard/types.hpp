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

// 64-bit intron-chain fingerprint (see fingerprint.hpp).
using Fingerprint = uint64_t;

}  // namespace panisoguard
