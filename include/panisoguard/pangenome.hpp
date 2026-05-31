#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// A set of splice junctions that are present on at least one haplotype PATH of a
// pangenome graph (e.g. HPRC), produced externally by vg/rpvg and supplied to
// PanIsoGuard as a coordinate file. This is the file/subprocess form of the
// pangenome tier: PanIsoGuard does not traverse the graph itself, it adjudicates
// using the graph's pre-extracted junction set (in-process gbwtgraph/GBZ remains
// behind -DWITH_PANGENOME_LIB).
//
// Used by the pangenome (reference-bias) rescue axis: a junction that is novel
// relative to the LINEAR reference but realizable on a pangenome haplotype path is
// consistent with reference bias rather than new splicing. This is independent
// (non-circular) evidence ONLY when the junction set is extracted from population
// assemblies; PanIsoGuard cannot verify the supplied file's provenance, so the
// rescue is gated by --pangenome-provenance and the same circularity firewall as
// the variant axis (a coordinate match shows the junction is possible on a path,
// it does not by itself prove the splice is used).
struct PangenomeJunctionRecord {
  std::string chrom;
  Junction intron;             // normalized 0-based half-open
  Strand strand = Strand::kUnknown;
  int n_haplotypes = 1;        // number of graph haplotypes carrying it (>= 1)
};

class PangenomeJunctions {
 public:
  void add(PangenomeJunctionRecord rec);

  // Exact match on (chrom, strand, intron.start, intron.end). Strand-aware, like
  // SjTable / the reference Catalog. Returns nullptr if absent.
  const PangenomeJunctionRecord* find_exact(const std::string& chrom, Strand strand,
                                            const Junction& intron) const;

  // True if the junction is present on >= min_haplotypes graph haplotypes.
  bool supports(const std::string& chrom, Strand strand, const Junction& intron,
                int min_haplotypes) const;

  std::size_t size() const { return records_.size(); }

 private:
  static std::string key(const std::string& chrom, Strand strand, const Junction& j);

  std::vector<PangenomeJunctionRecord> records_;
  std::unordered_map<std::string, std::size_t> by_coord_;
};

// Parse a pangenome junction file. Tab-separated, no header. Columns:
//   chrom  intron_start  intron_end  strand  [n_haplotypes]
// intron_start/intron_end are 1-based inclusive (first/last intronic base, the
// same convention as STAR SJ.tab and GTF introns); strand is "+"/"-" (1/2 also
// accepted). The optional 5th column is the supporting-haplotype count (default
// 1). Blank lines and lines beginning with '#' are ignored. Throws
// std::runtime_error on open failure or malformed rows.
PangenomeJunctions read_pangenome_junctions(const std::string& path);

}  // namespace panisoguard
