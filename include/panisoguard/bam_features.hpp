#pragma once

#include <string>

#include "panisoguard/types.hpp"

namespace panisoguard {

// Read-level features at one novel junction, computed over reads whose alignment
// contains an N (ref-skip) op exactly matching the junction interval. High
// fractions of low-MAPQ / supplementary / indel-adjacent spanning
// reads indicate a mapping/chimera/variant artifact rather than a clean splice.
//
// Strand: a read whose aligner wrote the transcript strand (STAR/HISAT2 XS:A, minimap2 ts:A)
// counts only for a junction or chain on that strand; a read without the tag counts for
// either strand. A junction or chain of unknown strand accepts every read.
struct JunctionBamFeatures {
  bool evaluated = false;     // false if the contig is absent from the BAM header
  int n_spanning = 0;         // reads with an N op == this junction
  int n_low_mapq = 0;
  int n_supplementary = 0;    // SECONDARY or SUPPLEMENTARY flagged
  int n_indel_near = 0;       // I/D within junction_window_bp of either boundary

  double frac_low_mapq() const     { return n_spanning ? double(n_low_mapq) / n_spanning : 0.0; }
  double frac_supplementary() const{ return n_spanning ? double(n_supplementary) / n_spanning : 0.0; }
  double frac_indel_near() const   { return n_spanning ? double(n_indel_near) / n_spanning : 0.0; }
};

// Reads carrying an isoform's whole intron chain: the chain appears as a contiguous run
// of the read's N ops in one primary alignment (n_full); among the reads at the chain's
// first donor site, those whose own intron list is a contiguous part of the chain are
// consistent but too short to carry it (n_partial). Exact coordinates, like the junction
// features: a wobbling ONT junction does not count.
struct ChainBamSupport {
  bool evaluated = false;  // false if the contig is absent from the BAM header
  int n_full = 0;
  int n_partial = 0;
};

struct BamFeatureParams {
  int min_mapq = 20;
  int junction_window_bp = 10;
};

// RAII wrapper over an indexed BAM/CRAM (pImpl keeps htslib out of the header). Both
// queries throw std::runtime_error when htslib reports a read error (a truncated file is
// an error, not a region with no reads).
class BamReader {
 public:
  explicit BamReader(const std::string& path, const std::string& reference = "");
  ~BamReader();
  BamReader(const BamReader&) = delete;
  BamReader& operator=(const BamReader&) = delete;

  JunctionBamFeatures features_at_junction(const std::string& chrom, const Junction& intron, Strand strand,
                                           const BamFeatureParams& p) const;
  // Primary alignments only (secondary/supplementary skipped); monoexonic chains -> not evaluated.
  ChainBamSupport chain_support(const IntronChain& chain) const;

 private:
  struct Impl;
  Impl* impl_;
};

}  // namespace panisoguard
