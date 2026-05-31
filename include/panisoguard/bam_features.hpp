#pragma once

#include <string>

#include "panisoguard/types.hpp"

namespace panisoguard {

// Read-level features at one novel junction, computed over reads whose alignment
// contains an N (ref-skip) op exactly matching the junction interval. High
// fractions of low-MAPQ / supplementary / soft-clipped / indel-adjacent spanning
// reads indicate a mapping/chimera/variant artifact rather than a clean splice.
struct JunctionBamFeatures {
  bool evaluated = false;     // false if the contig is absent from the BAM header
  int n_spanning = 0;         // reads with an N op == this junction
  int n_low_mapq = 0;
  int n_supplementary = 0;    // SECONDARY or SUPPLEMENTARY flagged
  int n_softclip = 0;         // terminal soft-clip >= softclip_min_bp
  int n_indel_near = 0;       // I/D within junction_window_bp of either boundary

  double frac_low_mapq() const     { return n_spanning ? double(n_low_mapq) / n_spanning : 0.0; }
  double frac_supplementary() const{ return n_spanning ? double(n_supplementary) / n_spanning : 0.0; }
  double frac_softclip() const     { return n_spanning ? double(n_softclip) / n_spanning : 0.0; }
  double frac_indel_near() const   { return n_spanning ? double(n_indel_near) / n_spanning : 0.0; }
};

struct BamFeatureParams {
  int min_mapq = 20;
  int softclip_min_bp = 20;
  int junction_window_bp = 10;
};

// RAII wrapper over an indexed BAM/CRAM (pImpl keeps htslib out of the header).
class BamReader {
 public:
  explicit BamReader(const std::string& path, const std::string& reference = "");
  ~BamReader();
  BamReader(const BamReader&) = delete;
  BamReader& operator=(const BamReader&) = delete;

  JunctionBamFeatures features_at_junction(const std::string& chrom, const Junction& intron,
                                           const BamFeatureParams& p) const;

 private:
  struct Impl;
  Impl* impl_;
};

}  // namespace panisoguard
