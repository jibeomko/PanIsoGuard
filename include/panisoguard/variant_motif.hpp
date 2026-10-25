#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "panisoguard/types.hpp"

namespace panisoguard {

// htslib faidx wrapper (pImpl): fetch uppercase reference subsequences.
class FastaFetcher {
 public:
  explicit FastaFetcher(const std::string& path);  // requires/loads a .fai
  ~FastaFetcher();
  FastaFetcher(const FastaFetcher&) = delete;
  FastaFetcher& operator=(const FastaFetcher&) = delete;

  bool has_chrom(const std::string& chrom) const;
  // Uppercase sequence for 0-based half-open [start, end); "" if unavailable.
  std::string fetch(const std::string& chrom, int64_t start, int64_t end) const;

 private:
  struct Impl;
  Impl* impl_;
};

// The transcript-oriented donor/acceptor dinucleotides of an intron and whether
// they form a canonical splice motif (GT-AG, GC-AG, or AT-AC).
struct SpliceMotif {
  bool evaluable = false;
  bool canonical = false;
  std::string motif;  // e.g. "GT-AG"
};
SpliceMotif splice_motif(const FastaFetcher& fa, const std::string& chrom, const Junction& intron,
                         Strand strand);

// Per-junction variant verdict comparing the reference to >=1 personalized
// haplotype FASTA.
enum class VariantVerdict {
  kNotEvaluable,
  kNone,       // motif status unchanged between reference and haplotypes
  kCreated,    // non-canonical on reference, canonical on >=1 haplotype (reference bias)
  kDisrupted,  // canonical on reference, non-canonical on all haplotypes
};

// The motif of one intron on the reference and on every haplotype, with the verdict.
struct JunctionMotifs {
  SpliceMotif ref;
  std::vector<SpliceMotif> haps;  // in --reference-haplotype order
  VariantVerdict verdict = VariantVerdict::kNotEvaluable;
};

class HaplotypeProvider {
 public:
  HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                    std::vector<std::shared_ptr<FastaFetcher>> haplotypes);
  JunctionMotifs motifs(const std::string& chrom, const Junction& intron, Strand strand) const;
  std::size_t n_haplotypes() const { return haps_.size(); }

 private:
  std::shared_ptr<FastaFetcher> ref_;
  std::vector<std::shared_ptr<FastaFetcher>> haps_;
};

}  // namespace panisoguard
