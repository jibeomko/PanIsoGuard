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
// The same from the two genomic dinucleotides already in hand (5' end of the intron, 3' end).
SpliceMotif splice_motif_from(const std::string& left, const std::string& right, Strand strand);

// What the alleles applied inside one window said about phase and completeness: an unphased
// heterozygous allele was used; an allele was missing (`.`, or a genotype shorter than the
// sample's ploidy, in which case the base is unknown, not reference); the phase set (FORMAT PS)
// of the phased heterozygous alleles, kDefaultPhaseSet for phased genotypes without PS,
// kNoPhaseSet when no heterozygous allele was applied, and `mixed_sets` when two sets met.
constexpr long kNoPhaseSet = -1;
constexpr long kDefaultPhaseSet = -2;
struct WindowPhase {
  bool unphased = false;
  bool missing = false;
  long phase_set = kNoPhaseSet;
  bool mixed_sets = false;
  void add_set(long ps) { if (phase_set == kNoPhaseSet) phase_set = ps; else if (phase_set != ps) mixed_sets = true; }
};

// One sample's haplotypes read from an indexed VCF/BCF, applied to the reference on demand.
// The sample's alleles are taken by position in the genotype (first allele = haplotype 1,
// second = haplotype 2), exactly as `bcftools consensus -H 1` / `-H 2` build a FASTA, so a
// run from the VCF equals a run from those FASTAs. Only same-length substitutions (SNVs,
// MNPs) are applied; indels and symbolic alleles are skipped and counted. A missing allele
// (`.`) or a genotype shorter than the ploidy leaves the base unknown (N) on that haplotype,
// so the motif there is not evaluable rather than reference. An unphased heterozygous
// genotype is applied by allele order and counted, and the junction is marked so the
// per-haplotype claim can be weighed; the phase set of a phased allele is kept so that a
// claim spanning several junctions can require one set. The ploidy is the largest genotype in
// the file (read once at open). A contig missing from the VCF header is not evaluable; a
// contig present but without records is reference.
class VcfHaplotypes {
 public:
  // `sample` may be empty when the file holds one sample. Throws when the file cannot be
  // opened, has no index, has no sample, or names several samples and none is chosen.
  explicit VcfHaplotypes(const std::string& path, const std::string& sample = "");
  ~VcfHaplotypes();
  VcfHaplotypes(const VcfHaplotypes&) = delete;
  VcfHaplotypes& operator=(const VcfHaplotypes&) = delete;

  std::size_t n_haplotypes() const;   // the sample's ploidy
  bool has_chrom(const std::string& chrom) const;
  const std::string& path() const;
  const std::string& sample() const;

  // The sample's haplotypes over 0-based half-open [start, end), given the reference bases of
  // that window. Empty when the contig is not in the VCF header. `phase` accumulates what the
  // applied alleles said (see WindowPhase); pass one object for both windows of a junction.
  std::vector<std::string> overlay(const std::string& chrom, int64_t start, int64_t end,
                                   const std::string& ref_window, WindowPhase& phase) const;

  long n_unphased_het_applied() const;   // over every window queried so far
  long n_non_snv_skipped() const;
  long n_missing_alleles() const;        // missing or absent alleles that left a motif base unknown

 private:
  struct Impl;
  Impl* impl_;
};

// Per-junction variant verdict comparing the reference to >=1 personalized
// haplotype FASTA.
enum class VariantVerdict {
  kNotEvaluable,
  kNone,       // motif status unchanged between reference and haplotypes
  kCreated,    // non-canonical on reference, canonical on >=1 haplotype (reference bias)
  kDisrupted,  // canonical on reference, non-canonical on all haplotypes
};

// Whether the haplotype alleles behind a junction's motifs are placed with certainty: unknown
// for FASTA input or when no variant touched the motif; phased (a phased or homozygous
// genotype); unphased (a heterozygous genotype applied by allele order).
enum class MotifPhase { kUnknown, kPhased, kUnphased };

// The motif of one intron on the reference and on every haplotype, with the verdict. For VCF
// input, `phase_set` is the phase set of the heterozygous alleles that shaped the motifs
// (kNoPhaseSet when none did, kDefaultPhaseSet for phased genotypes without PS) and
// `phase_mixed` says two sets met at this junction.
struct JunctionMotifs {
  SpliceMotif ref;
  std::vector<SpliceMotif> haps;  // in --reference-haplotype / --haplotype-vcf order
  VariantVerdict verdict = VariantVerdict::kNotEvaluable;
  MotifPhase phase = MotifPhase::kUnknown;
  long phase_set = kNoPhaseSet;
  bool phase_mixed = false;
};

class HaplotypeProvider {
 public:
  HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                    std::vector<std::shared_ptr<FastaFetcher>> haplotypes);
  HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                    std::vector<std::shared_ptr<VcfHaplotypes>> vcfs);
  JunctionMotifs motifs(const std::string& chrom, const Junction& intron, Strand strand) const;
  std::size_t n_haplotypes() const;

 private:
  std::shared_ptr<FastaFetcher> ref_;
  std::vector<std::shared_ptr<FastaFetcher>> haps_;
  std::vector<std::shared_ptr<VcfHaplotypes>> vcfs_;
};

}  // namespace panisoguard
