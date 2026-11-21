#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
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
  int64_t length(const std::string& chrom) const;   // -1 when the contig is absent
  const std::string& path() const;

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
  int n_het = 0;                 // records that place alleles by haplotype index (heterozygous or partly missing)
                                 // and shaped a window base: two or more need one phase set
  long phase_set = kNoPhaseSet;
  bool mixed_sets = false;
  std::vector<std::string> counted;   // records already counted (position, alleles and genotype): one MNP over both windows is
                                      // one record, two records at one position are two
  std::vector<std::pair<std::string, std::string>> seen;   // (position + alleles, genotype) of the records read: a second record with
                                                           // the same alleles and another genotype contradicts the first
  void add_set(long ps) { if (phase_set == kNoPhaseSet) phase_set = ps; else if (phase_set != ps) mixed_sets = true; }
};

// The sample's records near one side of a junction, as each haplotype carries them. Every allele
// is an edit of the reference: a same-length substitution, an insertion, a deletion, another
// replacement, or "unknown" (a missing or undefined allele, a record without GT, a record whose
// REF is not the reference, two records that give one haplotype different alleles at one base).
// Insertions and deletions are moved as far towards the intron as the reference sequence allows,
// so the exon is read as the reference has it wherever that is possible.
struct VcfLocal {
  struct Edit {
    int64_t q0 = 0, q1 = 0;   // the reference interval [q0, q1) it replaces; q0 == q1 is an insertion before q0
    std::string ins;          // the bases that stand there instead
    int rec = -1;             // index into recs
    char kind = 's';          // s substitution, i insertion, d deletion, c other replacement, u unknown
    std::vector<int> also;    // further records folded into this edit (agreeing substitutions, or a clash)
    bool clash = false;       // two records gave this haplotype different alleles here
    bool star = false;        // unknown because of a `*` allele: the deletion that removes the base is another record
  };
  struct Rec {
    std::string id;           // position, alleles and what the genotype states
    bool index_dependent = false;   // heterozygous or partly missing: its alleles are placed by haplotype index
    std::string group;              // of an index-dependent record: its phase set, or the record itself when unphased
    bool phased = true;
    long set = -2;            // kDefaultPhaseSet
    int unknown_slots = 0;    // haplotypes it leaves unknown (missing or undefined allele; all of them without GT)
    bool conflict = false;    // the same site written again with another genotype
    bool ref_mismatch = false;
    int skipped = 0;          // symbolic alleles that cannot be applied
    bool indel = false;
    bool het_snv = false;
  };
  bool known_contig = false;
  int64_t lo = 0, hi = 0;
  std::string ref;                        // reference bases of [lo, hi)
  std::vector<std::vector<Edit>> edits;   // per haplotype: sorted, not overlapping
  std::vector<Rec> recs;
};

// One haplotype's two bases next to the last exonic base on one side of an intron.
struct HapSide {
  bool known = false;
  std::string two;             // the bases, intron-side reading
  bool boundary_insertion = false;   // an insertion sits exactly between exon and intron: `two` reads it as
  std::string two_exonic;            // intronic, `two_exonic` as part of the exon; only the reads can tell
  bool exonic_known = false;
  int insertion_len = 0;
  bool indel = false;          // an insertion or deletion shaped the bases
  std::vector<int> used;       // records the reading rests on
  std::vector<int> used_exonic;
};

// One heterozygous single-base site of the sample: the base each haplotype carries.
struct PhaseSite {
  int64_t pos = 0;
  std::string bases;           // one base per haplotype, in haplotype order
};

// One sample's haplotypes read from an indexed VCF/BCF, applied to the reference on demand.
// The sample's alleles are taken by position in the genotype (first allele = haplotype 1,
// second = haplotype 2), as `bcftools consensus -H 1` / `-H 2` build a FASTA; for a file of
// substitutions a run from the VCF equals a run from those FASTAs. Insertions and deletions
// are applied too (the FASTA route cannot follow them). What the file leaves open stays
// unknown (N), never reference: a missing allele, a genotype shorter than the ploidy, an allele
// index the record does not define, a record without GT, a REF that is not the reference, two
// records that give one haplotype different alleles at one base. An unphased heterozygous
// genotype is applied by allele order and counted, and the junction is marked so the
// per-haplotype claim can be weighed; the phase set of a phased allele is kept so that a claim
// spanning several records can require one set. The ploidy is the largest genotype in the file
// (read once at open). A contig missing from the VCF header is not evaluable; a contig present
// but without records is reference.
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

  // The records overlapping 0-based [lo, hi) as edits per haplotype. `ref_window` is the
  // reference sequence of [lo, hi); `reference` serves the bases beyond it that moving an indel
  // may need (may be null). `toward` is +1 when the intron lies to the right of the window's
  // anchor, -1 when to the left: indels are moved that way; 0 leaves them where the file has them.
  // `flipped` names groups (Rec::group) whose records are applied with the haplotypes in the
  // other order: the file does not say how two groups lie to each other, so each such order is a
  // reading of the same file.
  VcfLocal local(const std::string& chrom, int64_t lo, int64_t hi, const std::string& ref_window,
                 const FastaFetcher* reference, int toward, const std::vector<std::string>* flipped = nullptr) const;

  // The sample's haplotypes over 0-based half-open [start, end) base by base, given the
  // reference bases of that window: substitutions applied; a base that is unknown, deleted or
  // replaced by a sequence of another length reads N. Empty when the contig is not in the VCF
  // header. `phase` accumulates what the applied alleles said (see WindowPhase).
  std::vector<std::string> overlay(const std::string& chrom, int64_t start, int64_t end,
                                   const std::string& ref_window, WindowPhase& phase) const;

  // Single-base sites of [lo, hi) the file can vouch for: heterozygous ones whose alleles are
  // phased in `set` (phased genotypes without PS are in kDefaultPhaseSet; kNoPhaseSet: none), and
  // homozygous differences from the reference. `blocked` receives the span of every other
  // record, where the file's haplotypes cannot be compared base by base.
  void snv_sites(const std::string& chrom, int64_t lo, int64_t hi, const std::string& ref_window, long set,
                 std::vector<PhaseSite>& sites, std::vector<std::pair<int64_t, int64_t>>& blocked) const;

  // Counted over the records a junction's verdict rested on (count_used).
  void count_used(const VcfLocal& local, const std::vector<int>& used) const;
  void count_deletion_match() const;
  long n_unphased_het_applied() const;
  long n_unapplied_alleles() const;      // symbolic or `*` alleles that cannot be applied
  long n_indels_applied() const;
  long n_missing_alleles() const;        // missing or undefined alleles that left a motif base unknown
  long n_conflicting_records() const;    // sites written twice with another genotype, or two alleles on one base
  long n_ref_mismatches() const;         // records whose REF is not the reference
  long n_deletions_equal_to_junctions() const;

 private:
  struct Impl;
  Impl* impl_;
};

// The two bases after `anchor` (toward > 0) or before it (toward < 0) on haplotype `h`.
HapSide hap_side(const VcfLocal& local, std::size_t h, int64_t anchor, int toward);
// Whether haplotype `h` carries a deletion equal to the intron [start, end) (the same sequence
// results wherever the reference lets the deletion be written). `used` receives its record.
bool hap_deletion_equals(const VcfLocal& local, std::size_t h, int64_t start, int64_t end,
                         const FastaFetcher& reference, const std::string& chrom, int* used);

// Per-junction variant verdict comparing the reference to the sample's haplotypes.
enum class VariantVerdict {
  kNotEvaluable,
  kNone,       // motif status unchanged between reference and haplotypes
  kCreated,    // non-canonical on reference, canonical on >=1 haplotype (reference bias)
  kDisrupted,  // canonical on reference, non-canonical on all haplotypes
  kDeleted,    // a haplotype carries a deletion equal to the intron: the reads lack those bases
               // because the genome does, no splice is needed (reference bias)
};

// Whether the haplotype alleles behind a junction's motifs are placed with certainty: unknown
// for FASTA input or when no variant touched the motif; phased (a phased or homozygous
// genotype); unphased (a heterozygous genotype applied by allele order).
enum class MotifPhase { kUnknown, kPhased, kUnphased };

// What the junction reads say about an insertion that sits exactly at an exon boundary: they
// carry it (it is part of the exon), they do not (it is intronic), or they do not settle it.
// The oracle is asked per side (left = the 5' genomic end of the intron) and inserted length.
enum class BoundaryReads { kUnknown, kWithout, kWith };
using BoundaryOracle = std::function<BoundaryReads(bool left, int insertion_len)>;

// The motif of one intron on the reference and on every haplotype, with the verdict. For VCF
// input, `phase_set` is the phase set of the heterozygous alleles that shaped the motifs
// (kNoPhaseSet when none did, kDefaultPhaseSet for phased genotypes without PS) and
// `phase_mixed` says two sets met at this junction.
struct JunctionMotifs {
  SpliceMotif ref;
  std::vector<SpliceMotif> haps;  // in --reference-haplotype / --haplotype-vcf order
  std::vector<bool> deleted;      // per haplotype: carries a deletion equal to the intron
  VariantVerdict verdict = VariantVerdict::kNotEvaluable;
  MotifPhase phase = MotifPhase::kUnknown;
  long phase_set = kNoPhaseSet;
  bool phase_mixed = false;
  bool phase_unlinked = false;    // what the haplotypes carry here depends on how >= 2 groups of heterozygous
                                  // records lie to each other, which the file does not say (unphased, or
                                  // different sets): no haplotype's motif is known
  bool indel = false;             // an insertion or deletion shaped a haplotype's motif
  int boundary_insertion_left = 0;    // length of an insertion exactly at the left / right exon boundary that
  int boundary_insertion_right = 0;   // the reads did not place (0: none, or placed)
  std::vector<long> file_set;     // per haplotype source: the phase set its heterozygous records at this junction
                                  // share (kNoPhaseSet: none used; kUnusableSet: unphased or mixed)
  std::vector<int> source_of;     // per haplotype: its source (a FASTA, or one VCF file), an index into file_set
};
constexpr long kUnusableSet = -3;

class HaplotypeProvider {
 public:
  HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                    std::vector<std::shared_ptr<FastaFetcher>> haplotypes);
  HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                    std::vector<std::shared_ptr<VcfHaplotypes>> vcfs);
  // `reads` is asked only when a haplotype has an insertion exactly at an exon boundary; without
  // it such a haplotype's motif stays unknown.
  JunctionMotifs motifs(const std::string& chrom, const Junction& intron, Strand strand,
                        const BoundaryOracle& reads = nullptr) const;
  std::size_t n_haplotypes() const;
  // Single-base sites of [lo, hi) at which the haplotypes differ and every haplotype's base is
  // known, usable to tell which haplotype a read comes from for the junction `jm` describes.
  std::vector<PhaseSite> phase_sites(const std::string& chrom, int64_t lo, int64_t hi,
                                     const JunctionMotifs& jm) const;
  long n_length_mismatches() const { return n_length_mismatch_; }   // haplotype FASTA contigs of another length

 private:
  std::shared_ptr<FastaFetcher> ref_;
  std::vector<std::shared_ptr<FastaFetcher>> haps_;
  std::vector<std::shared_ptr<VcfHaplotypes>> vcfs_;
  mutable long n_length_mismatch_ = 0;
  mutable std::vector<std::string> warned_;
};

}  // namespace panisoguard
