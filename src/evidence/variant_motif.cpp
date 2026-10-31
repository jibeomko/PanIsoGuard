#include "panisoguard/variant_motif.hpp"

#include <cctype>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#include <htslib/faidx.h>

namespace panisoguard {

struct FastaFetcher::Impl {
  faidx_t* fai = nullptr;
};

FastaFetcher::FastaFetcher(const std::string& path) : impl_(new Impl) {
  impl_->fai = fai_load(path.c_str());  // creates/loads <path>.fai
  if (impl_->fai == nullptr) {
    delete impl_;
    throw std::runtime_error("cannot load FASTA index for: " + path);
  }
}

FastaFetcher::~FastaFetcher() {
  if (impl_ != nullptr) {
    if (impl_->fai) fai_destroy(impl_->fai);
    delete impl_;
  }
}

bool FastaFetcher::has_chrom(const std::string& chrom) const {
  return faidx_has_seq(impl_->fai, chrom.c_str()) != 0;
}

std::string FastaFetcher::fetch(const std::string& chrom, int64_t start, int64_t end) const {
  if (end <= start || start < 0) return "";
  if (!faidx_has_seq(impl_->fai, chrom.c_str())) return "";
  hts_pos_t len = 0;
  // faidx uses 0-based inclusive coordinates [beg, end]; free() guaranteed via RAII.
  // Use the 64-bit fetch so coordinates past INT_MAX (chromosomes > ~2.1 Gbp) are not
  // silently truncated to a wrong (in-range) position.
  std::unique_ptr<char, void (*)(void*)> seq(
      faidx_fetch_seq64(impl_->fai, chrom.c_str(), static_cast<hts_pos_t>(start),
                        static_cast<hts_pos_t>(end - 1), &len),
      std::free);
  if (!seq || len <= 0) return "";
  std::string out(seq.get(), static_cast<std::size_t>(len));
  for (char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return out;
}

namespace {

char comp(char b) {
  switch (b) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default:  return 'N';
  }
}
std::string revcomp(const std::string& s) {
  std::string r(s.rbegin(), s.rend());
  for (char& c : r) c = comp(c);
  return r;
}

// Canonical pairs of (donor, acceptor) dinucleotides in transcript orientation.
bool is_canonical_pair(const std::string& donor, const std::string& acceptor) {
  return (donor == "GT" && acceptor == "AG") || (donor == "GC" && acceptor == "AG") ||
         (donor == "AT" && acceptor == "AC");
}

// A dinucleotide with an N or any other ambiguity code cannot be judged canonical or not.
bool is_acgt(const std::string& s) {
  for (char c : s)
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
  return true;
}

}  // namespace

SpliceMotif splice_motif(const FastaFetcher& fa, const std::string& chrom, const Junction& intron,
                         Strand strand) {
  SpliceMotif m;
  if (strand == Strand::kUnknown) return m;  // donor and acceptor cannot be told apart: not evaluable
  const std::string left = fa.fetch(chrom, intron.start, intron.start + 2);   // genomic 5' of intron
  const std::string right = fa.fetch(chrom, intron.end - 2, intron.end);      // genomic 3' of intron
  if (left.size() != 2 || right.size() != 2 || !is_acgt(left) || !is_acgt(right)) return m;  // missing or ambiguous
  m.evaluable = true;

  std::string donor, acceptor;
  if (strand == Strand::kMinus) {
    // On the minus strand the transcript donor is the genomic right end (revcomp).
    donor = revcomp(right);
    acceptor = revcomp(left);
  } else {
    donor = left;
    acceptor = right;
  }
  m.motif = donor + "-" + acceptor;
  m.canonical = is_canonical_pair(donor, acceptor);
  return m;
}

HaplotypeProvider::HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                                     std::vector<std::shared_ptr<FastaFetcher>> haplotypes)
    : ref_(std::move(reference)), haps_(std::move(haplotypes)) {}

JunctionMotifs HaplotypeProvider::motifs(const std::string& chrom, const Junction& intron,
                                         Strand strand) const {
  JunctionMotifs m;
  m.ref = splice_motif(*ref_, chrom, intron, strand);
  m.haps.reserve(haps_.size());
  for (const auto& hap : haps_) m.haps.push_back(splice_motif(*hap, chrom, intron, strand));
  if (!m.ref.evaluable || m.haps.empty()) return m;
  bool all_hap_evaluable = true, any_hap_canonical = false;
  for (const auto& h : m.haps) {
    if (!h.evaluable) all_hap_evaluable = false;
    else if (h.canonical) any_hap_canonical = true;
  }
  // One canonical haplotype is positive evidence on its own; "none" and "disrupted" are
  // statements about every haplotype, so an unreadable haplotype leaves them undecided.
  if (!m.ref.canonical) {
    m.verdict = any_hap_canonical ? VariantVerdict::kCreated
              : all_hap_evaluable ? VariantVerdict::kNone : VariantVerdict::kNotEvaluable;
  } else {
    m.verdict = any_hap_canonical ? VariantVerdict::kNone
              : all_hap_evaluable ? VariantVerdict::kDisrupted : VariantVerdict::kNotEvaluable;
  }
  return m;
}

}  // namespace panisoguard
