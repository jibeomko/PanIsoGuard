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
  return splice_motif_from(fa.fetch(chrom, intron.start, intron.start + 2),   // genomic 5' of intron
                           fa.fetch(chrom, intron.end - 2, intron.end), strand);  // genomic 3' of intron
}

SpliceMotif splice_motif_from(const std::string& left, const std::string& right, Strand strand) {
  SpliceMotif m;
  if (strand == Strand::kUnknown) return m;  // donor and acceptor cannot be told apart: not evaluable
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

HaplotypeProvider::HaplotypeProvider(std::shared_ptr<FastaFetcher> reference,
                                     std::vector<std::shared_ptr<VcfHaplotypes>> vcfs)
    : ref_(std::move(reference)), vcfs_(std::move(vcfs)) {}

std::size_t HaplotypeProvider::n_haplotypes() const {
  std::size_t n = haps_.size();
  for (const auto& v : vcfs_) n += v->n_haplotypes();
  return n;
}

JunctionMotifs HaplotypeProvider::motifs(const std::string& chrom, const Junction& intron,
                                         Strand strand) const {
  JunctionMotifs m;
  const std::string left = ref_->fetch(chrom, intron.start, intron.start + 2);
  const std::string right = ref_->fetch(chrom, intron.end - 2, intron.end);
  m.ref = splice_motif_from(left, right, strand);
  m.haps.reserve(n_haplotypes());
  for (const auto& hap : haps_) m.haps.push_back(splice_motif(*hap, chrom, intron, strand));
  for (const auto& vcf : vcfs_) {
    WindowPhase wp;   // one record for both windows: the phase sets of a junction are read together
    const std::vector<std::string> lefts = left.size() == 2 ? vcf->overlay(chrom, intron.start, intron.start + 2, left, wp) : std::vector<std::string>{};
    const std::vector<std::string> rights = right.size() == 2 ? vcf->overlay(chrom, intron.end - 2, intron.end, right, wp) : std::vector<std::string>{};
    if (lefts.size() != vcf->n_haplotypes() || rights.size() != vcf->n_haplotypes()) {
      m.haps.insert(m.haps.end(), vcf->n_haplotypes(), SpliceMotif{});   // contig unknown to the VCF: not evaluable
      continue;
    }
    if (wp.n_het >= 2 && (wp.unphased || wp.mixed_sets)) {
      // two heterozygous records shape this motif and the file does not say they lie on the same
      // haplotype: the per-haplotype motifs are compositions by allele order, not observations
      m.haps.insert(m.haps.end(), vcf->n_haplotypes(), SpliceMotif{});
      m.phase_unlinked = true;
      m.phase = MotifPhase::kUnphased;
      continue;
    }
    for (std::size_t h = 0; h < lefts.size(); ++h) m.haps.push_back(splice_motif_from(lefts[h], rights[h], strand));
    const bool changed = [&] {
      for (std::size_t h = 0; h < lefts.size(); ++h) if (lefts[h] != left || rights[h] != right) return true;
      return false;
    }();
    if (wp.unphased) m.phase = MotifPhase::kUnphased;
    else if (changed && m.phase == MotifPhase::kUnknown) m.phase = MotifPhase::kPhased;
    if (wp.phase_set != kNoPhaseSet) {
      if (m.phase_set == kNoPhaseSet) m.phase_set = wp.phase_set; else if (m.phase_set != wp.phase_set) m.phase_mixed = true;
    }
    if (wp.mixed_sets) m.phase_mixed = true;
  }
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
