#include "panisoguard/variant_motif.hpp"

#include <cctype>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>

#include <htslib/faidx.h>

namespace panisoguard {

struct FastaFetcher::Impl {
  faidx_t* fai = nullptr;
  std::string path;
};

FastaFetcher::FastaFetcher(const std::string& path) : impl_(new Impl) {
  impl_->path = path;
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

int64_t FastaFetcher::length(const std::string& chrom) const {
  if (!faidx_has_seq(impl_->fai, chrom.c_str())) return -1;
  return static_cast<int64_t>(faidx_seq_len64(impl_->fai, chrom.c_str()));
}

const std::string& FastaFetcher::path() const { return impl_->path; }

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

// ---- one haplotype near a boundary, base by base ------------------------------------------

namespace {

// One base of a haplotype: where it stands on the reference (-1 for an inserted base), the edit
// that put it there (-1: the reference), a deletion lying just before it, and for an inserted
// base the reference position it precedes.
struct Token {
  char base = 'N';
  int64_t ref = -1;
  int edit = -1;
  int gap = -1;
  int64_t ins_point = -1;
};

std::vector<Token> tokens_of(const VcfLocal& L, std::size_t h) {
  std::vector<Token> out;
  const std::vector<VcfLocal::Edit>& ed = L.edits[h];
  int64_t cur = L.lo;
  int gap = -1;
  auto ref_until = [&](int64_t to) {
    for (; cur < to && cur < L.hi; ++cur) {
      out.push_back(Token{static_cast<char>(std::toupper(static_cast<unsigned char>(L.ref[static_cast<std::size_t>(cur - L.lo)]))), cur, -1, gap, -1});
      gap = -1;
    }
  };
  for (std::size_t i = 0; i < ed.size(); ++i) {
    const VcfLocal::Edit& e = ed[i];
    const int ei = static_cast<int>(i);
    ref_until(e.q0);
    const int64_t a = std::max(e.q0, L.lo), b = std::min(e.q1, L.hi);
    switch (e.kind) {
      case 's':
        // a base the allele leaves as the reference has it is not this record's doing
        for (int64_t p = std::max(a, cur); p < b; ++p) {
          const char c = e.ins[static_cast<std::size_t>(p - e.q0)];
          const char r = static_cast<char>(std::toupper(static_cast<unsigned char>(L.ref[static_cast<std::size_t>(p - L.lo)])));
          out.push_back(Token{c, p, c == r ? -1 : ei, gap, -1});
          gap = -1;
        }
        break;
      case 'u':
        for (int64_t p = std::max(a, cur); p < b; ++p) { out.push_back(Token{'N', p, ei, gap, -1}); gap = -1; }
        break;
      case 'i':
        if (e.q0 >= L.lo && e.q0 <= L.hi)
          for (const char c : e.ins) { out.push_back(Token{c, -1, ei, gap, e.q0}); gap = -1; }
        break;
      case 'c':
        for (const char c : e.ins) { out.push_back(Token{c, -1, ei, gap, -1}); gap = -1; }
        break;
      default:   // 'd'
        gap = ei;
        break;
    }
    if (b > cur) cur = b;
  }
  ref_until(L.hi);
  return out;
}

bool is_base(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

void use_edit(const VcfLocal& L, std::size_t h, int ei, std::vector<int>& used, bool& indel) {
  if (ei < 0) return;
  const VcfLocal::Edit& e = L.edits[h][static_cast<std::size_t>(ei)];
  used.push_back(e.rec);
  used.insert(used.end(), e.also.begin(), e.also.end());
  if (e.kind == 'i' || e.kind == 'd' || e.kind == 'c') indel = true;
}

}  // namespace

HapSide hap_side(const VcfLocal& L, std::size_t h, int64_t anchor, int toward) {
  HapSide s;
  if (!L.known_contig || h >= L.edits.size() || L.ref.empty()) return s;
  const std::vector<Token> t = tokens_of(L, h);
  std::size_t at = t.size();
  for (std::size_t i = 0; i < t.size(); ++i) if (t[i].ref == anchor) { at = i; break; }
  if (at == t.size()) return s;   // the last exonic base is not on this haplotype (deleted or replaced)
  if (t[at].edit >= 0 && L.edits[h][static_cast<std::size_t>(t[at].edit)].kind == 'u') {
    bool unused = false;
    use_edit(L, h, t[at].edit, s.used, unused);   // what stands at the exon's last base is not known
    return s;
  }
  const int step = toward > 0 ? 1 : -1;
  // the two bases starting `skip` tokens beyond the anchor, with the records they rest on
  auto read = [&](std::size_t skip, std::string& two, std::vector<int>& used, bool& indel) {
    const int64_t first = static_cast<int64_t>(at) + step * static_cast<int64_t>(1 + skip);
    const int64_t second = first + step;
    if (first < 0 || second < 0 || first >= static_cast<int64_t>(t.size()) || second >= static_cast<int64_t>(t.size())) return false;
    // every token from the anchor to the second base: its edit, and the deletion just before it
    const int64_t a = std::min<int64_t>(static_cast<int64_t>(at), second), b = std::max<int64_t>(static_cast<int64_t>(at), second);
    for (int64_t i = a; i <= b; ++i) {
      const Token& k = t[static_cast<std::size_t>(i)];
      if (i != static_cast<int64_t>(at)) use_edit(L, h, k.edit, used, indel);
      if (i > a) use_edit(L, h, k.gap, used, indel);   // a deletion between this token and the one before it
    }
    const Token& x = t[static_cast<std::size_t>(first)];
    const Token& y = t[static_cast<std::size_t>(second)];
    if (!is_base(x.base) || !is_base(y.base)) return false;
    two = toward > 0 ? std::string{x.base, y.base} : std::string{y.base, x.base};
    return true;
  };
  // an insertion exactly between the exon's last base and the intron
  const int64_t nb = static_cast<int64_t>(at) + step;
  std::size_t ins_len = 0;
  if (nb >= 0 && nb < static_cast<int64_t>(t.size())) {
    const Token& n = t[static_cast<std::size_t>(nb)];
    if (n.ref < 0 && n.ins_point == (toward > 0 ? anchor + 1 : anchor) && n.gap < 0) {
      for (int64_t i = nb; i >= 0 && i < static_cast<int64_t>(t.size()) && t[static_cast<std::size_t>(i)].edit == n.edit && t[static_cast<std::size_t>(i)].ref < 0; i += step) ++ins_len;
    }
  }
  s.known = read(0, s.two, s.used, s.indel);
  if (ins_len > 0) {
    s.boundary_insertion = true;
    s.insertion_len = static_cast<int>(ins_len);
    bool indel = false;
    s.exonic_known = read(ins_len, s.two_exonic, s.used_exonic, indel);
  }
  return s;
}

bool hap_deletion_equals(const VcfLocal& L, std::size_t h, int64_t start, int64_t end,
                         const FastaFetcher& reference, const std::string& chrom, int* used) {
  if (!L.known_contig || h >= L.edits.size()) return false;
  const int64_t len = end - start;
  for (const VcfLocal::Edit& e : L.edits[h]) {
    if (e.kind != 'd' || e.q1 - e.q0 != len) continue;
    // the same sequence results when the bases between the two placements equal those `len` further on
    const int64_t a = std::min(e.q0, start), b = std::max(e.q0, start);
    if (b - a > 256) continue;
    if (b > a && reference.fetch(chrom, a, b) != reference.fetch(chrom, a + len, b + len)) continue;
    if (b > a && reference.fetch(chrom, a, b).empty()) continue;
    if (used != nullptr) *used = e.rec;
    return true;
  }
  return false;
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

namespace {
constexpr int64_t kSideWindow = 64;   // reference bases read on either side of a boundary
constexpr std::size_t kMaxGroups = 6; // linkage groups whose placements are all tried at one junction
}

JunctionMotifs HaplotypeProvider::motifs(const std::string& chrom, const Junction& intron, Strand strand,
                                         const BoundaryOracle& reads) const {
  JunctionMotifs m;
  const std::string left = ref_->fetch(chrom, intron.start, intron.start + 2);
  const std::string right = ref_->fetch(chrom, intron.end - 2, intron.end);
  m.ref = splice_motif_from(left, right, strand);
  m.haps.reserve(n_haplotypes());
  const int64_t ref_len = ref_->length(chrom);
  for (const auto& hap : haps_) {
    // a haplotype FASTA is read at the reference's coordinates: with another contig length they
    // are not the same coordinates, and nothing read there means anything
    const int64_t hap_len = hap->length(chrom);
    if (ref_len >= 0 && hap_len >= 0 && hap_len != ref_len) {
      const std::string key = hap->path() + "\t" + chrom;
      if (std::find(warned_.begin(), warned_.end(), key) == warned_.end()) {
        warned_.push_back(key);
        ++n_length_mismatch_;
        std::fprintf(stderr, "warning: %s: %s is %lld bp, the reference's is %lld bp; this haplotype is not evaluable there "
                             "(a haplotype FASTA must keep the reference's coordinates; give indels with --haplotype-vcf)\n",
                     hap->path().c_str(), chrom.c_str(), static_cast<long long>(hap_len), static_cast<long long>(ref_len));
      }
      m.haps.push_back(SpliceMotif{});
    } else {
      m.haps.push_back(splice_motif(*hap, chrom, intron, strand));
    }
    m.deleted.push_back(false);
    m.file_set.push_back(kNoPhaseSet);
    m.source_of.push_back(static_cast<int>(m.file_set.size()) - 1);
  }
  for (const auto& vcf : vcfs_) {
    const std::size_t nh = vcf->n_haplotypes();
    auto not_evaluable = [&](long set) {
      m.haps.insert(m.haps.end(), nh, SpliceMotif{});
      m.deleted.insert(m.deleted.end(), nh, false);
      m.file_set.push_back(set);
      m.source_of.insert(m.source_of.end(), nh, static_cast<int>(m.file_set.size()) - 1);
    };
    if (!vcf->has_chrom(chrom) || ref_len < 0 || intron.start < 1 || intron.end >= ref_len) { not_evaluable(kNoPhaseSet); continue; }
    // the bases next to the last exonic base on either side: intron.start - 1 and intron.end
    const int64_t l_lo = std::max<int64_t>(0, intron.start - 1 - kSideWindow), l_hi = std::min(ref_len, intron.start + 2 + kSideWindow);
    const int64_t r_lo = std::max<int64_t>(0, intron.end - 2 - kSideWindow), r_hi = std::min(ref_len, intron.end + 1 + kSideWindow);
    const std::string l_ref = ref_->fetch(chrom, l_lo, l_hi), r_ref = ref_->fetch(chrom, r_lo, r_hi);

    // What the file's haplotypes carry at this junction, with the groups in `flipped` turned over.
    struct Reading {
      VcfLocal L, R;
      std::vector<SpliceMotif> motifs;
      std::vector<bool> deleted;
      std::vector<int> used_l, used_r;
      bool indel = false;
      int unplaced_left = 0, unplaced_right = 0;
      std::string key;   // per haplotype: the four bases, DEL, or NA
    };
    auto read = [&](const std::vector<std::string>* flipped) {
      Reading g;
      g.L = vcf->local(chrom, l_lo, l_hi, l_ref, ref_.get(), +1, flipped);
      g.R = vcf->local(chrom, r_lo, r_hi, r_ref, ref_.get(), -1, flipped);
      g.motifs.assign(nh, SpliceMotif{});
      g.deleted.assign(nh, false);
      for (std::size_t h = 0; h < nh; ++h) {
        int del_rec = -1;
        if (hap_deletion_equals(g.L, h, intron.start, intron.end, *ref_, chrom, &del_rec)) {
          g.deleted[h] = true;
          g.used_l.push_back(del_rec);
          g.key += "DEL,";
          continue;
        }
        HapSide ls = hap_side(g.L, h, intron.start - 1, +1), rs = hap_side(g.R, h, intron.end, -1);
        // an insertion exactly at an exon boundary is intronic when the junction reads do not carry
        // it and part of the exon when they do; without reads the base next to the exon is not known
        auto place = [&](HapSide& s, bool at_left, int& unplaced) {
          if (!s.boundary_insertion) return;
          const BoundaryReads said = reads ? reads(at_left, s.insertion_len) : BoundaryReads::kUnknown;
          if (said == BoundaryReads::kWith) { s.two = s.two_exonic; s.known = s.exonic_known; s.used = s.used_exonic; }
          else if (said == BoundaryReads::kUnknown) { s.known = false; unplaced = std::max(unplaced, s.insertion_len); }
        };
        place(ls, true, g.unplaced_left);
        place(rs, false, g.unplaced_right);
        g.used_l.insert(g.used_l.end(), ls.used.begin(), ls.used.end());
        g.used_r.insert(g.used_r.end(), rs.used.begin(), rs.used.end());
        if (ls.known && rs.known) {
          g.motifs[h] = splice_motif_from(ls.two, rs.two, strand);
          if (ls.indel || rs.indel) g.indel = true;
          g.key += ls.two + rs.two + ",";
        } else {
          g.key += "NA,";
        }
      }
      return g;
    };
    Reading base = read(nullptr);

    // the records the file's own writing rests on, each once (a record may reach both sides)
    std::vector<const VcfLocal::Rec*> recs;
    std::vector<int> count_l, count_r;
    auto collect = [&](const VcfLocal& loc, std::vector<int>& used, std::vector<int>& counted) {
      std::sort(used.begin(), used.end());
      used.erase(std::unique(used.begin(), used.end()), used.end());
      for (const int r : used) {
        if (r < 0) continue;
        const VcfLocal::Rec& rec = loc.recs[static_cast<std::size_t>(r)];
        bool seen = false;
        for (const VcfLocal::Rec* o : recs) if (o->id == rec.id) seen = true;
        if (seen) continue;
        recs.push_back(&rec);
        counted.push_back(r);
      }
    };
    collect(base.L, base.used_l, count_l);
    collect(base.R, base.used_r, count_r);
    vcf->count_used(base.L, count_l);
    vcf->count_used(base.R, count_r);

    // Heterozygous (or partly missing) records are placed by haplotype index, and the file links
    // two of them only inside one phase set. A group matters here when turning it over changes
    // what the haplotypes carry under some placement of the other groups; with two groups that
    // matter, the bases above are a composition by allele order, not something the file says.
    std::vector<std::string> groups, matter;
    const VcfLocal::Rec* of_group = nullptr;
    for (const VcfLocal* loc : {&base.L, &base.R})
      for (const VcfLocal::Rec& rec : loc->recs)
        if (!rec.group.empty() && std::find(groups.begin(), groups.end(), rec.group) == groups.end()) groups.push_back(rec.group);
    bool open = false;
    if (nh == 2 && groups.size() <= kMaxGroups) {
      const unsigned n = 1u << groups.size();
      std::vector<std::string> keys(n);
      keys[0] = base.key;
      std::vector<std::string> flipped;
      for (unsigned o = 1; o < n; ++o) {
        flipped.clear();
        for (std::size_t g = 0; g < groups.size(); ++g) if ((o >> g) & 1u) flipped.push_back(groups[g]);
        keys[o] = read(&flipped).key;
      }
      for (std::size_t g = 0; g < groups.size(); ++g)
        for (unsigned o = 0; o < n; ++o)
          if (keys[o] != keys[o ^ (1u << g)]) { matter.push_back(groups[g]); break; }
    } else if (groups.size() >= 2) {
      open = true;   // more groups than are tried, or more than two haplotypes to permute
    } else {
      for (const VcfLocal::Rec* rec : recs) if (!rec->group.empty()) matter = {rec->group};
    }
    if (open || matter.size() >= 2) {
      not_evaluable(kUnusableSet);
      m.phase_unlinked = true;
      m.phase = MotifPhase::kUnphased;
      continue;
    }
    if (!matter.empty())
      for (const VcfLocal* loc : {&base.L, &base.R})
        for (const VcfLocal::Rec& rec : loc->recs) if (rec.group == matter[0]) of_group = &rec;
    const bool unphased = of_group != nullptr && !of_group->phased;
    const long set = of_group != nullptr && of_group->phased ? of_group->set : kNoPhaseSet;
    for (std::size_t h = 0; h < nh; ++h) {
      m.haps.push_back(base.motifs[h]);
      m.deleted.push_back(base.deleted[h]);
    }
    if (std::find(base.deleted.begin(), base.deleted.end(), true) != base.deleted.end()) vcf->count_deletion_match();
    if (base.indel) m.indel = true;
    m.boundary_insertion_left = std::max(m.boundary_insertion_left, base.unplaced_left);
    m.boundary_insertion_right = std::max(m.boundary_insertion_right, base.unplaced_right);
    m.file_set.push_back(of_group == nullptr ? kNoPhaseSet : unphased ? kUnusableSet : set);
    m.source_of.insert(m.source_of.end(), nh, static_cast<int>(m.file_set.size()) - 1);
    if (unphased) m.phase = MotifPhase::kUnphased;
    else if (!recs.empty() && m.phase == MotifPhase::kUnknown) m.phase = MotifPhase::kPhased;
    if (set != kNoPhaseSet) {
      if (m.phase_set == kNoPhaseSet) m.phase_set = set; else if (m.phase_set != set) m.phase_mixed = true;
    }
  }
  if (m.haps.empty()) return m;
  // a haplotype that lacks the intron's bases explains the junction whatever the motifs are
  if (std::find(m.deleted.begin(), m.deleted.end(), true) != m.deleted.end()) {
    m.verdict = VariantVerdict::kDeleted;
    return m;
  }
  if (!m.ref.evaluable) return m;
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

std::vector<PhaseSite> HaplotypeProvider::phase_sites(const std::string& chrom, int64_t lo, int64_t hi,
                                                      const JunctionMotifs& jm) const {
  std::vector<PhaseSite> out;
  const int64_t ref_len = ref_->length(chrom);
  if (ref_len < 0) return out;
  lo = std::max<int64_t>(0, lo);
  hi = std::min(hi, ref_len);
  if (hi <= lo) return out;
  const std::size_t nh = n_haplotypes();
  auto differs = [](const std::string& b) {
    for (const char c : b) if (!is_base(c)) return false;
    for (const char c : b) if (c != b[0]) return true;
    return false;
  };
  if (!haps_.empty()) {
    std::vector<std::string> seq;
    for (const auto& hap : haps_) {
      if (hap->length(chrom) != ref_len) return {};
      seq.push_back(hap->fetch(chrom, lo, hi));
      if (static_cast<int64_t>(seq.back().size()) != hi - lo) return {};
    }
    for (int64_t p = lo; p < hi; ++p) {
      PhaseSite s;
      s.pos = p;
      for (const auto& q : seq) s.bases += q[static_cast<std::size_t>(p - lo)];
      if (differs(s.bases)) out.push_back(std::move(s));
    }
    return out;
  }
  const std::string refw = ref_->fetch(chrom, lo, hi);
  if (static_cast<int64_t>(refw.size()) != hi - lo) return out;
  std::map<int64_t, std::string> at;
  std::vector<std::pair<int64_t, int64_t>> blocked;
  std::size_t first = 0;
  for (std::size_t f = 0; f < vcfs_.size(); ++f) {
    const std::size_t n = vcfs_[f]->n_haplotypes();
    const long set = f < jm.file_set.size() ? jm.file_set[f] : kNoPhaseSet;
    if (set == kUnusableSet || !vcfs_[f]->has_chrom(chrom)) return {};
    std::vector<PhaseSite> sites;
    vcfs_[f]->snv_sites(chrom, lo, hi, refw, set, sites, blocked);
    for (const PhaseSite& s : sites) {
      auto it = at.find(s.pos);
      if (it == at.end()) it = at.emplace(s.pos, std::string(nh, refw[static_cast<std::size_t>(s.pos - lo)])).first;
      for (std::size_t h = 0; h < n; ++h) it->second[first + h] = s.bases[h];
    }
    first += n;
  }
  for (const auto& kv : at) {
    bool skip = false;
    for (const auto& b : blocked) if (kv.first >= b.first && kv.first < b.second) { skip = true; break; }
    if (skip || !differs(kv.second)) continue;
    PhaseSite s;
    s.pos = kv.first;
    s.bases = kv.second;
    out.push_back(std::move(s));
  }
  return out;
}

}  // namespace panisoguard
