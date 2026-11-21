#include "panisoguard/variant_motif.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <stdexcept>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

namespace panisoguard {

namespace {

inline char up(char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }
inline bool acgt(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }
std::string upper(const char* s) {
  std::string out(s);
  for (char& c : out) c = up(c);
  return out;
}

constexpr int kMaxShift = 256;   // how far an insertion or deletion is moved through a repeat

// One record of the file as read: its alleles, and the allele index of every haplotype
// (-1: unknown, i.e. missing, absent from a shorter genotype, or an index the record lacks).
struct RawRecord {
  int64_t pos = 0, end = 0;            // reference span [pos, end)
  std::vector<std::string> alleles;    // upper case
  std::vector<int> allele;             // per haplotype
  bool has_gt = false, phased = true, het = false;
  bool sized = false;                  // a symbolic allele whose end the record gives
  int n_known = 0, n_unknown = 0;
  long set = kDefaultPhaseSet;
  std::string key, statement;
};

// The records of one window as read from the file.
struct RawWindow {
  std::string chrom;
  int64_t lo = -1, hi = -1;
  std::vector<RawRecord> raws;
};

}  // namespace

struct VcfHaplotypes::Impl {
  bcf_srs_t* sr = nullptr;      // owns the file, its index and header
  bcf_hdr_t* hdr = nullptr;     // borrowed from sr
  std::string path, sample;
  int n_samples = 0;
  int sample_idx = 0;
  int ploidy = 2;
  mutable int* gt = nullptr;    // bcf_get_genotypes buffer, grown by htslib
  mutable int gt_cap = 0;
  mutable long n_unphased = 0;
  mutable long n_non_snv = 0;
  mutable long n_indel = 0;
  mutable long n_missing = 0;
  mutable long n_conflicting = 0;
  mutable long n_ref_mismatch = 0;
  mutable long n_deletion_match = 0;
  mutable int32_t* ps = nullptr;   // bcf_get_format_int32 buffer for PS and END
  mutable int ps_cap = 0;
  mutable char** ps_str = nullptr;   // bcf_get_format_string buffer for a PS declared as a String
  mutable int ps_str_cap = 0;
  mutable std::unordered_map<std::string, long> ps_names;   // the sets such a file names, in the order met
  // the two windows read last (a junction's left and right side are built several times)
  mutable RawWindow cache[2];
  mutable int cache_next = 0;
  const std::vector<RawRecord>& window(const std::string& chrom, int rid, int64_t lo, int64_t hi, std::size_t nh) const;
  bool read_failed() const { return sr->errnum == vcf_parse_error || sr->errnum == bcf_read_error; }
  // The set a String PS names: a number is the set of that number, as in an Integer field; any
  // other value (some call sets write names, such as PATMAT) is a set of its own.
  long named_set(const std::string& v) const {
    if (v.size() <= 9 && std::all_of(v.begin(), v.end(), [](unsigned char c) { return std::isdigit(c) != 0; })) return std::stol(v);
    auto it = ps_names.find(v);
    if (it == ps_names.end()) it = ps_names.emplace(v, (1L << 40) + static_cast<long>(ps_names.size())).first;
    return it->second;
  }
};

VcfHaplotypes::VcfHaplotypes(const std::string& path, const std::string& sample) : impl_(new Impl) {
  impl_->path = path;
  impl_->sample = sample;
  // Format and index are checked up front: a synced reader whose add_reader failed cannot be
  // destroyed safely in every htslib release, so the failing cases never reach it.
  {
    htsFile* fp = hts_open(path.c_str(), "r");
    if (fp == nullptr) { delete impl_; throw std::runtime_error("cannot open --haplotype-vcf " + path); }
    const htsFormat fmt = *hts_get_format(fp);
    hts_close(fp);
    if (fmt.format != vcf && fmt.format != bcf) {
      delete impl_;
      throw std::runtime_error("--haplotype-vcf " + path + " is not a VCF or BCF file");
    }
    bool indexed = false;
    if (fmt.format == bcf) {
      hts_idx_t* idx = hts_idx_load3(path.c_str(), nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
      indexed = idx != nullptr;
      if (idx) hts_idx_destroy(idx);
    } else {
      tbx_t* tbx = tbx_index_load3(path.c_str(), nullptr, HTS_IDX_SILENT_FAIL);
      indexed = tbx != nullptr;
      if (tbx) tbx_destroy(tbx);
    }
    if (!indexed) {
      delete impl_;
      throw std::runtime_error("--haplotype-vcf " + path + " has no index; bgzip it and run tabix -p vcf (or bcftools index)");
    }
  }
  impl_->sr = bcf_sr_init();
  bcf_sr_set_opt(impl_->sr, BCF_SR_REQUIRE_IDX);
  if (!bcf_sr_add_reader(impl_->sr, path.c_str())) {
    const std::string why = bcf_sr_strerror(impl_->sr->errnum);
    impl_->sr = nullptr;   // not destroyed: see above
    delete impl_;
    throw std::runtime_error("cannot read --haplotype-vcf " + path + " (" + why + ")");
  }
  // A BGZF file without its end-of-file block was cut short, and the records after the cut
  // would read as reference.
  if (impl_->sr->errnum == no_eof) {
    bcf_sr_destroy(impl_->sr); delete impl_;
    throw std::runtime_error("--haplotype-vcf " + path + " has no BGZF end-of-file block: the file is truncated "
                             "(write it again with bgzip or bcftools view)");
  }
  impl_->hdr = bcf_sr_get_header(impl_->sr, 0);
  impl_->n_samples = bcf_hdr_nsamples(impl_->hdr);
  if (impl_->n_samples == 0) {
    bcf_sr_destroy(impl_->sr); delete impl_;
    throw std::runtime_error("--haplotype-vcf " + path + " has no sample column; genotypes are needed");
  }
  if (!sample.empty()) {
    impl_->sample_idx = bcf_hdr_id2int(impl_->hdr, BCF_DT_SAMPLE, sample.c_str());
    if (impl_->sample_idx < 0) {
      bcf_sr_destroy(impl_->sr); delete impl_;
      throw std::runtime_error("--vcf-sample " + sample + " is not in " + path);
    }
  } else if (impl_->n_samples > 1) {
    const int n = impl_->n_samples;   // read before the Impl is freed
    bcf_sr_destroy(impl_->sr); delete impl_;
    throw std::runtime_error("--haplotype-vcf " + path + " holds " + std::to_string(n) + " samples; name one with --vcf-sample");
  }
  // The sample's ploidy: the largest genotype in the file, read once here, so a haploid
  // record anywhere (a male X, a mitochondrial contig) neither hides a haplotype nor invents
  // one. A file without records keeps the default 2.
  int ploidy = 0;
  while (bcf_sr_next_line(impl_->sr)) {
    bcf1_t* rec = bcf_sr_get_line(impl_->sr, 0);
    const int ngt = bcf_get_genotypes(impl_->hdr, rec, &impl_->gt, &impl_->gt_cap);
    if (ngt <= 0) continue;
    const int per = ngt / impl_->n_samples;
    int p = 0;
    for (int h = 0; h < per; ++h) {
      if (impl_->gt[impl_->sample_idx * per + h] == bcf_int32_vector_end) break;
      ++p;
    }
    if (p > ploidy) ploidy = p;
  }
  if (ploidy > 0) impl_->ploidy = ploidy;
  // A record htslib cannot read ends this scan and would end every later read at the same row:
  // the records after it would read as reference. FORMAT/PS is read as the Integer the
  // specification defines, or as the String some call sets declare (and htslib assumes when the
  // header does not declare it); as another type its values could not be read, and every phased
  // genotype would fall into one set.
  std::string bad;
  if (impl_->read_failed())
    bad = std::string(" holds a record that cannot be read (") + bcf_sr_strerror(impl_->sr->errnum) + "; htslib names it above)";
  const int ps_id = bcf_hdr_id2int(impl_->hdr, BCF_DT_ID, "PS");
  if (bad.empty() && ps_id >= 0 && bcf_hdr_idinfo_exists(impl_->hdr, BCF_HL_FMT, ps_id) &&
      bcf_hdr_id2type(impl_->hdr, BCF_HL_FMT, ps_id) != BCF_HT_INT && bcf_hdr_id2type(impl_->hdr, BCF_HL_FMT, ps_id) != BCF_HT_STR)
    bad = ": FORMAT/PS is declared as neither an Integer nor a String, so its phase sets cannot be read";
  if (!bad.empty()) {
    if (impl_->gt) free(impl_->gt);
    bcf_sr_destroy(impl_->sr); delete impl_;
    throw std::runtime_error("--haplotype-vcf " + path + bad);
  }
}

VcfHaplotypes::~VcfHaplotypes() {
  if (impl_ != nullptr) {
    if (impl_->gt) free(impl_->gt);
    if (impl_->ps) free(impl_->ps);
    if (impl_->ps_str) { free(impl_->ps_str[0]); free(impl_->ps_str); }
    if (impl_->sr) bcf_sr_destroy(impl_->sr);
    delete impl_;
  }
}

std::size_t VcfHaplotypes::n_haplotypes() const { return static_cast<std::size_t>(impl_->ploidy); }
bool VcfHaplotypes::has_chrom(const std::string& chrom) const { return bcf_hdr_name2id(impl_->hdr, chrom.c_str()) >= 0; }
const std::string& VcfHaplotypes::path() const { return impl_->path; }
const std::string& VcfHaplotypes::sample() const { return impl_->sample; }
long VcfHaplotypes::n_unphased_het_applied() const { return impl_->n_unphased; }
long VcfHaplotypes::n_unapplied_alleles() const { return impl_->n_non_snv; }
long VcfHaplotypes::n_indels_applied() const { return impl_->n_indel; }
long VcfHaplotypes::n_missing_alleles() const { return impl_->n_missing; }
long VcfHaplotypes::n_conflicting_records() const { return impl_->n_conflicting; }
long VcfHaplotypes::n_ref_mismatches() const { return impl_->n_ref_mismatch; }
long VcfHaplotypes::n_deletions_equal_to_junctions() const { return impl_->n_deletion_match; }
void VcfHaplotypes::count_deletion_match() const { ++impl_->n_deletion_match; }

void VcfHaplotypes::count_used(const VcfLocal& local, const std::vector<int>& used) const {
  for (const int r : used) {
    if (r < 0 || r >= static_cast<int>(local.recs.size())) continue;
    const VcfLocal::Rec& rec = local.recs[static_cast<std::size_t>(r)];
    impl_->n_missing += rec.unknown_slots;
    impl_->n_non_snv += rec.skipped;
    if (rec.conflict) ++impl_->n_conflicting;
    if (rec.ref_mismatch) ++impl_->n_ref_mismatch;
    if (rec.indel) ++impl_->n_indel;
    if (rec.index_dependent && !rec.phased) ++impl_->n_unphased;
  }
}

const std::vector<RawRecord>& VcfHaplotypes::Impl::window(const std::string& chrom, int rid, int64_t lo, int64_t hi,
                                                          std::size_t nh) const {
  for (const RawWindow& w : cache) if (w.lo == lo && w.hi == hi && w.chrom == chrom) return w.raws;
  RawWindow& w = cache[cache_next];
  cache_next ^= 1;
  w.chrom = chrom; w.lo = lo; w.hi = hi;
  std::vector<RawRecord>& raws = w.raws;
  raws.clear();
  // A seek that fails names a contig the index has no records for: the window is reference.
  if (bcf_sr_seek(sr, chrom.c_str(), static_cast<hts_pos_t>(lo > 0 ? lo - 1 : 0)) < 0) return raws;
  while (bcf_sr_next_line(sr)) {
    bcf1_t* rec = bcf_sr_get_line(sr, 0);
    if (rec->rid != rid || rec->pos >= hi) break;
    bcf_unpack(rec, BCF_UN_STR);
    RawRecord r;
    r.pos = rec->pos;
    r.end = rec->pos + static_cast<int64_t>(std::strlen(rec->d.allele[0]));
    for (int a = 0; a < rec->n_allele; ++a) r.alleles.push_back(upper(rec->d.allele[a]));
    // a symbolic allele says where it ends in INFO/END (a deletion also by its length); one that
    // does not leaves the bases after it open
    bool symbolic = false, sym_del = false;
    for (int a = 1; a < rec->n_allele; ++a) {
      const std::string& alt = r.alleles[static_cast<std::size_t>(a)];
      if (!alt.empty() && alt[0] == '<') symbolic = true;
      if (alt.rfind("<DEL", 0) == 0) sym_del = true;
    }
    if (symbolic) {
      bool ended = false;
      int n = bcf_get_info_int32(hdr, rec, "END", &ps, &ps_cap);
      if (n > 0 && ps[0] != bcf_int32_missing && ps[0] > rec->pos + 1) { r.end = ps[0]; ended = true; }
      else if (sym_del) {
        n = bcf_get_info_int32(hdr, rec, "SVLEN", &ps, &ps_cap);
        if (n > 0 && ps[0] != bcf_int32_missing && ps[0] != 0) {
          r.end = rec->pos + 1 + (ps[0] < 0 ? -static_cast<int64_t>(ps[0]) : static_cast<int64_t>(ps[0]));
          ended = true;
        }
      }
      r.sized = ended;
      if (!ended) r.end = rec->pos + 1 + kMaxShift;
    }
    if (r.end <= lo) continue;   // ends before the window
    r.key = std::to_string(r.pos);
    for (std::size_t a = 0; a < r.alleles.size(); ++a) { r.key += a ? ',' : ':'; r.key += r.alleles[a]; }
    const int ngt = bcf_get_genotypes(hdr, rec, &gt, &gt_cap);
    r.allele.assign(nh, -1);
    if (ngt > 0) {
      r.has_gt = true;
      const int per = ngt / n_samples;
      const int* g = gt + sample_idx * per;
      bool ended = false, slash = false;
      int first = -1;
      std::vector<std::string> said;
      for (int h = 0; h < static_cast<int>(nh); ++h) {
        if (h >= per || g[h] == bcf_int32_vector_end) ended = true;
        if (!ended && h > 0 && !bcf_gt_is_phased(g[h])) { r.phased = false; slash = true; }   // the separator before a missing allele counts too
        if (ended) { ++r.n_unknown; continue; }
        if (bcf_gt_is_missing(g[h])) { ++r.n_unknown; said.push_back("."); continue; }
        const int a = bcf_gt_allele(g[h]);
        said.push_back(std::to_string(a));
        if (a >= rec->n_allele) { ++r.n_unknown; continue; }   // an allele the record does not define
        r.allele[static_cast<std::size_t>(h)] = a;
        ++r.n_known;
        if (first < 0) first = a; else if (a != first) r.het = true;
      }
      const int nps = bcf_get_format_int32(hdr, rec, "PS", &ps, &ps_cap);
      if (nps > 0) {
        const int32_t v = ps[sample_idx * (nps / n_samples)];
        if (v != bcf_int32_missing && v != bcf_int32_vector_end) r.set = v;
      } else if (nps == -2 && bcf_get_format_string(hdr, rec, "PS", &ps_str, &ps_str_cap) > 0) {   // declared as a String
        const std::string v = ps_str[sample_idx];
        if (!v.empty() && v != ".") r.set = named_set(v);
      }
      // What the genotype states: a phased one its alleles in order and its phase set; an unphased
      // one its alleles in no order (0/1 and 1/0 say the same, and a PS next to a slash says nothing).
      if (slash) std::sort(said.begin(), said.end());
      for (std::size_t i = 0; i < said.size(); ++i) { if (i) r.statement += slash ? '/' : '|'; r.statement += said[i]; }
      if (!slash) r.statement += ":" + std::to_string(r.set);
    }
    raws.push_back(std::move(r));
  }
  if (read_failed()) {
    w.lo = w.hi = -1;   // not a window to serve again
    throw std::runtime_error("error reading --haplotype-vcf " + path + " at " + chrom + ":" + std::to_string(lo + 1) +
                             " (" + bcf_sr_strerror(sr->errnum) + ")");
  }
  return raws;
}

VcfLocal VcfHaplotypes::local(const std::string& chrom, int64_t lo, int64_t hi, const std::string& ref_window,
                              const FastaFetcher* reference, int toward, const std::vector<std::string>* flipped) const {
  VcfLocal L;
  const int rid = bcf_hdr_name2id(impl_->hdr, chrom.c_str());
  if (rid < 0) return L;
  L.known_contig = true;
  L.lo = lo; L.hi = hi; L.ref = ref_window;
  const std::size_t nh = n_haplotypes();
  L.edits.assign(nh, {});
  if (hi <= lo || static_cast<int64_t>(ref_window.size()) != hi - lo) { L.hi = lo; L.ref.clear(); return L; }
  const std::vector<RawRecord>& raws = impl_->window(chrom, rid, lo, hi, nh);

  // ---- a site written more than once is one statement, or a contradiction
  std::unordered_map<std::string, std::size_t> first_of;
  std::vector<char> drop(raws.size(), 0), contradicted(raws.size(), 0);
  for (std::size_t i = 0; i < raws.size(); ++i) {
    const auto it = first_of.find(raws[i].key);
    if (it == first_of.end()) { first_of.emplace(raws[i].key, i); continue; }
    drop[i] = 1;
    if (raws[it->second].statement != raws[i].statement) contradicted[it->second] = 1;
  }

  // ---- every allele as an edit of the reference
  for (std::size_t i = 0; i < raws.size(); ++i) {
    if (drop[i]) continue;
    const RawRecord& r = raws[i];
    VcfLocal::Rec rec;
    rec.id = r.key + "#" + r.statement;
    rec.phased = r.phased;
    rec.set = r.set;
    const int ri = static_cast<int>(L.recs.size());
    const std::string& ref = r.alleles[0];
    bool unknown_all = false;
    if (contradicted[i]) { rec.conflict = true; unknown_all = true; }
    else if (!r.has_gt) { rec.unknown_slots = static_cast<int>(nh); unknown_all = true; }
    // a REF that is not the reference: the record describes another genome
    for (int64_t p = std::max(r.pos, lo); p < std::min<int64_t>(r.pos + static_cast<int64_t>(ref.size()), hi) && !unknown_all; ++p) {
      const char a = ref[static_cast<std::size_t>(p - r.pos)], b = up(L.ref[static_cast<std::size_t>(p - lo)]);
      if (acgt(a) && acgt(b) && a != b) { rec.ref_mismatch = true; unknown_all = true; }
    }
    auto unknown = [&](std::size_t h) { L.edits[h].push_back(VcfLocal::Edit{r.pos, r.end, "", ri, 'u', {}, false, false}); };
    if (unknown_all) {
      for (std::size_t h = 0; h < nh; ++h) unknown(h);
      L.recs.push_back(std::move(rec));
      continue;
    }
    rec.index_dependent = r.n_known >= 1 && (r.het || r.n_unknown >= 1);
    rec.unknown_slots = r.n_unknown;
    // the records that are placed together: a phase set, or one unphased record on its own
    if (rec.index_dependent) rec.group = r.phased ? "PS:" + std::to_string(r.set) : "U:" + rec.id;
    const bool flip = flipped != nullptr && rec.index_dependent &&
                      std::find(flipped->begin(), flipped->end(), rec.group) != flipped->end();
    bool simple = ref.size() == 1;
    for (std::size_t h = 0; h < nh; ++h) {
      const int a = r.allele[flip ? nh - 1 - h : h];
      if (a < 0) { unknown(h); continue; }
      if (a == 0) continue;
      const std::string& alt = r.alleles[static_cast<std::size_t>(a)];
      if (alt == "*") {
        // removed by a deletion that another record describes; when the file does not carry that
        // record, what stands here is not known
        simple = false;
        L.edits[h].push_back(VcfLocal::Edit{r.pos, r.pos + static_cast<int64_t>(ref.size()), "", ri, 'u', {}, false, true});
        continue;
      }
      if (alt.empty() || alt[0] == '<' || alt.find('[') != std::string::npos || alt.find(']') != std::string::npos || alt[0] == '.' || alt.back() == '.') {
        simple = false;
        if (alt.rfind("<DEL", 0) == 0 && r.sized) {
          L.edits[h].push_back(VcfLocal::Edit{r.pos + 1, r.end, "", ri, 'd', {}, false, false});
          rec.indel = true;
        } else {
          ++rec.skipped;   // an allele without a sequence: what stands here is not known
          unknown(h);
        }
        continue;
      }
      if (alt.size() == ref.size()) {
        if (alt.size() != 1 || !acgt(alt[0])) simple = false;
        L.edits[h].push_back(VcfLocal::Edit{r.pos, r.pos + static_cast<int64_t>(ref.size()), alt, ri, 's', {}, false, false});
        continue;
      }
      simple = false;
      std::size_t k = 0;
      while (k < ref.size() && k < alt.size() && ref[k] == alt[k]) ++k;
      std::string rr = ref.substr(k), aa = alt.substr(k);
      std::size_t t = 0;
      while (t < rr.size() && t < aa.size() && rr[rr.size() - 1 - t] == aa[aa.size() - 1 - t]) ++t;
      rr.resize(rr.size() - t); aa.resize(aa.size() - t);
      const int64_t q0 = r.pos + static_cast<int64_t>(k);
      L.edits[h].push_back(VcfLocal::Edit{q0, q0 + static_cast<int64_t>(rr.size()), aa, ri, rr.empty() ? 'i' : aa.empty() ? 'd' : 'c', {}, false, false});
      rec.indel = true;
    }
    rec.het_snv = simple && r.het && r.n_unknown == 0;
    L.recs.push_back(std::move(rec));
  }

  auto refbase = [&](int64_t p) -> char {
    if (p < 0) return 'N';
    if (p >= lo && p < hi) return up(L.ref[static_cast<std::size_t>(p - lo)]);
    if (reference == nullptr) return 'N';
    const std::string b = reference->fetch(chrom, p, p + 1);
    return b.empty() ? 'N' : b[0];
  };

  for (std::size_t h = 0; h < nh; ++h) {
    std::vector<VcfLocal::Edit>& ed = L.edits[h];
    auto by_start = [](const VcfLocal::Edit& a, const VcfLocal::Edit& b) { return a.q0 != b.q0 ? a.q0 < b.q0 : a.q1 < b.q1; };
    std::stable_sort(ed.begin(), ed.end(), by_start);
    // ---- two insertions at one point: the same bases twice are one statement; different bases
    // leave what stands there unknown, and with it the base on either side
    {
      std::vector<VcfLocal::Edit> one;
      bool inserts = false;   // one.back() is what the insertions at its point came to
      int64_t point = 0;
      for (VcfLocal::Edit& e : ed) {
        if (e.kind == 'i' && inserts && point == e.q0) {
          VcfLocal::Edit& p = one.back();
          p.also.push_back(e.rec);
          if (p.kind == 'i' && p.ins == e.ins) continue;
          p.kind = 'u'; p.clash = true; p.ins.clear();
          p.q0 = std::max<int64_t>(0, point - 1); p.q1 = point + 1;
          continue;
        }
        inserts = e.kind == 'i';
        point = e.q0;
        one.push_back(std::move(e));
      }
      ed.swap(one);
      std::stable_sort(ed.begin(), ed.end(), by_start);   // an edit made unknown starts one base earlier
    }
    // ---- two records on one base: the same substitution twice is one, anything else is unknown
    std::vector<VcfLocal::Edit> out;
    for (VcfLocal::Edit& e : ed) {
      if (!out.empty()) {
        VcfLocal::Edit& p = out.back();
        if (e.q0 < p.q1) {
          // a `*` allele inside the deletion or replacement that removes its bases: that record speaks
          if (p.star != e.star) {
            const VcfLocal::Edit& star = p.star ? p : e;
            const VcfLocal::Edit& other = p.star ? e : p;
            if ((other.kind == 'd' || other.kind == 'c') && other.q0 <= star.q0 && star.q1 <= other.q1) {
              if (p.star) p = e;
              continue;
            }
          }
          bool agree = p.kind == 's' && e.kind == 's';
          for (int64_t x = e.q0; agree && x < std::min(p.q1, e.q1); ++x)
            if (p.ins[static_cast<std::size_t>(x - p.q0)] != e.ins[static_cast<std::size_t>(x - e.q0)]) agree = false;
          p.also.push_back(e.rec);
          p.also.insert(p.also.end(), e.also.begin(), e.also.end());
          if (agree) {
            if (e.q1 > p.q1) { p.ins += e.ins.substr(static_cast<std::size_t>(p.q1 - e.q0)); p.q1 = e.q1; }
          } else {
            if (p.kind != 'u' && e.kind != 'u') p.clash = true;
            p.clash = p.clash || e.clash;
            p.kind = 'u'; p.ins.clear(); p.star = false;
            p.q1 = std::max(p.q1, e.q1);
          }
          continue;
        }
      }
      out.push_back(std::move(e));
    }
    ed.swap(out);
    // ---- insertions and deletions move towards the intron as far as the reference allows
    auto touched = [&](std::size_t self, int64_t p) {
      for (std::size_t j = 0; j < ed.size(); ++j) {
        if (j == self) continue;
        const VcfLocal::Edit& o = ed[j];
        if (o.q0 == o.q1 ? (o.q0 == p || o.q0 == p + 1) : (o.q0 <= p && p < o.q1)) return true;
      }
      return false;
    };
    for (std::size_t i = 0; i < ed.size() && toward != 0; ++i) {
      VcfLocal::Edit& e = ed[i];
      if (e.kind != 'i' && e.kind != 'd') continue;
      for (int step = 0; step < kMaxShift; ++step) {
        if (toward > 0) {
          const char a = refbase(e.q0);
          if (!acgt(a) || touched(i, e.q0)) break;
          if (e.kind == 'd') { if (a != refbase(e.q1) || touched(i, e.q1)) break; }
          else { if (a != e.ins.front()) break; e.ins = e.ins.substr(1) + e.ins.substr(0, 1); }
          ++e.q0; ++e.q1;
        } else {
          if (e.q0 == 0) break;
          const char a = refbase(e.q0 - 1);
          if (!acgt(a) || touched(i, e.q0 - 1)) break;
          if (e.kind == 'd') { if (a != refbase(e.q1 - 1) || touched(i, e.q1 - 1)) break; }
          else { if (a != e.ins.back()) break; e.ins = e.ins.substr(e.ins.size() - 1) + e.ins.substr(0, e.ins.size() - 1); }
          --e.q0; --e.q1;
        }
      }
    }
    std::stable_sort(ed.begin(), ed.end(), by_start);
  }
  // a record whose edit became part of a clash is a contradiction of the file; a `*` left on its
  // own is an allele that could not be applied
  for (std::size_t h = 0; h < nh; ++h)
    for (const VcfLocal::Edit& e : L.edits[h]) {
      if (e.clash) L.recs[static_cast<std::size_t>(e.rec)].conflict = true;
      if (e.star) ++L.recs[static_cast<std::size_t>(e.rec)].skipped;
    }
  return L;
}

std::vector<std::string> VcfHaplotypes::overlay(const std::string& chrom, int64_t start, int64_t end,
                                                const std::string& ref_window, WindowPhase& phase) const {
  if (bcf_hdr_name2id(impl_->hdr, chrom.c_str()) < 0) return {};
  std::vector<std::string> haps(n_haplotypes(), ref_window);
  if (end <= start || static_cast<int64_t>(ref_window.size()) != end - start) return haps;
  const VcfLocal L = local(chrom, start, end, ref_window, nullptr, 0);   // base by base: no indel has to move
  std::vector<int> used;
  for (std::size_t h = 0; h < haps.size(); ++h) {
    for (const VcfLocal::Edit& e : L.edits[h]) {
      const int64_t a = std::max(e.q0, start), b = std::min(e.q1, end);
      if (e.kind == 'i') continue;   // between two bases: no base of the window changes
      if (b <= a) continue;
      bool changed = false;   // an allele that leaves the window's bases as the reference has them did not shape it
      for (int64_t p = a; p < b; ++p) {
        const char c = e.kind == 's' ? e.ins[static_cast<std::size_t>(p - e.q0)] : 'N';
        if (c != up(ref_window[static_cast<std::size_t>(p - start)])) changed = true;
        haps[h][static_cast<std::size_t>(p - start)] = c;
      }
      if (!changed) continue;
      used.push_back(e.rec);
      used.insert(used.end(), e.also.begin(), e.also.end());
    }
  }
  std::sort(used.begin(), used.end());
  used.erase(std::unique(used.begin(), used.end()), used.end());
  for (const int r : used) {
    const VcfLocal::Rec& rec = L.recs[static_cast<std::size_t>(r)];
    if (rec.unknown_slots > 0 || rec.conflict || rec.ref_mismatch) phase.missing = true;
    if (!rec.index_dependent) continue;
    if (std::find(phase.counted.begin(), phase.counted.end(), rec.id) != phase.counted.end()) continue;
    phase.counted.push_back(rec.id);   // a record spanning both windows of a junction is read twice but is one record
    ++phase.n_het;
    if (!rec.phased) phase.unphased = true; else phase.add_set(rec.set);
  }
  count_used(L, used);
  return haps;
}

void VcfHaplotypes::snv_sites(const std::string& chrom, int64_t lo, int64_t hi, const std::string& ref_window, long set,
                              std::vector<PhaseSite>& sites, std::vector<std::pair<int64_t, int64_t>>& blocked) const {
  const VcfLocal L = local(chrom, lo, hi, ref_window, nullptr, 0);
  if (!L.known_contig || L.ref.empty()) return;
  const std::size_t nh = n_haplotypes();
  std::map<int64_t, std::string> at;   // position -> one base per haplotype
  for (std::size_t h = 0; h < nh; ++h) {
    for (const VcfLocal::Edit& e : L.edits[h]) {
      const VcfLocal::Rec& rec = L.recs[static_cast<std::size_t>(e.rec)];
      const bool single = e.kind == 's' && e.q1 - e.q0 == 1 && acgt(e.ins[0]) && e.also.empty() && !e.clash;
      // a heterozygous site counts only when its alleles are phased in the wanted set
      const bool linked = !rec.index_dependent || (rec.het_snv && rec.phased && set != kNoPhaseSet && rec.set == set);
      if (!single || !linked || rec.unknown_slots > 0 || rec.indel || rec.skipped > 0) {
        blocked.emplace_back(e.q0 > 0 && e.q0 == e.q1 ? e.q0 - 1 : e.q0, std::max(e.q1, e.q0 + 1));
        continue;
      }
      if (e.q0 < lo || e.q0 >= hi) continue;
      auto it = at.find(e.q0);
      if (it == at.end()) it = at.emplace(e.q0, std::string(nh, up(L.ref[static_cast<std::size_t>(e.q0 - lo)]))).first;
      it->second[h] = e.ins[0];
    }
  }
  for (const auto& kv : at) {
    PhaseSite s;
    s.pos = kv.first;
    s.bases = kv.second;
    sites.push_back(std::move(s));
  }
}

}  // namespace panisoguard
