#include "panisoguard/variant_motif.hpp"

#include <cctype>
#include <algorithm>
#include <cstring>
#include <stdexcept>

#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

namespace panisoguard {

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
  mutable long n_missing = 0;
  mutable long n_conflicting = 0;
  mutable int32_t* ps = nullptr;   // bcf_get_format_int32 buffer for PS
  mutable int ps_cap = 0;
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
}

VcfHaplotypes::~VcfHaplotypes() {
  if (impl_ != nullptr) {
    if (impl_->gt) free(impl_->gt);
    if (impl_->ps) free(impl_->ps);
    if (impl_->sr) bcf_sr_destroy(impl_->sr);
    delete impl_;
  }
}

std::size_t VcfHaplotypes::n_haplotypes() const { return static_cast<std::size_t>(impl_->ploidy); }
bool VcfHaplotypes::has_chrom(const std::string& chrom) const { return bcf_hdr_name2id(impl_->hdr, chrom.c_str()) >= 0; }
const std::string& VcfHaplotypes::path() const { return impl_->path; }
const std::string& VcfHaplotypes::sample() const { return impl_->sample; }
long VcfHaplotypes::n_unphased_het_applied() const { return impl_->n_unphased; }
long VcfHaplotypes::n_non_snv_skipped() const { return impl_->n_non_snv; }
long VcfHaplotypes::n_missing_alleles() const { return impl_->n_missing; }
long VcfHaplotypes::n_conflicting_records() const { return impl_->n_conflicting; }

std::vector<std::string> VcfHaplotypes::overlay(const std::string& chrom, int64_t start, int64_t end,
                                                const std::string& ref_window, WindowPhase& phase) const {
  const int rid = bcf_hdr_name2id(impl_->hdr, chrom.c_str());
  if (rid < 0) return {};
  std::vector<std::string> haps(n_haplotypes(), ref_window);
  if (end <= start || static_cast<int64_t>(ref_window.size()) != end - start) return haps;
  // A seek that fails names a contig the index has no records for: the window is reference.
  if (bcf_sr_seek(impl_->sr, chrom.c_str(), static_cast<hts_pos_t>(start > 0 ? start - 1 : 0)) < 0) return haps;
  while (bcf_sr_next_line(impl_->sr)) {
    bcf1_t* rec = bcf_sr_get_line(impl_->sr, 0);
    if (rec->rid != rid || rec->pos >= end) break;
    bcf_unpack(rec, BCF_UN_STR);
    const int reflen = static_cast<int>(std::strlen(rec->d.allele[0]));
    if (rec->pos + reflen <= start) continue;   // ends before the window
    // the window bases this record covers (an MNP may start before the window)
    const int64_t k0 = std::max<int64_t>(start, rec->pos), k1 = std::min<int64_t>(end, rec->pos + reflen);
    const int ngt = bcf_get_genotypes(impl_->hdr, rec, &impl_->gt, &impl_->gt_cap);
    if (ngt <= 0) {
      // a record with no GT says a variant may be here but not which haplotype carries it: unknown bases
      for (auto& h : haps) for (int64_t pos = k0; pos < k1; ++pos) h[static_cast<std::size_t>(pos - start)] = 'N';
      if (k1 > k0) { phase.missing = true; impl_->n_missing += static_cast<long>(haps.size()); }
      continue;
    }
    const int per = ngt / impl_->n_samples;
    const int* g = impl_->gt + impl_->sample_idx * per;
    // the record's identity: position and alleles; its genotype (alleles, separators, phase set) apart
    std::string record_key = std::to_string(rec->pos);
    for (int a = 0; a < rec->n_allele; ++a) { record_key += a ? ',' : ':'; record_key += rec->d.allele[a]; }
    std::string gt_key;
    for (int h = 0; h < per; ++h) {
      if (g[h] == bcf_int32_vector_end) break;
      if (h) gt_key += bcf_gt_is_phased(g[h]) ? '|' : '/';
      gt_key += bcf_gt_is_missing(g[h]) ? "." : std::to_string(bcf_gt_allele(g[h]));
    }
    {
      const int nps = bcf_get_format_int32(impl_->hdr, rec, "PS", &impl_->ps, &impl_->ps_cap);
      if (nps > 0) {
        const int32_t v = impl_->ps[impl_->sample_idx * (nps / impl_->n_samples)];
        if (v != bcf_int32_missing && v != bcf_int32_vector_end) gt_key += ":" + std::to_string(v);
      }
    }
    bool conflict = false, duplicate = false;
    for (const auto& s : phase.seen)
      if (s.first == record_key) { if (s.second != gt_key) conflict = true; else duplicate = true; }
    if (conflict) {
      // the same site with another genotype: the file contradicts itself here, so the bases are unknown
      for (auto& h : haps) for (int64_t pos = k0; pos < k1; ++pos) h[static_cast<std::size_t>(pos - start)] = 'N';
      if (k1 > k0) { phase.missing = true; ++impl_->n_conflicting; }
      continue;
    }
    if (!duplicate) phase.seen.emplace_back(record_key, gt_key);
    bool applied = false, het = false, phased = true, ended = false, wrote_n = false;
    int first_allele = -1, n_known = 0, n_missing_here = 0;
    for (int h = 0; h < static_cast<int>(haps.size()); ++h) {
      if (h >= per || g[h] == bcf_int32_vector_end) ended = true;
      if (!ended && h > 0 && !bcf_gt_is_phased(g[h])) phased = false;   // the separator before a missing allele counts too
      if (ended || bcf_gt_is_missing(g[h])) {
        // no allele for this haplotype at a motif base: unknown, not reference
        ++n_missing_here;
        for (int64_t pos = k0; pos < k1; ++pos) { haps[h][static_cast<std::size_t>(pos - start)] = 'N'; wrote_n = true; }
        if (k1 > k0) { phase.missing = true; ++impl_->n_missing; }
        continue;
      }
      ++n_known;
      const int a = bcf_gt_allele(g[h]);
      if (first_allele < 0) first_allele = a; else if (a != first_allele) het = true;
      if (a <= 0 || a >= rec->n_allele) continue;
      const char* alt = rec->d.allele[a];
      const int altlen = static_cast<int>(std::strlen(alt));
      if (altlen != reflen || alt[0] == '<' || alt[0] == '*') { ++impl_->n_non_snv; continue; }
      for (int64_t pos = k0; pos < k1; ++pos) {
        haps[h][static_cast<std::size_t>(pos - start)] = static_cast<char>(std::toupper(static_cast<unsigned char>(alt[pos - rec->pos])));
        applied = true;
      }
    }
    // A record whose haplotypes did not all receive the same known allele (heterozygous, or partly
    // missing) places its alleles by haplotype index; when two such records shape one motif, the
    // composition is an observation only if they are phased together. A record that changed no base
    // and left no base unknown (a skipped indel) places nothing.
    const bool index_dependent = n_known >= 1 && (het || n_missing_here >= 1);
    const std::string full_key = record_key + "#" + gt_key;
    const bool counted_before = std::find(phase.counted.begin(), phase.counted.end(), full_key) != phase.counted.end();
    if (index_dependent && (applied || wrote_n) && !counted_before) {
      phase.counted.push_back(full_key);   // a record spanning both windows of a junction is read twice but is one record
      ++phase.n_het;
      if (!phased) { ++impl_->n_unphased; phase.unphased = true; }
      else {
        // the phase set of this genotype: FORMAT PS, else the file's default set
        long set = kDefaultPhaseSet;
        const int nps = bcf_get_format_int32(impl_->hdr, rec, "PS", &impl_->ps, &impl_->ps_cap);
        if (nps > 0) {
          const int per_ps = nps / impl_->n_samples;
          const int32_t v = impl_->ps[impl_->sample_idx * per_ps];
          if (v != bcf_int32_missing && v != bcf_int32_vector_end) set = v;
        }
        phase.add_set(set);
      }
    }
  }
  return haps;
}

}  // namespace panisoguard
