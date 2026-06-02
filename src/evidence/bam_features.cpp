#include "panisoguard/bam_features.hpp"

#include <cstdint>
#include <cstdlib>
#include <stdexcept>

#include <htslib/hts.h>
#include <htslib/sam.h>

namespace panisoguard {

struct BamReader::Impl {
  htsFile* fp = nullptr;
  sam_hdr_t* hdr = nullptr;
  hts_idx_t* idx = nullptr;
  std::string path;
};

BamReader::BamReader(const std::string& path, const std::string& reference) : impl_(new Impl) {
  impl_->path = path;
  impl_->fp = sam_open(path.c_str(), "r");
  if (impl_->fp == nullptr) {
    delete impl_;
    throw std::runtime_error("cannot open BAM/CRAM: " + path);
  }
  if (!reference.empty()) {
    hts_set_fai_filename(impl_->fp, reference.c_str());  // needed to decode CRAM
  }
  impl_->hdr = sam_hdr_read(impl_->fp);
  if (impl_->hdr == nullptr) {
    sam_close(impl_->fp);
    delete impl_;
    throw std::runtime_error("cannot read BAM/CRAM header: " + path);
  }
  impl_->idx = sam_index_load(impl_->fp, path.c_str());
  if (impl_->idx == nullptr) {
    sam_hdr_destroy(impl_->hdr);
    sam_close(impl_->fp);
    delete impl_;
    throw std::runtime_error("cannot load BAM/CRAM index (.bai/.csi) for: " + path);
  }
}

BamReader::~BamReader() {
  if (impl_ != nullptr) {
    if (impl_->idx) hts_idx_destroy(impl_->idx);
    if (impl_->hdr) sam_hdr_destroy(impl_->hdr);
    if (impl_->fp) sam_close(impl_->fp);
    delete impl_;
  }
}

namespace {
inline bool near(int64_t pos, const Junction& j, int w) {
  return (pos >= j.start - w && pos <= j.start + w) || (pos >= j.end - w && pos <= j.end + w);
}
}  // namespace

JunctionBamFeatures BamReader::features_at_junction(const std::string& chrom,
                                                    const Junction& intron,
                                                    const BamFeatureParams& p) const {
  JunctionBamFeatures f;
  const int tid = sam_hdr_name2tid(impl_->hdr, chrom.c_str());
  if (tid < 0) return f;  // contig absent -> not evaluated
  f.evaluated = true;

  const int64_t beg = intron.start - p.junction_window_bp > 0 ? intron.start - p.junction_window_bp : 0;
  const int64_t end = intron.end + p.junction_window_bp;
  hts_itr_t* itr = sam_itr_queryi(impl_->idx, tid, beg, end);
  if (itr == nullptr) return f;

  bam1_t* b = bam_init1();
  if (b == nullptr) {  // allocation failure: clean up the iterator, do not dereference null
    hts_itr_destroy(itr);
    return f;
  }
  while (sam_itr_next(impl_->fp, itr, b) >= 0) {
    if (b->core.flag & BAM_FUNMAP) continue;

    const uint32_t* cig = bam_get_cigar(b);
    const int ncig = b->core.n_cigar;
    int64_t rpos = b->core.pos;  // 0-based leftmost reference coordinate
    bool spanning = false, softclip = false, indel_near = false;

    for (int i = 0; i < ncig; ++i) {
      const int op = bam_cigar_op(cig[i]);
      const int len = bam_cigar_oplen(cig[i]);
      switch (op) {
        case BAM_CSOFT_CLIP:
          if (len >= p.softclip_min_bp && (i == 0 || i == ncig - 1)) softclip = true;
          break;
        case BAM_CREF_SKIP:  // N: an intron
          if (rpos == intron.start && rpos + len == intron.end) spanning = true;
          rpos += len;
          break;
        case BAM_CDEL:
          if (near(rpos, intron, p.junction_window_bp)) indel_near = true;
          rpos += len;
          break;
        case BAM_CINS:
          if (near(rpos, intron, p.junction_window_bp)) indel_near = true;
          break;  // insertion consumes no reference
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
          rpos += len;
          break;
        default:  // hard clip / pad consume neither here
          break;
      }
    }

    if (!spanning) continue;
    ++f.n_spanning;
    if (b->core.qual < p.min_mapq) ++f.n_low_mapq;
    if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ++f.n_supplementary;
    if (softclip) ++f.n_softclip;
    if (indel_near) ++f.n_indel_near;
  }
  bam_destroy1(b);
  hts_itr_destroy(itr);
  return f;
}

}  // namespace panisoguard
