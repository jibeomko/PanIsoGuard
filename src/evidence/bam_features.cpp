#include "panisoguard/bam_features.hpp"

#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

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

// The transcript strand the aligner assigned to a spliced read. STAR and HISAT2 write the
// genomic strand as XS:A; minimap2 writes ts:A relative to the read, so a reverse-strand
// alignment flips it. No tag: unknown.
Strand read_strand(const bam1_t* b) {
  if (const uint8_t* xs = bam_aux_get(b, "XS"); xs != nullptr && *xs == 'A') return parse_strand(bam_aux2A(xs));
  if (const uint8_t* ts = bam_aux_get(b, "ts"); ts != nullptr && *ts == 'A') {
    const char t = bam_aux2A(ts);
    if (t != '+' && t != '-') return Strand::kUnknown;
    const bool reverse = (b->core.flag & BAM_FREVERSE) != 0;
    return ((t == '+') != reverse) ? Strand::kPlus : Strand::kMinus;
  }
  return Strand::kUnknown;
}

// A read counts for `wanted` unless the aligner put it on the other strand.
bool strand_compatible(const bam1_t* b, Strand wanted) {
  if (wanted == Strand::kUnknown) return true;
  const Strand s = read_strand(b);
  return s == Strand::kUnknown || s == wanted;
}

hts_itr_t* query(hts_idx_t* idx, const std::string& path, int tid, int64_t beg, int64_t end) {
  hts_itr_t* itr = sam_itr_queryi(idx, tid, beg, end);
  if (itr == nullptr) throw std::runtime_error("cannot query the index of " + path);
  return itr;
}

bam1_t* new_record() {
  bam1_t* b = bam_init1();
  if (b == nullptr) throw std::runtime_error("out of memory reading a BAM record");
  return b;
}

// After an iterator loop: -1 is the end of the region, anything lower an htslib read error.
void check_read(int ret, const std::string& path, const std::string& chrom) {
  if (ret < -1) throw std::runtime_error("error reading " + path + " at " + chrom + " (htslib code " + std::to_string(ret) + ")");
}

}  // namespace

JunctionBamFeatures BamReader::features_at_junction(const std::string& chrom, const Junction& intron,
                                                    Strand strand, const BamFeatureParams& p) const {
  JunctionBamFeatures f;
  const int tid = sam_hdr_name2tid(impl_->hdr, chrom.c_str());
  if (tid < 0) return f;  // contig absent -> not evaluated

  const int64_t beg = intron.start - p.junction_window_bp > 0 ? intron.start - p.junction_window_bp : 0;
  const int64_t end = intron.end + p.junction_window_bp;
  hts_itr_t* itr = query(impl_->idx, impl_->path, tid, beg, end);
  bam1_t* b = new_record();
  int ret;
  while ((ret = sam_itr_next(impl_->fp, itr, b)) >= 0) {
    if (b->core.flag & BAM_FUNMAP) continue;
    if (!strand_compatible(b, strand)) continue;

    const uint32_t* cig = bam_get_cigar(b);
    const int ncig = b->core.n_cigar;
    int64_t rpos = b->core.pos;  // 0-based leftmost reference coordinate
    bool spanning = false, indel_near = false;

    for (int i = 0; i < ncig; ++i) {
      const int op = bam_cigar_op(cig[i]);
      const int len = bam_cigar_oplen(cig[i]);
      switch (op) {
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
        default:  // soft/hard clip and pad consume no reference
          break;
      }
    }

    if (!spanning) continue;
    ++f.n_spanning;
    if (b->core.qual < p.min_mapq) ++f.n_low_mapq;
    if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ++f.n_supplementary;
    if (indel_near) ++f.n_indel_near;
  }
  bam_destroy1(b);
  hts_itr_destroy(itr);
  check_read(ret, impl_->path, chrom);
  f.evaluated = true;
  return f;
}

ChainBamSupport BamReader::chain_support(const IntronChain& chain) const {
  ChainBamSupport s;
  if (chain.monoexonic()) return s;
  const int tid = sam_hdr_name2tid(impl_->hdr, chain.chrom.c_str());
  if (tid < 0) return s;
  // A read carrying the chain aligns across the first intron's donor base, so fetching the
  // reads at that one position is enough (and far cheaper than the whole isoform span in a
  // deep, isoform-rich locus). n_partial is counted over the same reads.
  const int64_t beg = chain.introns.front().start > 0 ? chain.introns.front().start - 1 : 0;
  const int64_t end = chain.introns.front().start + 1;
  hts_itr_t* itr = query(impl_->idx, impl_->path, tid, beg, end);
  bam1_t* b = new_record();
  std::vector<Junction> introns;
  int ret;
  while ((ret = sam_itr_next(impl_->fp, itr, b)) >= 0) {
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
    if (!strand_compatible(b, chain.strand)) continue;
    introns.clear();
    const uint32_t* cig = bam_get_cigar(b);
    int64_t rpos = b->core.pos;
    for (int i = 0; i < static_cast<int>(b->core.n_cigar); ++i) {
      const int op = bam_cigar_op(cig[i]);
      const int len = bam_cigar_oplen(cig[i]);
      if (op == BAM_CREF_SKIP) { introns.push_back(Junction{rpos, rpos + len}); rpos += len; }
      else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) rpos += len;
    }
    if (introns.empty()) continue;
    if (contiguous_subchain(chain.introns, introns)) ++s.n_full;
    else if (contiguous_subchain(introns, chain.introns)) ++s.n_partial;
  }
  bam_destroy1(b);
  hts_itr_destroy(itr);
  check_read(ret, impl_->path, chain.chrom);
  s.evaluated = true;
  return s;
}

}  // namespace panisoguard
