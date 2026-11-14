#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "panisoguard/fingerprint.hpp"
#include "panisoguard/types.hpp"

namespace panisoguard {

// A transcript reconstructed from GTF exon lines.
struct Transcript {
  std::string id;
  IntronChain chain;
  std::string gene_id;  // GTF gene_id ("" when the source has none, e.g. BED12)
  int64_t start = 0;    // transcript span, 0-based half-open (first exon start .. last exon end)
  int64_t end = 0;
  std::vector<std::pair<int64_t, int64_t>> exons;  // 0-based half-open, sorted
};

// Parse a GTF into transcripts (groups `exon` lines by transcript_id, derives the
// intron chain). Used for both the reference annotation and caller GTFs.
std::vector<Transcript> read_gtf_transcripts(const std::string& path);

// Extract the (unquoted) value of an attribute from a GTF column-9 string,
// e.g. extract_gtf_attr(attrs, "transcript_id"). Returns "" if not found.
std::string extract_gtf_attr(const std::string& attrs, const std::string& key);

// Reference annotation model: the known intron chains (by fingerprint; `combine` marks a
// chain known or novel with it) and the known introns (the adjudicator classifies each
// junction of a caller chain as known or novel with it). Both strand-aware.
class Catalog {
 public:
  void add_transcript(const Transcript& t);

  bool has_chain(Fingerprint fp) const { return chains_.count(fp) != 0; }
  bool has_intron(const std::string& chrom, Strand strand, const Junction& j) const {
    return introns_.count(junction_key(chrom, strand, j)) != 0;
  }
  // An annotated splice site: a position where an annotated intron starts (resp. ends) on
  // that strand. SQANTI3's novel-in-catalog is "known sites, new combination".
  bool has_intron_start(const std::string& chrom, Strand strand, int64_t pos) const;
  bool has_intron_end(const std::string& chrom, Strand strand, int64_t pos) const;
  // An annotated transcript on the same chromosome and strand carries `chain` as a
  // contiguous part of its intron chain, and none of its other introns lies inside
  // [start, end) (a retained intron makes the isoform novel-in-catalog, as in SQANTI3).
  bool has_subchain(const IntronChain& chain, int64_t start, int64_t end) const;
  // Mono-exonic isoform [start, end) against the annotated exons on `strand`: "fsm" when
  // inside an exon of a mono-exonic transcript, "ism" inside an exon of a multi-exonic one,
  // "genic" overlapping an exon otherwise, "genic_intron" inside a transcript span only, "".
  std::string mono_exon_category(const std::string& chrom, Strand strand, int64_t start, int64_t end) const;
  // An annotated transcript span overlaps [start, end) on `strand`, or on either strand.
  bool overlaps(const std::string& chrom, Strand strand, int64_t start, int64_t end) const;
  bool overlaps_any_strand(const std::string& chrom, int64_t start, int64_t end) const;

 private:
  std::unordered_set<Fingerprint> chains_;
  std::unordered_set<std::string> introns_;
  std::unordered_set<std::string> intron_starts_, intron_ends_;           // "chrom|strand:pos"
  std::vector<Transcript> txs_;                                           // every annotated transcript
  std::unordered_map<std::string, std::vector<std::size_t>> txs_by_intron_;     // junction key -> txs_ indices
  std::unordered_map<std::string, std::vector<std::size_t>> txs_by_locus_;      // "chrom|strand" -> txs_ indices
};

Catalog build_catalog_from_gtf(const std::string& path);

}  // namespace panisoguard
