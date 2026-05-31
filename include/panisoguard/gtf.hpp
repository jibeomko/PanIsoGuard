#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "panisoguard/fingerprint.hpp"
#include "panisoguard/types.hpp"

namespace panisoguard {

// A transcript reconstructed from GTF exon lines.
struct Transcript {
  std::string id;
  IntronChain chain;
};

// Parse a GTF into transcripts (groups `exon` lines by transcript_id, derives the
// intron chain). Used for both the reference annotation and caller GTFs.
std::vector<Transcript> read_gtf_transcripts(const std::string& path);

// Extract the (unquoted) value of an attribute from a GTF column-9 string,
// e.g. extract_gtf_attr(attrs, "transcript_id"). Returns "" if not found.
std::string extract_gtf_attr(const std::string& attrs, const std::string& key);

// Reference annotation model: known intron chains (by fingerprint), known
// introns, and known splice sites. Feeds the NIC (all sites known, novel
// combination) vs NNC (>=1 novel site) boundary and circularity detection.
class Catalog {
 public:
  void add_transcript(const Transcript& t);

  bool has_chain(Fingerprint fp) const { return chains_.count(fp) != 0; }
  bool has_intron(const std::string& chrom, Strand strand, const Junction& j) const {
    return introns_.count(intron_key(chrom, strand, j)) != 0;
  }
  bool has_site(const std::string& chrom, Strand strand, int64_t pos) const {
    return sites_.count(site_key(chrom, strand, pos)) != 0;
  }

  std::size_t n_chains() const { return chains_.size(); }
  std::size_t n_introns() const { return introns_.size(); }
  std::size_t n_sites() const { return sites_.size(); }

 private:
  static std::string intron_key(const std::string& chrom, Strand strand, const Junction& j);
  static std::string site_key(const std::string& chrom, Strand strand, int64_t pos);

  std::unordered_set<Fingerprint> chains_;
  std::unordered_set<std::string> introns_;
  std::unordered_set<std::string> sites_;
};

Catalog build_catalog_from_gtf(const std::string& path);

}  // namespace panisoguard
