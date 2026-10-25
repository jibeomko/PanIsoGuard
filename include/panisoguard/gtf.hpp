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
  std::string gene_id;  // GTF gene_id ("" when the source has none, e.g. BED12)
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

 private:
  std::unordered_set<Fingerprint> chains_;
  std::unordered_set<std::string> introns_;
};

Catalog build_catalog_from_gtf(const std::string& path);

}  // namespace panisoguard
