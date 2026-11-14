#include "panisoguard/classify.hpp"

#include <cstdio>

#include "panisoguard/fingerprint.hpp"

namespace panisoguard {

namespace {

std::string category_of(const CallerIsoform& iso, const Catalog& cat) {
  const IntronChain& c = iso.chain;
  if (c.monoexonic()) {
    const std::string m = cat.mono_exon_category(c.chrom, c.strand, iso.start, iso.end);
    if (!m.empty()) return m;
    return cat.overlaps_any_strand(c.chrom, iso.start, iso.end) ? "antisense" : "intergenic";
  }
  if (cat.has_chain(fingerprint_intron_chain(c))) return "full-splice_match";
  if (cat.has_subchain(c, iso.start, iso.end)) return "incomplete-splice_match";
  if (!cat.overlaps(c.chrom, c.strand, iso.start, iso.end))
    return cat.overlaps_any_strand(c.chrom, iso.start, iso.end) ? "antisense" : "intergenic";
  bool all_sites_known = true;
  for (const auto& j : c.introns)
    if (!cat.has_intron_start(c.chrom, c.strand, j.start) || !cat.has_intron_end(c.chrom, c.strand, j.end)) { all_sites_known = false; break; }
  return all_sites_known ? "novel_in_catalog" : "novel_not_in_catalog";
}

std::string all_canonical_of(const CallerIsoform& iso, const FastaFetcher& ref) {
  if (iso.chain.monoexonic()) return "NA";
  for (const auto& j : iso.chain.introns) {
    const SpliceMotif m = splice_motif(ref, iso.chain.chrom, j, iso.chain.strand);
    if (!m.evaluable) return "NA";
    if (!m.canonical) return "non_canonical";
  }
  return "canonical";
}

// Percent A in the 20 genomic bases downstream of the 3' end, read on the transcript strand
// (SQANTI3's perc_A_downstream_TTS with its default window).
double perc_a_of(const CallerIsoform& iso, const FastaFetcher& ref) {
  constexpr int64_t kWindow = 20;
  if (iso.end <= iso.start || iso.chain.strand == Strand::kUnknown) return SqantiRecord::kNaN;
  std::string seq;
  char base;
  if (iso.chain.strand == Strand::kPlus) { seq = ref.fetch(iso.chain.chrom, iso.end, iso.end + kWindow); base = 'A'; }
  else { seq = ref.fetch(iso.chain.chrom, iso.start - kWindow, iso.start); base = 'T'; }
  if (seq.empty()) return SqantiRecord::kNaN;
  int n = 0;
  for (char ch : seq) if (ch == base) ++n;
  return 100.0 * n / static_cast<double>(seq.size());
}

}  // namespace

SqantiTable classify_isoforms(const std::vector<CallerIsoform>& isoforms, const Catalog& catalog,
                              const FastaFetcher* reference) {
  SqantiTable table;
  table.n_columns = 8;
  table.records.reserve(isoforms.size());
  for (const auto& iso : isoforms) {
    SqantiRecord r;
    r.isoform = iso.id;
    r.chrom = iso.chain.chrom;
    r.strand = iso.chain.strand;
    r.structural_category = category_of(iso, catalog);
    r.RTS_stage = "NA";
    r.polyA_motif_found = "NA";
    r.all_canonical = reference != nullptr ? all_canonical_of(iso, *reference) : "NA";
    r.perc_A_downstream_TTS = reference != nullptr ? perc_a_of(iso, *reference) : SqantiRecord::kNaN;
    table.records.push_back(std::move(r));
  }
  return table;
}

}  // namespace panisoguard
