#include "panisoguard/decoys.hpp"

#include <algorithm>
#include <map>
#include <random>
#include <unordered_set>

#include "panisoguard/fingerprint.hpp"

namespace panisoguard {

std::vector<Decoy> make_decoys(const std::vector<Transcript>& transcripts, int max_per_gene,
                               unsigned seed) {
  // std::map keeps gene order deterministic across runs.
  std::map<std::string, std::vector<const Transcript*>> genes;
  for (const auto& t : transcripts) {
    if (t.gene_id.empty() || t.chain.introns.size() < 2) continue;
    genes[t.chain.chrom + '|' + strand_char(t.chain.strand) + '|' + t.gene_id].push_back(&t);
  }
  std::mt19937 rng(seed);
  std::vector<Decoy> out;
  for (const auto& kv : genes) {
    const auto& g = kv.second;
    if (g.size() < 2) continue;
    std::vector<Decoy> cand;
    std::unordered_set<Fingerprint> seen;
    for (const Transcript* a : g) {
      for (const Transcript* b : g) {
        if (a == b) continue;
        const auto& ia = a->chain.introns;
        const auto& ib = b->chain.introns;
        for (std::size_t i = 1; i < ia.size(); ++i) {           // i == 0 would give a suffix of b
          for (std::size_t j = 0; j + 1 < ib.size(); ++j) {     // j == last would give a prefix of a
            if (!(ia[i] == ib[j])) continue;
            Decoy d;
            d.gene_id = a->gene_id;
            d.parent_a = a->id;
            d.parent_b = b->id;
            d.chain.chrom = a->chain.chrom;
            d.chain.strand = a->chain.strand;
            d.chain.introns.assign(ia.begin(), ia.begin() + static_cast<std::ptrdiff_t>(i) + 1);
            d.chain.introns.insert(d.chain.introns.end(), ib.begin() + static_cast<std::ptrdiff_t>(j) + 1, ib.end());
            const Fingerprint fp = fingerprint_intron_chain(d.chain);
            if (seen.count(fp)) continue;
            bool piece_of_known = false;   // annotated, or a contiguous part of an annotated chain
            for (const Transcript* t : g) {
              if (contiguous_subchain(d.chain.introns, t->chain.introns)) { piece_of_known = true; break; }
            }
            if (piece_of_known) continue;
            seen.insert(fp);
            cand.push_back(std::move(d));
          }
        }
      }
    }
    if (max_per_gene > 0 && cand.size() > static_cast<std::size_t>(max_per_gene)) {
      std::shuffle(cand.begin(), cand.end(), rng);
      cand.resize(static_cast<std::size_t>(max_per_gene));
    }
    for (auto& d : cand) out.push_back(std::move(d));
  }
  return out;
}

}  // namespace panisoguard
