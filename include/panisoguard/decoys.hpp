#pragma once

#include <string>
#include <vector>

#include "panisoguard/gtf.hpp"
#include "panisoguard/types.hpp"

namespace panisoguard {

// A novel-in-catalog decoy: the 5' introns of one annotated transcript joined to the 3'
// introns of another transcript of the same gene at an intron the two share, kept only when
// the result is neither an annotated chain nor a contiguous piece of one. Every junction is
// known and only the combination is new: the class of chain that whole-chain read support
// lifts and that no short-read junction test can reach. Counting the reads that carry such
// chains in a BAM gives a per-dataset estimate of how often the lift fires on a chain the
// annotation does not know (an upper bound on its false-lift rate: some decoys are real,
// unannotated isoforms).
struct Decoy {
  std::string gene_id;
  std::string parent_a;  // supplies the introns up to and including the shared one
  std::string parent_b;  // supplies the introns after it
  IntronChain chain;
};

// Decoys from every ordered pair of multi-intron transcripts of one gene (same chrom and
// strand, same gene_id); at most `max_per_gene` per gene, sampled with `seed` when a gene
// yields more. Deterministic for a given input and seed.
std::vector<Decoy> make_decoys(const std::vector<Transcript>& transcripts, int max_per_gene,
                               unsigned seed);

}  // namespace panisoguard
