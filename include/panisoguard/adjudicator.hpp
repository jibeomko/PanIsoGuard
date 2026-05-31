#pragma once

#include <map>
#include <string>
#include <vector>

#include "panisoguard/bam_features.hpp"
#include "panisoguard/gtf.hpp"   // Catalog
#include "panisoguard/rules.hpp"
#include "panisoguard/sj_tab.hpp"
#include "panisoguard/sqanti.hpp"
#include "panisoguard/types.hpp"
#include "panisoguard/variant_motif.hpp"
#include "panisoguard/verdict.hpp"

namespace panisoguard {

struct AdjudicationResult {
  std::string isoform_id;
  std::string chrom;
  Strand strand = Strand::kUnknown;
  std::string structural_category;
  EvidenceVector evidence;
  Verdict verdict;
};

struct AdjudicateInputs {
  const SqantiTable* sqanti = nullptr;                       // required
  const std::map<std::string, IntronChain>* chains = nullptr;  // isoform_id -> chain (required)
  const Catalog* catalog = nullptr;                          // optional (needed for novel-junction axes)
  const SjTable* sj = nullptr;                               // optional (short-read corroboration)
  const BamReader* bam = nullptr;                            // optional (mapping axis)
  const HaplotypeProvider* haplotype = nullptr;              // optional (variant/reference-bias axis)
  bool haplotype_circular = false;                           // true if haplotype provenance is RNA-derived/unknown
};

// Adjudicate every SQANTI record: assemble its EvidenceVector (joining the caller
// chain, classifying novel junctions vs the catalog, and corroborating them
// against short-read SJ.tab) and apply the rule engine.
std::vector<AdjudicationResult> adjudicate(const AdjudicateInputs& in, const RuleEngine& engine);

}  // namespace panisoguard
