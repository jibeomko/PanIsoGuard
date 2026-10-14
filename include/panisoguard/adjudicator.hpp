#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "panisoguard/bam_features.hpp"
#include "panisoguard/consensus.hpp"  // CallerSupport
#include "panisoguard/gtf.hpp"   // Catalog
#include "panisoguard/pangenome.hpp"
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
  // SQANTI3 QC descriptors passed through verbatim (consumed as priors, NEVER
  // recomputed by PanIsoGuard) and emitted as `bio_flags`. They are biological-
  // plausibility annotations only and do NOT enter the rule engine / verdict.
  double bio_dist_to_CAGE_peak = SqantiRecord::kNaN;
  double bio_dist_to_polyA_site = SqantiRecord::kNaN;
  std::string bio_predicted_NMD;
  std::string bio_polyA_motif_found;
};

struct AdjudicateInputs {
  const SqantiTable* sqanti = nullptr;                       // required
  const std::map<std::string, IntronChain>* chains = nullptr;  // isoform_id -> chain (required)
  const Catalog* catalog = nullptr;                          // optional (needed for novel-junction axes)
  const SjTable* sj = nullptr;                               // optional (short-read corroboration)
  const BamReader* bam = nullptr;                            // optional (mapping axis)
  const HaplotypeProvider* haplotype = nullptr;              // optional (variant/reference-bias axis)
  bool haplotype_circular = false;                           // true if haplotype provenance is RNA-derived/unknown
  const PangenomeJunctions* pangenome = nullptr;             // optional (pangenome reference-bias axis)
  bool pangenome_circular = false;                           // true if junction-set provenance is sample-derived/unknown
  const CallerSupport* caller_support = nullptr;  // optional: `combine` matrix (consensus axis), joined by chain fingerprint
};

// One junction of the adjudicated isoforms with the per-junction evidence the isoform
// verdict aggregates (<prefix>.junctions.tsv). A reference-bias junction inside an isoform
// that was not rescued (another novel junction unexplained) stays visible here.
struct JunctionRow {
  std::string chrom;
  Junction intron;                       // 0-based half-open, as everywhere internally
  Strand strand = Strand::kUnknown;
  std::vector<std::string> isoforms;     // every adjudicated isoform using it
  bool known = false;                    // in the reference catalog
  bool sj_evaluable = false;  int sj_n_uniq = 0;  bool sr_supported = false;
  // novel junctions only (the axes below classify novelty vs the reference):
  bool variant_evaluable = false;  std::string ref_motif;  std::vector<std::string> hap_motifs;
  std::string variant = "NA";           // created | disrupted | none | NA
  std::vector<int> canonical_on;         // 1-based haplotype indices whose motif is canonical
  bool pangenome_evaluable = false;  bool pangenome_supported = false;
  int pangenome_n_haplotypes = 0;  std::vector<std::string> pangenome_paths;
  bool bam_evaluable = false;  JunctionBamFeatures bam;
};

// Adjudicate every SQANTI record: assemble its EvidenceVector (joining the caller
// chain, classifying novel junctions vs the catalog, and corroborating them
// against short-read SJ.tab) and apply the rule engine. If `junctions` is given it
// receives one row per distinct junction of the adjudicated isoforms, in genomic order.
std::vector<AdjudicationResult> adjudicate(const AdjudicateInputs& in, const RuleEngine& engine,
                                           std::vector<JunctionRow>* junctions = nullptr);

}  // namespace panisoguard
