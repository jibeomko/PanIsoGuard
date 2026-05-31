#include "panisoguard/adjudicator.hpp"

#include <stdexcept>

namespace panisoguard {
namespace {

// Assemble the evidence for one SQANTI record.
EvidenceVector build_evidence(const SqantiRecord& r,
                              const std::map<std::string, IntronChain>& chains,
                              const Catalog* catalog, const SjTable* sj,
                              const RuleConfig& cfg) {
  EvidenceVector ev;
  ev.structural_category = r.structural_category;
  ev.is_novel = r.is_novel();

  // QC priors (consumed as-is).
  ev.rts_evaluable = !r.RTS_stage.empty();
  ev.rts_stage = (r.RTS_stage == "TRUE" || r.RTS_stage == "true" || r.RTS_stage == "True");
  ev.canon_evaluable = !r.all_canonical.empty();
  ev.noncanonical = (r.all_canonical == "non_canonical");
  ev.percA_evaluable = !sqanti_is_na(r.perc_A_downstream_TTS);
  ev.perc_A_downstream_TTS = r.perc_A_downstream_TTS;

  // Short-read corroboration of novel junctions requires the caller chain AND a
  // catalog (to know which junctions are novel) AND an SJ table.
  auto it = chains.find(r.isoform);
  ev.chain_available = (it != chains.end());
  if (ev.chain_available && catalog != nullptr && sj != nullptr) {
    const IntronChain& chain = it->second;
    int n_novel = 0, n_supported = 0;
    for (const auto& intron : chain.introns) {
      if (catalog->has_intron(chain.chrom, chain.strand, intron)) continue;  // known junction
      ++n_novel;
      const SjRecord* hit = sj->find_exact(chain.chrom, intron);
      const bool supported = hit != nullptr && hit->n_uniq >= cfg.sj_min_uniq_reads &&
                             (!cfg.sj_require_canonical_motif || hit->canonical());
      if (supported) ++n_supported;
    }
    ev.sj_evaluable = true;
    ev.n_novel_junctions = n_novel;
    ev.n_novel_jx_sr_supported = n_supported;
  }
  return ev;
}

}  // namespace

std::vector<AdjudicationResult> adjudicate(const AdjudicateInputs& in, const RuleEngine& engine) {
  if (in.sqanti == nullptr || in.chains == nullptr) {
    throw std::invalid_argument("adjudicate: sqanti and chains inputs are required");
  }

  std::vector<AdjudicationResult> out;
  out.reserve(in.sqanti->records.size());
  for (const auto& r : in.sqanti->records) {
    AdjudicationResult res;
    res.isoform_id = r.isoform;
    res.chrom = r.chrom;
    res.strand = r.strand;
    res.structural_category = r.structural_category;
    res.evidence = build_evidence(r, *in.chains, in.catalog, in.sj, engine.config());
    res.verdict = engine.evaluate(res.evidence);
    out.push_back(std::move(res));
  }
  return out;
}

}  // namespace panisoguard
