#include "panisoguard/adjudicator.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

namespace panisoguard {
namespace {

std::string jx_key(const std::string& chrom, Strand strand, const Junction& j) {
  // Strand-aware to match SjTable::key / Catalog: two opposite-strand isoforms
  // sharing identical intron coordinates (antisense/overlapping loci) must not
  // collide in the BAM-feature cache.
  return chrom + '|' + strand_char(strand) + ':' + std::to_string(j.start) + '-' +
         std::to_string(j.end);
}

// Assemble the evidence for one SQANTI record. `bam_cache` memoizes per-junction
// BAM features across isoforms that share a junction.
EvidenceVector build_evidence(const SqantiRecord& r,
                              const std::map<std::string, IntronChain>& chains,
                              const Catalog* catalog, const SjTable* sj, const BamReader* bam,
                              const HaplotypeProvider* haplotype, bool haplotype_circular,
                              const PangenomeJunctions* pangenome, bool pangenome_circular,
                              const RuleConfig& cfg,
                              std::unordered_map<std::string, JunctionBamFeatures>& bam_cache) {
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

  // Novel-junction axes (short-read corroboration + BAM mapping features) require
  // the caller chain AND a catalog to know which junctions are novel.
  auto it = chains.find(r.isoform);
  ev.chain_available = (it != chains.end());
  if (ev.chain_available && catalog != nullptr) {
    const IntronChain& chain = it->second;
    BamFeatureParams bp;
    bp.min_mapq = cfg.bam_min_mapq;
    bp.softclip_min_bp = cfg.bam_softclip_min_bp;
    bp.junction_window_bp = cfg.bam_junction_window_bp;

    int n_novel = 0, n_supported = 0, n_pan = 0, n_variant_created = 0;
    bool any_bam_eval = false;
    bool any_variant_eval = false;
    bool any_pan_eval = false;
    for (const auto& intron : chain.introns) {
      if (catalog->has_intron(chain.chrom, chain.strand, intron)) continue;  // known junction
      ++n_novel;
      if (haplotype != nullptr) {
        const VariantVerdict vv = haplotype->classify(chain.chrom, intron, chain.strand);
        if (vv != VariantVerdict::kNotEvaluable) any_variant_eval = true;
        if (vv == VariantVerdict::kCreated) ++n_variant_created;
      }
      if (pangenome != nullptr) {
        any_pan_eval = true;
        if (pangenome->supports(chain.chrom, chain.strand, intron, cfg.pangenome_min_haplotypes)) {
          ++n_pan;
        }
      }
      if (sj != nullptr) {
        const SjRecord* hit = sj->find_exact(chain.chrom, chain.strand, intron);
        const bool supported = hit != nullptr && hit->n_uniq >= cfg.sj_min_uniq_reads &&
                               (!cfg.sj_require_canonical_motif || hit->canonical());
        if (supported) ++n_supported;
      }
      if (bam != nullptr) {
        const std::string key = jx_key(chain.chrom, chain.strand, intron);
        auto cit = bam_cache.find(key);
        if (cit == bam_cache.end()) {
          cit = bam_cache.emplace(key, bam->features_at_junction(chain.chrom, intron, bp)).first;
        }
        const JunctionBamFeatures& f = cit->second;
        if (f.evaluated) {
          any_bam_eval = true;
          ev.bam_n_spanning_total += f.n_spanning;
          if (f.n_spanning > 0) {
            if (f.frac_low_mapq() > ev.bam_max_frac_low_mapq) ev.bam_max_frac_low_mapq = f.frac_low_mapq();
            if (f.frac_supplementary() > ev.bam_max_frac_supplementary) ev.bam_max_frac_supplementary = f.frac_supplementary();
            if (f.frac_softclip() > ev.bam_max_frac_softclip) ev.bam_max_frac_softclip = f.frac_softclip();
            if (f.frac_indel_near() > ev.bam_max_frac_indel_near) ev.bam_max_frac_indel_near = f.frac_indel_near();
          }
        }
      }
    }
    if (sj != nullptr) {
      ev.sj_evaluable = true;
      ev.n_novel_junctions = n_novel;
      ev.n_novel_jx_sr_supported = n_supported;
    } else {
      ev.n_novel_junctions = n_novel;  // recorded even without SR for reference
    }
    ev.bam_evaluable = (bam != nullptr) && any_bam_eval;
    ev.variant_evaluable = (haplotype != nullptr) && any_variant_eval;
    // Rescue requires ALL of the isoform's novel junctions to be explained by the
    // reference-bias evidence -- a single explained junction does not make a
    // multi-novel-junction isoform "false novel" (the others may be genuinely new).
    ev.variant_rescue = (n_novel > 0 && n_variant_created == n_novel);
    ev.variant_circular = haplotype_circular;
    ev.pangenome_evaluable = (pangenome != nullptr) && any_pan_eval;
    ev.n_novel_jx_pangenome = n_pan;
    ev.pangenome_rescue = (n_novel > 0 && n_pan == n_novel);
    ev.pangenome_circular = pangenome_circular;
  }
  return ev;
}

}  // namespace

std::vector<AdjudicationResult> adjudicate(const AdjudicateInputs& in, const RuleEngine& engine) {
  if (in.sqanti == nullptr || in.chains == nullptr) {
    throw std::invalid_argument("adjudicate: sqanti and chains inputs are required");
  }

  std::unordered_map<std::string, JunctionBamFeatures> bam_cache;
  std::vector<AdjudicationResult> out;
  out.reserve(in.sqanti->records.size());
  for (const auto& r : in.sqanti->records) {
    AdjudicationResult res;
    res.isoform_id = r.isoform;
    res.chrom = r.chrom;
    res.strand = r.strand;
    res.structural_category = r.structural_category;
    res.evidence = build_evidence(r, *in.chains, in.catalog, in.sj, in.bam, in.haplotype,
                                  in.haplotype_circular, in.pangenome, in.pangenome_circular,
                                  engine.config(), bam_cache);
    res.verdict = engine.evaluate(res.evidence);
    out.push_back(std::move(res));
  }
  return out;
}

}  // namespace panisoguard
