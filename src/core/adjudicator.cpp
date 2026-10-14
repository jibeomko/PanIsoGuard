#include "panisoguard/adjudicator.hpp"
#include "panisoguard/fingerprint.hpp"

#include <algorithm>

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
                              const CallerSupport* caller_support,
                              const RuleConfig& cfg,
                              std::unordered_map<std::string, JunctionBamFeatures>& bam_cache) {
  EvidenceVector ev;
  ev.structural_category = r.structural_category;
  ev.is_novel = r.is_novel();

  // QC priors (consumed as-is). An NA column is not_evaluable, never "false".
  ev.rts_evaluable = !sqanti_is_na_token(r.RTS_stage);
  ev.rts_stage = (r.RTS_stage == "TRUE" || r.RTS_stage == "true" || r.RTS_stage == "True");
  ev.canon_evaluable = !sqanti_is_na_token(r.all_canonical);
  ev.noncanonical = (r.all_canonical == "non_canonical");
  ev.percA_evaluable = !sqanti_is_na(r.perc_A_downstream_TTS);
  ev.perc_A_downstream_TTS = r.perc_A_downstream_TTS;
  // Present only when SQANTI3 ran with --polyA_motif_list; refines the intra-priming sign.
  ev.polyA_motif_evaluable = !sqanti_is_na_token(r.polyA_motif_found);
  ev.polyA_motif_found = (r.polyA_motif_found == "TRUE" || r.polyA_motif_found == "True" ||
                          r.polyA_motif_found == "true");

  auto it = chains.find(r.isoform);
  ev.chain_available = (it != chains.end());
  ev.catalog_evaluable = ev.chain_available && catalog != nullptr;

  // Multi-caller consensus axis: how many independent callers recovered this isoform's
  // intron chain (from a `combine` matrix). Joined by the chain fingerprint when the
  // matrix carries one (two callers may reuse a native id for different chains); the
  // native id is the fallback for older matrices and monoexonic isoforms.
  if (caller_support != nullptr) {
    const int* n = nullptr;
    if (caller_support->has_fingerprints() && ev.chain_available && !it->second.monoexonic()) {
      auto cs = caller_support->by_chain.find(fingerprint_intron_chain(it->second));
      if (cs != caller_support->by_chain.end()) n = &cs->second;
    } else {
      auto cs = caller_support->by_native_id.find(r.isoform);
      if (cs != caller_support->by_native_id.end()) n = &cs->second;
    }
    if (n != nullptr) {
      ev.consensus_evaluable = true;
      ev.n_callers = *n;
    }
  }

  // Novel-junction axes (short-read corroboration + BAM mapping features) require
  // the caller chain AND a catalog to know which junctions are novel.
  if (ev.catalog_evaluable) {
    const IntronChain& chain = it->second;
    BamFeatureParams bp;
    bp.min_mapq = cfg.bam_min_mapq;
    bp.softclip_min_bp = cfg.bam_softclip_min_bp;
    bp.junction_window_bp = cfg.bam_junction_window_bp;

    int n_novel = 0, n_supported = 0, n_pan = 0, n_variant_created = 0;
    int n_known = 0, n_known_supported = 0, n_bam_gated = 0;
    bool any_bam_eval = false;
    bool any_variant_eval = false;
    bool any_pan_eval = false;
    // Haplotype/path consistency: a transcript is one molecule from one haplotype, so a
    // rescue needs ONE haplotype (or graph path) that explains every novel junction, not
    // a different one per junction. Intersect the per-junction sets.
    std::vector<bool> hap_common(haplotype != nullptr ? haplotype->n_haplotypes() : 0, true);
    std::vector<bool> hap_at;
    std::vector<std::string> path_common;
    bool path_ids_complete = true;  // every graph-supported junction listed its path ids
    bool first_path = true;
    auto sr_supported = [&](const auto& intron) {
      const SjRecord* hit = sj->find_exact(chain.chrom, chain.strand, intron);
      return hit != nullptr && hit->n_uniq >= cfg.sj_min_uniq_reads &&
             (!cfg.sj_require_canonical_motif || hit->canonical());
    };
    for (const auto& intron : chain.introns) {
      if (catalog->has_intron(chain.chrom, chain.strand, intron)) {
        // Known junction: its short-read coverage tells whether the locus is covered at all
        // (the rule engine treats "0 reads at the novel junction" as absence only then).
        ++n_known;
        if (sj != nullptr && sr_supported(intron)) ++n_known_supported;
        continue;
      }
      ++n_novel;
      if (haplotype != nullptr) {
        const VariantVerdict vv = haplotype->classify(chain.chrom, intron, chain.strand, &hap_at);
        if (vv != VariantVerdict::kNotEvaluable) any_variant_eval = true;
        // Only kCreated (non-canonical on reference, canonical on a haplotype) drives a
        // reference-bias rescue. kDisrupted (canonical on reference, broken on every
        // haplotype) is computed by HaplotypeProvider but intentionally NOT consumed
        // yet -- its adjudication policy is undefined (see docs/decision_engine.md,
        // "Future work: Variant-disrupted handling").
        if (vv == VariantVerdict::kCreated) {
          ++n_variant_created;
          for (std::size_t h = 0; h < hap_common.size(); ++h) hap_common[h] = hap_common[h] && hap_at[h];
        }
      }
      if (pangenome != nullptr) {
        any_pan_eval = true;
        if (pangenome->supports(chain.chrom, chain.strand, intron, cfg.pangenome_min_haplotypes)) {
          ++n_pan;
          const auto& ids = pangenome->find_exact(chain.chrom, chain.strand, intron)->haplotype_ids;
          if (ids.empty()) {
            path_ids_complete = false;
          } else if (first_path) {
            path_common = ids;
            first_path = false;
          } else {
            std::vector<std::string> keep;
            for (const auto& id : path_common)
              if (std::find(ids.begin(), ids.end(), id) != ids.end()) keep.push_back(id);
            path_common.swap(keep);
          }
        }
      }
      if (sj != nullptr && sr_supported(intron)) ++n_supported;
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
          // Fractions are only meaningful with enough reads: with 1-2 reads a single indel
          // is a fraction of 1.0. Junctions below bam_min_spanning do not feed the gate.
          if (f.n_spanning >= cfg.bam_min_spanning) {
            ++n_bam_gated;
            if (f.frac_low_mapq() > ev.bam_max_frac_low_mapq) ev.bam_max_frac_low_mapq = f.frac_low_mapq();
            if (f.frac_supplementary() > ev.bam_max_frac_supplementary) ev.bam_max_frac_supplementary = f.frac_supplementary();
            // All four fractions feed RuleEngine::evaluate's mapping_flag (the softclip gate
            // is disabled by default via its threshold) and are emitted to attribution.jsonl.
            if (f.frac_softclip() > ev.bam_max_frac_softclip) ev.bam_max_frac_softclip = f.frac_softclip();
            if (f.frac_indel_near() > ev.bam_max_frac_indel_near) ev.bam_max_frac_indel_near = f.frac_indel_near();
          }
        }
      }
    }
    ev.n_novel_junctions = n_novel;  // recorded even without SR for reference
    ev.n_known_junctions = n_known;
    if (sj != nullptr) {
      ev.sj_evaluable = true;
      ev.n_novel_jx_sr_supported = n_supported;
      ev.n_known_jx_sr_supported = n_known_supported;
    }
    ev.bam_evaluable = (bam != nullptr) && any_bam_eval;
    ev.bam_n_junctions_gated = n_bam_gated;
    ev.variant_evaluable = (haplotype != nullptr) && any_variant_eval;
    ev.n_novel_jx_variant_created = n_variant_created;
    ev.haplotype_consistent = n_variant_created > 0 &&
                              std::find(hap_common.begin(), hap_common.end(), true) != hap_common.end();
    // Rescue requires ALL of the isoform's novel junctions to be explained by the
    // reference-bias evidence -- a single explained junction does not make a
    // multi-novel-junction isoform "false novel" (the others may be genuinely new) --
    // and explained on ONE haplotype (the per-junction count alone would accept junction
    // 1 on haplotype A plus junction 2 on haplotype B).
    ev.variant_rescue = (n_novel > 0 && n_variant_created == n_novel && ev.haplotype_consistent);
    // Keep the circular-risk flag consistent with evaluability (it is only read by the
    // rescue branch, which already requires *_evaluable, so this is a clarity guard).
    ev.variant_circular = ev.variant_evaluable && haplotype_circular;
    ev.pangenome_evaluable = (pangenome != nullptr) && any_pan_eval;
    ev.n_novel_jx_pangenome = n_pan;
    // Path identity is optional in the junction file; when every supporting row carries
    // it, the rescue also needs one path shared by all novel junctions.
    ev.pangenome_path_evaluable = n_pan > 0 && path_ids_complete;
    ev.pangenome_path_consistent = ev.pangenome_path_evaluable && !path_common.empty();
    ev.pangenome_rescue = (n_novel > 0 && n_pan == n_novel &&
                           (!ev.pangenome_path_evaluable || ev.pangenome_path_consistent));
    ev.pangenome_circular = ev.pangenome_evaluable && pangenome_circular;
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
    // Pass-through SQANTI3 descriptors for bio_flags (not used in the verdict).
    res.bio_dist_to_CAGE_peak = r.dist_to_CAGE_peak;
    res.bio_dist_to_polyA_site = r.dist_to_polyA_site;
    res.bio_predicted_NMD = r.predicted_NMD;
    res.bio_polyA_motif_found = r.polyA_motif_found;
    res.evidence = build_evidence(r, *in.chains, in.catalog, in.sj, in.bam, in.haplotype,
                                  in.haplotype_circular, in.pangenome, in.pangenome_circular,
                                  in.caller_support, engine.config(), bam_cache);
    res.verdict = engine.evaluate(res.evidence);
    out.push_back(std::move(res));
  }
  return out;
}

}  // namespace panisoguard
