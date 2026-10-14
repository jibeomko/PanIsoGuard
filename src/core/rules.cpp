#include "panisoguard/rules.hpp"

#include <stdexcept>

#include <string>

#include "tomlplusplus/toml.hpp"

namespace panisoguard {

RuleEngine RuleEngine::from_toml(const std::string& path) {
  RuleEngine engine;
  RuleConfig& c = engine.cfg_;

  const toml::table tbl = toml::parse_file(path);  // throws toml::parse_error on failure

  if (auto meta = tbl["meta"].as_table()) {
    c.sqanti3_version_target = (*meta)["sqanti3_version_target"].value_or(c.sqanti3_version_target);
    c.ruleset_version = (*meta)["ruleset_version"].value_or(c.ruleset_version);
  }
  if (auto ns = tbl["axis_novelty_support"].as_table()) {
    c.sj_min_uniq_reads = static_cast<int>((*ns)["sj_min_uniq_reads"].value_or<int64_t>(c.sj_min_uniq_reads));
    c.sj_require_canonical_motif = (*ns)["sj_require_canonical_motif"].value_or(c.sj_require_canonical_motif);
    // Multi-caller agreement gate: a novel chain recovered by >= this many callers is
    // consensus-supported (consumed by evaluate() when --caller-support populates
    // consensus_evaluable/n_callers; inert otherwise).
    c.consensus_min_callers = static_cast<int>((*ns)["consensus_min_callers"].value_or<int64_t>(c.consensus_min_callers));
  }
  if (auto art = tbl["axis_artifact"].as_table()) {
    if (auto deg = (*art)["degradation"].as_table()) {
      c.perc_A_degradation_threshold =
          (*deg)["max_perc_A_downstream_TTS"].value_or(c.perc_A_degradation_threshold);
    }
    if (auto map = (*art)["mapping"].as_table()) {
      c.bam_min_mapq = static_cast<int>((*map)["min_mapq"].value_or<int64_t>(c.bam_min_mapq));
      c.bam_min_spanning = static_cast<int>((*map)["min_spanning"].value_or<int64_t>(c.bam_min_spanning));
      c.bam_softclip_min_bp = static_cast<int>((*map)["softclip_min_bp"].value_or<int64_t>(c.bam_softclip_min_bp));
      c.bam_junction_window_bp = static_cast<int>((*map)["junction_window_bp"].value_or<int64_t>(c.bam_junction_window_bp));
      c.bam_max_low_mapq_frac = (*map)["max_low_mapq_frac"].value_or(c.bam_max_low_mapq_frac);
      c.bam_max_supplementary_frac = (*map)["max_supplementary_frac"].value_or(c.bam_max_supplementary_frac);
      c.bam_max_indel_near_frac = (*map)["max_indel_near_frac"].value_or(c.bam_max_indel_near_frac);
      c.bam_max_softclip_frac = (*map)["max_softclip_frac"].value_or(c.bam_max_softclip_frac);
    }
  }
  if (auto pg = tbl["axis_pangenome"].as_table()) {
    c.pangenome_min_haplotypes =
        static_cast<int>((*pg)["min_haplotypes"].value_or<int64_t>(c.pangenome_min_haplotypes));
    if (c.pangenome_min_haplotypes < 1) c.pangenome_min_haplotypes = 1;  // a 0/negative gate would rescue everything
  }
  return engine;
}

void RuleEngine::apply_platform_preset(const std::string& platform) {
  cfg_.platform = platform;
  if (platform == "ont") {
    // The indel-near gate is calibrated on PacBio HiFi; on ONT it fires on ~90 % of junctions
    // (benchmark/gm12878_realdata) because of the chemistry's indel error, so it is off. The
    // low-MAPQ and supplementary gates do not depend on the error profile and stay on.
    cfg_.bam_max_indel_near_frac = 1.01;
  }
}

namespace {

// Short code for the mechanism inside a verdict_reason (the output name, e.g. mapping_or_repeat).
std::string mech_code(Mechanism m) { return to_string(m); }

}  // namespace

Verdict RuleEngine::evaluate(const EvidenceVector& ev) const {
  Verdict v;
  auto trace = [&](std::string s) { v.rule_trace.push_back(std::move(s)); };
  auto finish = [&](ConfidenceClass c, const char* verdict, std::string reason) {
    v.confidence = c;
    v.verdict = verdict;
    v.verdict_reason = std::move(reason);
    return v;
  };

  // SQANTI3 QC signs consumed as priors. Intra-priming needs BOTH an A-rich stretch after the
  // 3' end (perc_A) AND no polyA motif upstream of it: on real data, 20-30 % of perc_A-flagged
  // ends carry a polyA signal and are genuine terminations (benchmark/refbias_3way, HG03516).
  // When SQANTI3 was run without a motif list the column is absent and perc_A alone decides.
  const bool percA_high = cfg_.use_degradation && ev.percA_evaluable &&
                          ev.perc_A_downstream_TTS >= cfg_.perc_A_degradation_threshold;
  const bool motif_known = ev.polyA_motif_evaluable;
  const bool intrapriming = percA_high && !(motif_known && ev.polyA_motif_found);
  const bool rts = cfg_.use_rts && ev.rts_evaluable && ev.rts_stage;

  // --- 0. Pass-through categories: no novelty to adjudicate -----------------------
  if (ev.structural_category == "full-splice_match") {
    v.novelty_support = NoveltySupport::kSupported;
    trace("category=full-splice_match -> KNOWN");
    return finish(ConfidenceClass::kKnown, "known", "full_splice_match");
  }
  if (ev.structural_category == "incomplete-splice_match") {
    // A truncated known transcript. 61 % of real ISM carry the intra-priming sign (HG03516),
    // so the sign is reported in the reason, but the class does not change: the splice
    // structure is known; only the 3' end is unreliable, and that is not a splicing artifact.
    std::string reason = "ism";
    if (intrapriming) reason = "ism_intrapriming";
    else if (rts) reason = "ism_rts";
    trace("category=incomplete-splice_match -> UNCONFIRMED (" + reason + ")");
    return finish(ConfidenceClass::kUnconfirmed, "unconfirmed", reason);
  }
  if (!ev.is_novel) {
    trace("category=" + ev.structural_category + " (not a targeted novel class) -> UNRESOLVED");
    return finish(ConfidenceClass::kUnresolved, "not_adjudicated", "category_" + ev.structural_category);
  }

  // --- BAM mapping sign (used by the rescue trace and the mechanism axis) ------------
  // Fractions are taken only over novel junctions with >= bam_min_spanning reads: with one
  // or two reads a single indel makes the fraction 1.0 and would decide the verdict alone
  // (45 % of HG03516's mapping-artifact calls rested on < 3 reads before this gate).
  const bool mapping_gated = cfg_.use_mapping && ev.bam_evaluable && ev.bam_n_junctions_gated > 0;
  const bool mapping_sign = mapping_gated &&
                            (ev.bam_max_frac_low_mapq > cfg_.bam_max_low_mapq_frac ||
                             ev.bam_max_frac_supplementary > cfg_.bam_max_supplementary_frac ||
                             ev.bam_max_frac_indel_near > cfg_.bam_max_indel_near_frac ||
                             ev.bam_max_frac_softclip > cfg_.bam_max_softclip_frac);
  if (cfg_.use_mapping && ev.bam_evaluable && ev.bam_n_junctions_gated == 0) {
    trace("mapping axis not evaluated: fewer than " + std::to_string(cfg_.bam_min_spanning) +
          " spanning reads at every novel junction (" + std::to_string(ev.bam_n_spanning_total) + " in total)");
  }
  const std::string mapping_note = mapping_sign ? ";mapping_sign" : "";

  // --- 1. Reference-bias rescue: pangenome first, then the sample's haplotypes ------
  // Rescue needs ALL of the isoform's novel junctions explained (an isoform is "reference
  // bias" only if its whole novelty is); the per-junction count is in the evidence. Both
  // axes pass the circularity firewall: provenance that is not independent of the sample's
  // RNA holds the call instead of promoting it. The pangenome axis is checked first so an
  // independent graph rescues even when the haplotype provenance is circular-risk. A rescue
  // is not blocked by a mapping sign (that blocked a correct rescue in refbias_phasing), but
  // the sign is carried in the reason so the call can be reviewed.
  const bool all_pan = ev.n_novel_junctions > 0 && ev.n_novel_jx_pangenome == ev.n_novel_junctions;
  const bool all_created = ev.n_novel_junctions > 0 && ev.n_novel_jx_variant_created == ev.n_novel_junctions;
  if (cfg_.use_pangenome && ev.pangenome_evaluable && ev.pangenome_rescue) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    const std::string frac = std::to_string(ev.n_novel_jx_pangenome) + "/" + std::to_string(ev.n_novel_junctions);
    if (ev.pangenome_circular) {
      v.circularity_flag = true;
      trace("all " + frac + " novel junctions realizable on a pangenome path, but the junction-set "
            "provenance is sample-derived/unknown (circular-risk) -> held UNRESOLVED (not promoted)");
      return finish(ConfidenceClass::kUnresolved, "unconfirmed", "firewall_held");
    }
    trace("all " + frac + " novel junctions realizable on a pangenome haplotype path (independent "
          "population data) -> REFERENCE_BIAS (population_known)" + (mapping_sign ? "; mapping sign present" : ""));
    return finish(ConfidenceClass::kReferenceBias, "reference_bias", "pangenome_path" + mapping_note);
  }
  if (cfg_.use_variant && ev.variant_evaluable && ev.variant_rescue) {
    v.primary_mechanism = Mechanism::kVariant;
    if (ev.variant_circular) {
      v.circularity_flag = true;
      trace("variant creates canonical motif on haplotype but provenance is circular-risk "
            "(RNA-derived/unknown) -> held UNRESOLVED (not promoted)");
      return finish(ConfidenceClass::kUnresolved, "unconfirmed", "firewall_held");
    }
    trace("variant creates canonical splice motif on haplotype, non-canonical on reference -> "
          "REFERENCE_BIAS (reference bias)" + std::string(mapping_sign ? "; mapping sign present" : ""));
    return finish(ConfidenceClass::kReferenceBias, "reference_bias", "haplotype_variant" + mapping_note);
  }
  // Every novel junction is explained individually, but by different haplotypes / graph
  // paths: no single haplotype carries the chain, so the isoform is not reference bias of
  // one allele. The evidence conflicts -> held, not rescued and not called an artifact.
  if (cfg_.use_pangenome && ev.pangenome_evaluable && all_pan && ev.pangenome_path_evaluable &&
      !ev.pangenome_path_consistent) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    trace("all novel junctions are on pangenome paths, but no single path carries them all -> held UNRESOLVED");
    return finish(ConfidenceClass::kUnresolved, "unconfirmed", "path_inconsistent");
  }
  if (cfg_.use_variant && ev.variant_evaluable && all_created && !ev.haplotype_consistent) {
    v.primary_mechanism = Mechanism::kVariant;
    trace("every novel junction is variant-created, but on different haplotypes (no haplotype "
          "carries the whole chain) -> held UNRESOLVED");
    return finish(ConfidenceClass::kUnresolved, "unconfirmed", "haplotype_inconsistent");
  }

  // --- 2. Axis A: novelty support from short-read corroboration ----------------------
  const bool sj_ok = ev.sj_evaluable && cfg_.use_short_read;
  NoveltySupport sup;
  std::string unknown_reason;
  if (!sj_ok || !ev.chain_available || ev.n_novel_junctions == 0) {
    sup = NoveltySupport::kUnknown;
    if (!ev.chain_available) {
      unknown_reason = "caller_chain_absent";
      trace("novelty-support UNKNOWN reason=caller_chain_absent");
    } else if (!ev.sj_evaluable) {
      unknown_reason = "no_short_reads";
      trace("novelty-support UNKNOWN reason=short_read/catalog_axis_absent");
    } else {
      unknown_reason = "no_novel_junction";
      trace("novelty-support UNKNOWN reason=no_novel_junctions (e.g. NIC combinatorial novelty)");
    }
  } else if (ev.n_novel_jx_sr_supported == ev.n_novel_junctions) {
    sup = NoveltySupport::kSupported;
    trace("sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" +
          std::to_string(ev.n_novel_junctions) + " -> SUPPORTED");
  } else if (ev.n_novel_jx_sr_supported == 0) {
    // Absence of reads at the novel junctions is evidence against them only if the locus is
    // covered at all: when none of the isoform's KNOWN junctions has short-read support
    // either, the sample is too shallow there and the support stays UNKNOWN.
    if (ev.n_known_junctions > 0 && ev.n_known_jx_sr_supported == 0) {
      sup = NoveltySupport::kUnknown;
      unknown_reason = "no_coverage";
      trace("sj_support=0/" + std::to_string(ev.n_novel_junctions) + " but 0/" +
            std::to_string(ev.n_known_junctions) + " known junctions supported either -> no coverage, "
            "not absence: UNKNOWN reason=no_coverage");
    } else {
      sup = NoveltySupport::kUnsupported;
      trace("sj_support=0/" + std::to_string(ev.n_novel_junctions) + " -> UNSUPPORTED" +
            (ev.n_known_junctions > 0
                 ? " (" + std::to_string(ev.n_known_jx_sr_supported) + "/" + std::to_string(ev.n_known_junctions) +
                       " known junctions covered)"
                 : " (no known junctions to gauge coverage)"));
    }
  } else {
    sup = NoveltySupport::kPartial;
    trace("sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" +
          std::to_string(ev.n_novel_junctions) + " -> PARTIAL");
  }
  v.novelty_support = sup;

  // --- 3. Axis B: dominant artifact mechanism (priority order) -------------------------
  // Priority: mapping > noncanonical > rt_switch > degradation. Alignment reliability is
  // the most upstream concern: if the spanning reads may be misplaced, the downstream
  // motif/QC signs are moot. Only the first-matching mechanism is reported as primary.
  Mechanism mech = Mechanism::kNone;
  if (mapping_sign) {
    mech = Mechanism::kMapping;
    trace("bam mapping artifact (low_mapq_frac=" + std::to_string(ev.bam_max_frac_low_mapq) +
          ", supplementary_frac=" + std::to_string(ev.bam_max_frac_supplementary) +
          ", indel_near_frac=" + std::to_string(ev.bam_max_frac_indel_near) +
          ", softclip_frac=" + std::to_string(ev.bam_max_frac_softclip) +
          ", over " + std::to_string(ev.bam_n_junctions_gated) + " junction(s) with >= " +
          std::to_string(cfg_.bam_min_spanning) + " reads) -> mechanism=mapping_or_repeat");
  } else if (cfg_.use_noncanonical && ev.canon_evaluable && ev.noncanonical) {
    mech = Mechanism::kNoncanonical;
    trace("all_canonical=non_canonical -> mechanism=noncanonical");
  } else if (rts) {
    mech = Mechanism::kRtSwitch;
    trace("RTS_stage=TRUE -> mechanism=rt_switch");
  } else if (intrapriming) {
    mech = Mechanism::kDegradation;
    trace("perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) + " >= " +
          std::to_string(cfg_.perc_A_degradation_threshold) +
          (motif_known ? " and no polyA motif upstream" : " (polyA motif column absent: perc_A alone)") +
          " -> mechanism=degradation");
  } else {
    if (percA_high) {
      trace("perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) +
            " >= threshold but a polyA motif is found upstream -> not intra-priming");
    }
    trace("no artifact mechanism flagged -> mechanism=none");
  }
  v.primary_mechanism = mech;
  const std::string mc = mech_code(mech);

  // --- 4. Projection: (support x mechanism) -> class, verdict, reason -------------------
  ConfidenceClass cls;
  const char* verdict;
  std::string reason;
  if (sup == NoveltySupport::kUnknown) {
    // A strong mapping artifact is decisive even when short-read support is not
    // evaluable. Otherwise, multi-caller consensus can corroborate the novel chain:
    // >= consensus_min_callers independent callers recovering the same intron chain is
    // methodological (not experimental) support, so it promotes only to PROBABLE, never
    // CONFIRMED. With neither signal the call is held UNRESOLVED.
    if (mech == Mechanism::kMapping) {
      cls = ConfidenceClass::kArtifact; verdict = "artifact"; reason = "mapping_without_short_reads";
    } else if (cfg_.use_consensus && ev.consensus_evaluable && ev.n_callers >= cfg_.consensus_min_callers &&
               !(ev.chain_available && ev.catalog_evaluable)) {
      // Caller agreement can only corroborate a chain we could examine: without the
      // caller's chain or a reference catalog nothing about the isoform was verified, so
      // the consensus count is reported but does not promote.
      trace("caller-consensus n_callers=" + std::to_string(ev.n_callers) + " ignored: " +
            (ev.chain_available ? "no reference catalog" : "caller chain absent") + " -> UNRESOLVED");
      cls = ConfidenceClass::kUnresolved; verdict = "unconfirmed"; reason = unknown_reason;
    } else if (cfg_.use_consensus && ev.consensus_evaluable && ev.n_callers >= cfg_.consensus_min_callers) {
      trace("caller-consensus n_callers=" + std::to_string(ev.n_callers) + " >= " +
            std::to_string(cfg_.consensus_min_callers) +
            " corroborates the novel chain (no short-read axis) -> consensus-supported");
      if (mech == Mechanism::kNone) { cls = ConfidenceClass::kNovelProbable; verdict = "confirmed"; reason = "caller_consensus"; }
      else { cls = ConfidenceClass::kUnconfirmed; verdict = "unconfirmed"; reason = "consensus_with_" + mc; }
    } else {
      cls = ConfidenceClass::kUnresolved; verdict = "unconfirmed"; reason = unknown_reason;
    }
  } else if (sup == NoveltySupport::kSupported) {
    if (mech == Mechanism::kNone) { cls = ConfidenceClass::kNovelConfirmed; verdict = "confirmed"; reason = "all_junctions_supported"; }
    else { cls = ConfidenceClass::kNovelProbable; verdict = "confirmed"; reason = "supported_with_" + mc; }
  } else if (sup == NoveltySupport::kPartial) {
    // An uncorroborated novel junction leaves the chain unconfirmed, so PARTIAL never
    // promotes to a confident novel call (it did, and those were the extra false
    // positives in benchmark/sqanti3_filter_h2h).
    cls = ConfidenceClass::kUnconfirmed; verdict = "unconfirmed";
    reason = (mech == Mechanism::kNone) ? "partial_support" : "partial_support_with_" + mc;
  } else {  // UNSUPPORTED
    if (mech == Mechanism::kNone) { cls = ConfidenceClass::kUnconfirmed; verdict = "unconfirmed"; reason = "unsupported_no_sign"; }
    else { cls = ConfidenceClass::kArtifact; verdict = "artifact"; reason = "unsupported_" + mc; }
  }
  trace(std::string("project(") + to_string(sup) + "," + to_string(mech) + ") -> " + to_string(cls));

  // Reaching the grid means no reference-bias rescue fired (pangenome/variant
  // branches return early), so this verdict never rests on RNA-derived variants.
  v.circularity_flag = false;
  return finish(cls, verdict, reason);
}

RuleEngine RuleEngine::with_axis_disabled(const std::string& axis) const {
  RuleEngine e = *this;
  if (axis == "short_read")        e.cfg_.use_short_read = false;
  else if (axis == "mapping")      e.cfg_.use_mapping = false;
  else if (axis == "noncanonical") e.cfg_.use_noncanonical = false;
  else if (axis == "rt_switch")    e.cfg_.use_rts = false;
  else if (axis == "degradation")  e.cfg_.use_degradation = false;
  else if (axis == "variant")      e.cfg_.use_variant = false;
  else if (axis == "pangenome")    e.cfg_.use_pangenome = false;
  else if (axis == "consensus")    e.cfg_.use_consensus = false;
  else throw std::invalid_argument("unknown ablation axis '" + axis + "' (expected one of: short_read, mapping, "
                                   "noncanonical, rt_switch, degradation, variant, pangenome, consensus)");
  return e;
}

}  // namespace panisoguard
