#include "panisoguard/rules.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

#include "tomlplusplus/toml.hpp"

namespace panisoguard {

RuleEngine RuleEngine::from_toml(const std::string& path) {
  RuleEngine engine;
  RuleConfig& c = engine.cfg_;
  const toml::table tbl = toml::parse_file(path);  // throws toml::parse_error on failure

  // Every key is optional (a missing key keeps the default). A key that is present must have
  // the right type and a value inside its range; otherwise the run stops here instead of
  // running with a threshold that was silently clamped or converted.
  auto fail = [&](const std::string& key, const std::string& why) {
    throw std::runtime_error(path + ": " + key + " " + why);
  };
  auto get_int = [&](const toml::table* t, const char* section, const char* key, int64_t lo, int64_t hi, int& dst) {
    const toml::node* n = t == nullptr ? nullptr : t->get(key);
    if (n == nullptr) return;
    const std::string name = std::string(section) + "." + key;
    if (!n->is_integer()) fail(name, "must be an integer");
    const int64_t v = n->as_integer()->get();
    if (v < lo || v > hi) fail(name, "must be between " + std::to_string(lo) + " and " + std::to_string(hi));
    dst = static_cast<int>(v);
  };
  auto get_double = [&](const toml::table* t, const char* section, const char* key, double lo, double hi, double& dst) {
    const toml::node* n = t == nullptr ? nullptr : t->get(key);
    if (n == nullptr) return;
    const std::string name = std::string(section) + "." + key;
    if (!n->is_number()) fail(name, "must be a number");
    const double v = n->is_integer() ? static_cast<double>(n->as_integer()->get()) : n->as_floating_point()->get();
    if (!std::isfinite(v) || v < lo || v > hi) fail(name, "must be between " + std::to_string(lo) + " and " + std::to_string(hi));
    dst = v;
  };
  auto get_bool = [&](const toml::table* t, const char* section, const char* key, bool& dst) {
    const toml::node* n = t == nullptr ? nullptr : t->get(key);
    if (n == nullptr) return;
    if (!n->is_boolean()) fail(std::string(section) + "." + key, "must be true or false");
    dst = n->as_boolean()->get();
  };
  auto get_string = [&](const toml::table* t, const char* section, const char* key, std::string& dst) {
    const toml::node* n = t == nullptr ? nullptr : t->get(key);
    if (n == nullptr) return;
    if (!n->is_string()) fail(std::string(section) + "." + key, "must be a string");
    dst = n->as_string()->get();
  };

  const toml::table* meta = tbl["meta"].as_table();
  get_string(meta, "meta", "sqanti3_version_target", c.sqanti3_version_target);
  get_string(meta, "meta", "ruleset_version", c.ruleset_version);

  const toml::table* ns = tbl["axis_novelty_support"].as_table();
  get_int(ns, "axis_novelty_support", "sj_min_uniq_reads", 1, 1000000, c.sj_min_uniq_reads);
  get_bool(ns, "axis_novelty_support", "sj_require_canonical_motif", c.sj_require_canonical_motif);
  get_int(ns, "axis_novelty_support", "consensus_min_callers", 1, 1000, c.consensus_min_callers);

  const toml::table* art = tbl["axis_artifact"].as_table();
  const toml::table* deg = art == nullptr ? nullptr : (*art)["degradation"].as_table();
  get_double(deg, "axis_artifact.degradation", "max_perc_A_downstream_TTS", 0.0, 100.0, c.perc_A_degradation_threshold);
  const toml::table* map = art == nullptr ? nullptr : (*art)["mapping"].as_table();
  get_int(map, "axis_artifact.mapping", "min_mapq", 0, 255, c.bam_min_mapq);
  get_int(map, "axis_artifact.mapping", "min_spanning", 1, 1000000, c.bam_min_spanning);
  get_int(map, "axis_artifact.mapping", "junction_window_bp", 0, 10000, c.bam_junction_window_bp);
  // A fraction gate above 1.0 can never fire: that is how a gate is switched off.
  get_double(map, "axis_artifact.mapping", "max_low_mapq_frac", 0.0, 2.0, c.bam_max_low_mapq_frac);
  get_double(map, "axis_artifact.mapping", "max_supplementary_frac", 0.0, 2.0, c.bam_max_supplementary_frac);
  get_double(map, "axis_artifact.mapping", "max_indel_near_frac", 0.0, 2.0, c.bam_max_indel_near_frac);

  get_int(tbl["axis_chain"].as_table(), "axis_chain", "min_full_reads", 1, 1000000, c.chain_min_full_reads);
  get_int(tbl["axis_pangenome"].as_table(), "axis_pangenome", "min_haplotypes", 1, 1000000, c.pangenome_min_haplotypes);
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

std::string RuleEngine::thresholds_summary() const {
  const RuleConfig& c = cfg_;
  auto num = [](double x) { std::ostringstream o; o << x; return o.str(); };  // 0.5, not 0.500000
  return "sj_min_uniq_reads=" + std::to_string(c.sj_min_uniq_reads) +
         " sj_require_canonical_motif=" + (c.sj_require_canonical_motif ? "true" : "false") +
         " consensus_min_callers=" + std::to_string(c.consensus_min_callers) +
         " chain_min_full_reads=" + std::to_string(c.chain_min_full_reads) +
         " max_perc_A_downstream_TTS=" + num(c.perc_A_degradation_threshold) +
         " platform=" + c.platform +
         " min_mapq=" + std::to_string(c.bam_min_mapq) +
         " min_spanning=" + std::to_string(c.bam_min_spanning) +
         " junction_window_bp=" + std::to_string(c.bam_junction_window_bp) +
         " max_low_mapq_frac=" + num(c.bam_max_low_mapq_frac) +
         " max_supplementary_frac=" + num(c.bam_max_supplementary_frac) +
         " max_indel_near_frac=" + num(c.bam_max_indel_near_frac) +
         " pangenome_min_haplotypes=" + std::to_string(c.pangenome_min_haplotypes);
}

// ---------------------------------------------------------------------------------------
// evaluate(): five steps, each a function below. The order is the rulebook of
// docs/decision_engine.md: pass-through categories, the reference-bias rescue, the support
// axis, the mechanism axis, and the projection of the two onto a class. Every step appends
// to rule_trace, so a verdict can be read back from the record alone.
// ---------------------------------------------------------------------------------------
namespace {

void trace(Verdict& v, std::string s) { v.rule_trace.push_back(std::move(s)); }

void finish(Verdict& v, ConfidenceClass c, const char* verdict, std::string reason) {
  v.confidence = c;
  v.verdict = verdict;
  v.verdict_reason = std::move(reason);
}

// The signs, as booleans: SQANTI3's QC priors and the BAM mapping sign. Intra-priming needs
// BOTH an A-rich stretch after the 3' end (perc_A) AND no polyA motif upstream of it (many
// perc_A-flagged ends carry a motif and are genuine terminations; benchmark/refbias_3way).
// Without a motif column perc_A alone decides. The mapping fractions are taken only over novel
// junctions with >= bam_min_spanning reads (with one or two reads a single indel is a fraction
// of 1.0). A motif that is canonical on the reference only, non-canonical on every haplotype
// of the sample, is a non-canonical sign, under the same firewall as the rescue.
struct Signs {
  bool percA_high = false;
  bool motif_known = false;
  bool intrapriming = false;
  bool rts = false;
  bool mapping_sign = false;
  bool disrupted_sign = false;
};

Signs read_signs(const RuleConfig& cfg, const EvidenceVector& ev) {
  Signs s;
  s.percA_high = ev.percA_evaluable &&
                 ev.perc_A_downstream_TTS >= cfg.perc_A_degradation_threshold;
  s.motif_known = ev.polyA_motif_evaluable;
  s.intrapriming = s.percA_high && !(s.motif_known && ev.polyA_motif_found);
  s.rts = ev.rts_evaluable && ev.rts_stage;
  const bool mapping_gated = ev.bam_evaluable && ev.bam_n_junctions_gated > 0;
  s.mapping_sign = mapping_gated && (ev.bam_max_frac_low_mapq > cfg.bam_max_low_mapq_frac ||
                                     ev.bam_max_frac_supplementary > cfg.bam_max_supplementary_frac ||
                                     ev.bam_max_frac_indel_near > cfg.bam_max_indel_near_frac);
  s.disrupted_sign = ev.variant_evaluable &&
                     ev.n_novel_jx_variant_disrupted > 0 && !ev.variant_circular;
  return s;
}

// Signs that were seen but could not be used; said once, after the pass-through (which does
// not read them).
void trace_unused_signs(const RuleConfig& cfg, const EvidenceVector& ev, Verdict& v) {
  if (ev.bam_evaluable && ev.bam_n_junctions_gated == 0) {
    trace(v, "mapping axis not evaluated: fewer than " + std::to_string(cfg.bam_min_spanning) +
                 " spanning reads at every novel junction (" + std::to_string(ev.bam_n_spanning_total) + " in total)");
  }
  if (ev.variant_evaluable && ev.n_novel_jx_variant_disrupted > 0 && ev.variant_circular) {
    trace(v, std::to_string(ev.n_novel_jx_variant_disrupted) +
                 " novel junction(s) variant-disrupted on the haplotypes, but provenance is circular-risk -> sign ignored");
  }
}

// --- 0. Pass-through categories: no novelty to adjudicate. True when the verdict is final.
bool pass_through(const RuleConfig& cfg, const EvidenceVector& ev, const Signs& s, Verdict& v) {
  if (ev.structural_category == "full-splice_match") {
    v.novelty_support = NoveltySupport::kSupported;
    trace(v, "category=full-splice_match -> KNOWN");
    finish(v, ConfidenceClass::kKnown, "known", "full_splice_match");
    return true;
  }
  if (ev.structural_category == "incomplete-splice_match") {
    // A truncated known transcript: the splice structure is known, only the 3' end is
    // unreliable, so the intra-priming sign goes into the reason, not the class. An ISM whose
    // whole intron chain >= chain_min_full_reads long reads carry, with no intra-priming sign,
    // is a KNOWN_FRAGMENT: a known structure seen end to end, kept apart from the novel classes
    // because nothing novel was tested.
    const bool chain_ok = ev.chain_evaluable && ev.chain_n_full_reads >= cfg.chain_min_full_reads;
    if (chain_ok && !s.intrapriming) {
      trace(v, "category=incomplete-splice_match, " + std::to_string(ev.chain_n_full_reads) +
                   " whole-chain read(s), no intra-priming sign -> KNOWN_FRAGMENT");
      finish(v, ConfidenceClass::kKnownFragment, "known_fragment",
             std::string("ism_chain_support") + (s.rts ? ";rts_sign" : "") +
                 (ev.chain_independent ? ";independent_library" : ""));
      return true;
    }
    std::string reason = "ism";
    if (s.intrapriming) reason = "ism_intrapriming";
    else if (s.rts) reason = "ism_rts";
    trace(v, "category=incomplete-splice_match -> UNCONFIRMED (" + reason + ")");
    finish(v, ConfidenceClass::kUnconfirmed, "unconfirmed", reason);
    return true;
  }
  if (!ev.is_novel) {
    trace(v, "category=" + ev.structural_category + " (not a targeted novel class) -> UNRESOLVED");
    finish(v, ConfidenceClass::kUnresolved, "not_adjudicated", "category_" + ev.structural_category);
    return true;
  }
  return false;
}

// --- 1. Reference-bias rescue: the pangenome and the sample's haplotypes. True when final.
// Rescue needs ALL of the isoform's novel junctions explained (an isoform is "reference bias"
// only if its whole novelty is); the per-junction count is in the evidence. Both axes pass the
// circularity firewall: provenance that is not independent of the sample's RNA holds the call
// instead of promoting it, but only when neither axis offers an independent rescue. A rescue
// is not blocked by a mapping sign (that blocked a correct rescue in refbias_phasing); the
// sign is carried in the reason so the call can be reviewed.
bool reference_bias(const EvidenceVector& ev, const Signs& s, Verdict& v) {
  const std::string mapping_note = s.mapping_sign ? ";mapping_sign" : "";
  const std::string mapping_trace = s.mapping_sign ? "; mapping sign present" : "";
  const bool all_pan = ev.n_novel_junctions > 0 && ev.n_novel_jx_pangenome == ev.n_novel_junctions;
  const bool all_created = ev.n_novel_junctions > 0 && ev.n_novel_jx_variant_created == ev.n_novel_junctions;
  const bool pan_rescue = ev.pangenome_evaluable && ev.pangenome_rescue;
  const bool var_rescue = ev.variant_evaluable && ev.variant_rescue;
  const std::string pan_frac = std::to_string(ev.n_novel_jx_pangenome) + "/" + std::to_string(ev.n_novel_junctions);
  // Independent evidence on either axis promotes; a circular-risk rescue holds the call only
  // when no independent rescue applies.
  if (pan_rescue && !ev.pangenome_circular) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    trace(v, "all " + pan_frac + " novel junctions realizable on a pangenome haplotype path (independent "
                 "population data) -> REFERENCE_BIAS (population_known)" + mapping_trace);
    finish(v, ConfidenceClass::kReferenceBias, "reference_bias", "pangenome_path" + mapping_note);
    return true;
  }
  if (var_rescue && !ev.variant_circular) {
    v.primary_mechanism = Mechanism::kVariant;
    trace(v, "variant creates canonical splice motif on haplotype, non-canonical on reference -> "
             "REFERENCE_BIAS (reference bias)" + mapping_trace);
    finish(v, ConfidenceClass::kReferenceBias, "reference_bias", "haplotype_variant" + mapping_note);
    return true;
  }
  if (pan_rescue) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    v.circularity_flag = true;
    trace(v, "all " + pan_frac + " novel junctions realizable on a pangenome path, but the junction-set "
                 "provenance is sample-derived/unknown (circular-risk) -> held UNRESOLVED (not promoted)");
    finish(v, ConfidenceClass::kUnresolved, "unconfirmed", "firewall_held");
    return true;
  }
  if (var_rescue) {
    v.primary_mechanism = Mechanism::kVariant;
    v.circularity_flag = true;
    trace(v, "variant creates canonical motif on haplotype but provenance is circular-risk "
             "(RNA-derived/unknown) -> held UNRESOLVED (not promoted)");
    finish(v, ConfidenceClass::kUnresolved, "unconfirmed", "firewall_held");
    return true;
  }
  // Every novel junction is explained individually, but by different haplotypes / graph
  // paths: no single haplotype carries the chain, so the isoform is not reference bias of
  // one allele. The evidence conflicts -> held, not rescued and not called an artifact.
  if (ev.pangenome_evaluable && all_pan && ev.pangenome_path_evaluable &&
      !ev.pangenome_path_consistent) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    trace(v, "all novel junctions are on pangenome paths, but no single path carries them all -> held UNRESOLVED");
    finish(v, ConfidenceClass::kUnresolved, "unconfirmed", "path_inconsistent");
    return true;
  }
  if (ev.variant_evaluable && all_created && !ev.haplotype_consistent) {
    v.primary_mechanism = Mechanism::kVariant;
    trace(v, "every novel junction is variant-created, but on different haplotypes (no haplotype "
             "carries the whole chain) -> held UNRESOLVED");
    finish(v, ConfidenceClass::kUnresolved, "unconfirmed", "haplotype_inconsistent");
    return true;
  }
  return false;
}

// --- 2. Axis A: novelty support from short-read corroboration. `unknown_reason` names why
// the axis is UNKNOWN when it is.
NoveltySupport support_axis(const EvidenceVector& ev, Verdict& v, std::string& unknown_reason) {
  if (!ev.sj_evaluable || !ev.chain_available || ev.n_novel_junctions == 0) {
    if (!ev.chain_available) {
      unknown_reason = "caller_chain_absent";
      trace(v, "novelty-support UNKNOWN reason=caller_chain_absent");
    } else if (!ev.sj_evaluable) {
      unknown_reason = "no_short_reads";
      trace(v, "novelty-support UNKNOWN reason=short_read/catalog_axis_absent");
    } else {
      unknown_reason = "no_novel_junction";
      trace(v, "novelty-support UNKNOWN reason=no_novel_junctions (e.g. NIC combinatorial novelty)");
    }
    return NoveltySupport::kUnknown;
  }
  const std::string n = std::to_string(ev.n_novel_junctions);
  if (ev.n_novel_jx_sr_supported == ev.n_novel_junctions) {
    trace(v, "sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" + n + " -> SUPPORTED");
    return NoveltySupport::kSupported;
  }
  if (ev.n_novel_jx_sr_supported == 0) {
    // Absence of reads at the novel junctions is evidence against them only if the locus is
    // covered at all: when none of the isoform's KNOWN junctions has short-read support
    // either, the sample is too shallow there and the support stays UNKNOWN.
    if (ev.n_known_junctions > 0 && ev.n_known_jx_sr_supported == 0) {
      unknown_reason = "no_coverage";
      trace(v, "sj_support=0/" + n + " but 0/" + std::to_string(ev.n_known_junctions) +
                   " known junctions supported either -> no coverage, not absence: UNKNOWN reason=no_coverage");
      return NoveltySupport::kUnknown;
    }
    trace(v, "sj_support=0/" + n + " -> UNSUPPORTED" +
                 (ev.n_known_junctions > 0
                      ? " (" + std::to_string(ev.n_known_jx_sr_supported) + "/" + std::to_string(ev.n_known_junctions) +
                            " known junctions covered)"
                      : " (no known junctions to gauge coverage)"));
    return NoveltySupport::kUnsupported;
  }
  trace(v, "sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" + n + " -> PARTIAL");
  return NoveltySupport::kPartial;
}

// --- 3. Axis B: the dominant artifact mechanism. Priority: mapping > noncanonical > rt_switch
// > degradation. Alignment reliability is the most upstream concern: if the spanning reads may
// be misplaced, the downstream motif/QC signs are moot. Only the first match is reported.
Mechanism mechanism_axis(const RuleConfig& cfg, const EvidenceVector& ev, const Signs& s, Verdict& v) {
  if (s.mapping_sign) {
    trace(v, "bam mapping artifact (low_mapq_frac=" + std::to_string(ev.bam_max_frac_low_mapq) +
                 ", supplementary_frac=" + std::to_string(ev.bam_max_frac_supplementary) +
                 ", indel_near_frac=" + std::to_string(ev.bam_max_frac_indel_near) +
                 ", over " + std::to_string(ev.bam_n_junctions_gated) + " junction(s) with >= " +
                 std::to_string(cfg.bam_min_spanning) + " reads) -> mechanism=mapping_or_repeat");
    return Mechanism::kMapping;
  }
  if (ev.canon_evaluable && ev.noncanonical) {
    trace(v, "all_canonical=non_canonical -> mechanism=noncanonical");
    return Mechanism::kNoncanonical;
  }
  if (s.disrupted_sign) {
    // SQANTI3 judged the motif on the reference; on every haplotype of this sample it is
    // non-canonical, so for these reads the junction IS non-canonical.
    trace(v, std::to_string(ev.n_novel_jx_variant_disrupted) +
                 " novel junction(s) canonical on the reference but non-canonical on every haplotype "
                 "(variant-disrupted, independent provenance) -> mechanism=noncanonical");
    return Mechanism::kNoncanonical;
  }
  if (s.rts) {
    trace(v, "RTS_stage=TRUE -> mechanism=rt_switch");
    return Mechanism::kRtSwitch;
  }
  if (s.intrapriming) {
    trace(v, "perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) + " >= " +
                 std::to_string(cfg.perc_A_degradation_threshold) +
                 (s.motif_known ? " and no polyA motif upstream" : " (polyA motif column absent: perc_A alone)") +
                 " -> mechanism=degradation");
    return Mechanism::kDegradation;
  }
  if (s.percA_high) {
    trace(v, "perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) +
                 " >= threshold but a polyA motif is found upstream -> not intra-priming");
  }
  trace(v, "no artifact mechanism flagged -> mechanism=none");
  return Mechanism::kNone;
}

// --- 4. Projection: (support x mechanism) -> class, verdict, reason.
struct Call {
  ConfidenceClass cls;
  const char* verdict;
  std::string reason;
};

// Short-read support is not evaluable. A mapping sign holds the call; otherwise caller
// consensus or whole-chain long-read support, both methodological (not experimental)
// corroboration, lift it to PROBABLE, never CONFIRMED; with neither it is held UNRESOLVED.
// The RT-switch sign does not demote methodological support (benchmark/rts_sensitivity);
// the intra-priming sign still holds, as for ISM above.
Call project_unknown(const RuleConfig& cfg, const EvidenceVector& ev, Mechanism mech,
                     const std::string& unknown_reason, Verdict& v) {
  const std::string mc = to_string(mech);
  const bool weak_sign = mech == Mechanism::kRtSwitch;
  const std::string sign_note = weak_sign ? ";" + mc + "_sign" : "";
  const bool consensus = ev.consensus_evaluable && ev.n_callers >= cfg.consensus_min_callers;
  const bool chain_reads = ev.chain_evaluable && ev.chain_n_full_reads >= cfg.chain_min_full_reads;
  const bool examined = ev.chain_available && ev.catalog_evaluable;
  if (mech == Mechanism::kMapping) {
    // Without short reads a mapping sign alone is not decisive (a third of such calls were
    // genuine on SQANTI-SIM), so the call is held UNCONFIRMED.
    return {ConfidenceClass::kUnconfirmed, "unconfirmed", "mapping_without_short_reads"};
  }
  if (consensus && !examined) {
    // Caller agreement can only corroborate a chain we could examine: without the caller's
    // chain or a reference catalog nothing about the isoform was verified, so the consensus
    // count is reported but does not promote.
    trace(v, "caller-consensus n_callers=" + std::to_string(ev.n_callers) + " ignored: " +
                 (ev.chain_available ? "no reference catalog" : "caller chain absent") + " -> UNRESOLVED");
    return {ConfidenceClass::kUnresolved, "unconfirmed", unknown_reason};
  }
  if (consensus) {
    trace(v, "caller-consensus n_callers=" + std::to_string(ev.n_callers) + " >= " +
                 std::to_string(cfg.consensus_min_callers) +
                 " corroborates the novel chain (no short-read axis) -> consensus-supported");
    if (mech == Mechanism::kNone || weak_sign) return {ConfidenceClass::kNovelProbable, "confirmed", "caller_consensus" + sign_note};
    return {ConfidenceClass::kUnconfirmed, "unconfirmed", "consensus_with_" + mc};
  }
  if (chain_reads && examined) {
    // Whole-chain long-read support: reads that carry every intron of the chain in one
    // alignment. They are the caller's own reads, not independent evidence, and they separate
    // caller-level false chains, not alignment-level ones (benchmark/chain_support), so
    // PROBABLE, never CONFIRMED.
    trace(v, "whole-chain long-read support: " + std::to_string(ev.chain_n_full_reads) +
                 " read(s) carry the entire intron chain (>= " + std::to_string(cfg.chain_min_full_reads) +
                 (ev.chain_independent ? "; independent library" : "; same reads as the caller") + ") -> chain-supported");
    const std::string src = ev.chain_independent ? ";independent_library" : "";
    if (mech == Mechanism::kNone || weak_sign) return {ConfidenceClass::kNovelProbable, "confirmed", "chain_support" + sign_note + src};
    return {ConfidenceClass::kUnconfirmed, "unconfirmed", "chain_support_with_" + mc + src};
  }
  return {ConfidenceClass::kUnresolved, "unconfirmed", unknown_reason};
}

// Short-read support was evaluated: the grid of docs/decision_engine.md.
Call project_grid(NoveltySupport sup, Mechanism mech) {
  const std::string mc = to_string(mech);
  const bool clean = mech == Mechanism::kNone;
  switch (sup) {
    case NoveltySupport::kSupported:
      return clean ? Call{ConfidenceClass::kNovelConfirmed, "confirmed", "all_junctions_supported"}
                   : Call{ConfidenceClass::kNovelProbable, "confirmed", "supported_with_" + mc};
    case NoveltySupport::kPartial:
      // An uncorroborated novel junction leaves the chain unconfirmed, so PARTIAL never
      // promotes to a confident novel call (it did, and those were the extra false
      // positives in benchmark/sqanti3_filter_h2h).
      return {ConfidenceClass::kUnconfirmed, "unconfirmed", clean ? "partial_support" : "partial_support_with_" + mc};
    default:  // UNSUPPORTED
      return clean ? Call{ConfidenceClass::kUnconfirmed, "unconfirmed", "unsupported_no_sign"}
                   : Call{ConfidenceClass::kArtifact, "artifact", "unsupported_" + mc};
  }
}

void project(const RuleConfig& cfg, const EvidenceVector& ev, const Signs& s, NoveltySupport sup, Mechanism mech,
             const std::string& unknown_reason, Verdict& v) {
  Call c = sup == NoveltySupport::kUnknown ? project_unknown(cfg, ev, mech, unknown_reason, v)
                                           : project_grid(sup, mech);
  trace(v, std::string("project(") + to_string(sup) + "," + to_string(mech) + ") -> " + to_string(c.cls));
  // Reaching the grid means no reference-bias rescue fired, so this verdict never rests on
  // RNA-derived variants.
  v.circularity_flag = false;
  if (s.disrupted_sign && mech == Mechanism::kNoncanonical && !(ev.canon_evaluable && ev.noncanonical))
    c.reason += ";haplotype_disrupted";  // the sign came from the sample's haplotypes, not from SQANTI3
  finish(v, c.cls, c.verdict, std::move(c.reason));
}

}  // namespace

Verdict RuleEngine::evaluate(const EvidenceVector& ev) const {
  Verdict v;
  const Signs s = read_signs(cfg_, ev);
  if (pass_through(cfg_, ev, s, v)) return v;
  trace_unused_signs(cfg_, ev, v);
  if (reference_bias(ev, s, v)) return v;
  std::string unknown_reason;
  v.novelty_support = support_axis(ev, v, unknown_reason);
  v.primary_mechanism = mechanism_axis(cfg_, ev, s, v);
  project(cfg_, ev, s, v.novelty_support, v.primary_mechanism, unknown_reason, v);
  return v;
}

}  // namespace panisoguard
