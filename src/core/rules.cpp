#include "panisoguard/rules.hpp"

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
    c.consensus_min_callers = static_cast<int>((*ns)["consensus_min_callers"].value_or<int64_t>(c.consensus_min_callers));
  }
  if (auto art = tbl["axis_artifact"].as_table()) {
    if (auto deg = (*art)["degradation"].as_table()) {
      c.perc_A_degradation_threshold =
          (*deg)["max_perc_A_downstream_TTS"].value_or(c.perc_A_degradation_threshold);
    }
  }
  return engine;
}

Verdict RuleEngine::evaluate(const EvidenceVector& ev) const {
  Verdict v;
  auto trace = [&](std::string s) { v.rule_trace.push_back(std::move(s)); };

  // --- Pass-through for non-targeted categories --------------------------------
  if (ev.structural_category == "full-splice_match") {
    v.confidence = ConfidenceClass::kHighConfKnown;
    v.novelty_support = NoveltySupport::kSupported;
    trace("category=full-splice_match -> HIGH_CONF_KNOWN");
    return v;
  }
  if (ev.structural_category == "incomplete-splice_match") {
    v.confidence = ConfidenceClass::kLowConfPartial;
    trace("category=incomplete-splice_match -> LOW_CONF_PARTIAL");
    return v;
  }
  if (!ev.is_novel) {
    v.confidence = ConfidenceClass::kAmbiguous;
    trace("category=" + ev.structural_category + " (not a targeted novel class) -> AMBIGUOUS");
    return v;
  }

  // --- Axis A: novelty support from short-read corroboration -------------------
  NoveltySupport sup;
  if (!ev.sj_evaluable || !ev.chain_available || ev.n_novel_junctions == 0) {
    sup = NoveltySupport::kUnknown;
    trace("novelty-support not evaluable (sj/chain absent or no novel junctions) -> UNKNOWN");
  } else if (ev.n_novel_jx_sr_supported == ev.n_novel_junctions) {
    sup = NoveltySupport::kSupported;
    trace("sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" +
          std::to_string(ev.n_novel_junctions) + " -> SUPPORTED");
  } else if (ev.n_novel_jx_sr_supported == 0) {
    sup = NoveltySupport::kUnsupported;
    trace("sj_support=0/" + std::to_string(ev.n_novel_junctions) + " -> UNSUPPORTED");
  } else {
    sup = NoveltySupport::kPartial;
    trace("sj_support=" + std::to_string(ev.n_novel_jx_sr_supported) + "/" +
          std::to_string(ev.n_novel_junctions) + " -> PARTIAL");
  }
  v.novelty_support = sup;

  // --- Axis B: dominant artifact mechanism (priority order) --------------------
  Mechanism mech = Mechanism::kNone;
  if (ev.canon_evaluable && ev.noncanonical) {
    mech = Mechanism::kNoncanonical;
    trace("all_canonical=non_canonical -> mechanism=noncanonical");
  } else if (ev.rts_evaluable && ev.rts_stage) {
    mech = Mechanism::kRtSwitch;
    trace("RTS_stage=TRUE -> mechanism=rt_switch");
  } else if (ev.percA_evaluable && ev.perc_A_downstream_TTS > cfg_.perc_A_degradation_threshold) {
    mech = Mechanism::kDegradation;
    trace("perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) + " > " +
          std::to_string(cfg_.perc_A_degradation_threshold) + " -> mechanism=degradation");
  } else {
    trace("no artifact mechanism flagged -> mechanism=none");
  }
  v.primary_mechanism = mech;

  // --- Projection: (support x mechanism) -> confidence class -------------------
  ConfidenceClass cls;
  if (sup == NoveltySupport::kUnknown) {
    cls = ConfidenceClass::kAmbiguous;
  } else if (sup == NoveltySupport::kSupported) {
    cls = (mech == Mechanism::kNone) ? ConfidenceClass::kHighConfNovel
                                     : ConfidenceClass::kMediumConfNovel;
  } else if (sup == NoveltySupport::kPartial) {
    cls = (mech == Mechanism::kNone) ? ConfidenceClass::kMediumConfNovel
                                     : ConfidenceClass::kLowConfPartial;
  } else {  // UNSUPPORTED
    cls = (mech == Mechanism::kNone) ? ConfidenceClass::kLowConfPartial
                                     : ConfidenceClass::kArtifact;
  }
  v.confidence = cls;
  trace(std::string("project(") + to_string(sup) + "," + to_string(mech) + ") -> " +
        to_string(cls));

  // Tier-0 carries no variant/pangenome axis, so no rescue classes are reachable
  // and no verdict rests on RNA-derived variants.
  v.circularity_flag = false;
  return v;
}

}  // namespace panisoguard
