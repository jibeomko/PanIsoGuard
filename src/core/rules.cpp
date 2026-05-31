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
    if (auto map = (*art)["mapping"].as_table()) {
      c.bam_min_mapq = static_cast<int>((*map)["min_mapq"].value_or<int64_t>(c.bam_min_mapq));
      c.bam_softclip_min_bp = static_cast<int>((*map)["softclip_min_bp"].value_or<int64_t>(c.bam_softclip_min_bp));
      c.bam_junction_window_bp = static_cast<int>((*map)["junction_window_bp"].value_or<int64_t>(c.bam_junction_window_bp));
      c.bam_max_low_mapq_frac = (*map)["max_low_mapq_frac"].value_or(c.bam_max_low_mapq_frac);
      c.bam_max_supplementary_frac = (*map)["max_supplementary_frac"].value_or(c.bam_max_supplementary_frac);
    }
  }
  if (auto pg = tbl["axis_pangenome"].as_table()) {
    c.pangenome_min_haplotypes =
        static_cast<int>((*pg)["min_haplotypes"].value_or<int64_t>(c.pangenome_min_haplotypes));
    if (c.pangenome_min_haplotypes < 1) c.pangenome_min_haplotypes = 1;  // a 0/negative gate would rescue everything
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

  // --- Pangenome (reference-bias) rescue: highest priority ---------------------
  // When ALL of an isoform's novel-vs-linear-reference junctions are realizable on a
  // pangenome graph haplotype path, the apparent novelty is consistent with reference
  // bias rather than new splicing. This is independent evidence ONLY when the junction
  // set was extracted from population assemblies; PanIsoGuard cannot verify the
  // supplied file's provenance, so it applies the same circularity firewall as the
  // variant axis (sample-derived/unknown provenance is held, not promoted). It is
  // checked before the variant axis so an independent graph rescues even when the
  // sample's own haplotype provenance is circular-risk.
  if (cfg_.use_pangenome && ev.pangenome_evaluable && ev.pangenome_rescue) {
    v.primary_mechanism = Mechanism::kPopulationKnown;
    const std::string frac = std::to_string(ev.n_novel_jx_pangenome) + "/" +
                             std::to_string(ev.n_novel_junctions);
    if (ev.pangenome_circular) {
      v.confidence = ConfidenceClass::kAmbiguous;
      v.circularity_flag = true;
      trace("all " + frac + " novel junctions realizable on a pangenome path, but the "
            "junction-set provenance is sample-derived/unknown (circular-risk) -> held "
            "AMBIGUOUS (not promoted)");
    } else {
      v.confidence = ConfidenceClass::kPanRefRescuedFalseNovel;
      trace("all " + frac + " novel junctions realizable on a pangenome haplotype path "
            "(independent population data) -> PAN_REF_RESCUED_FALSE_NOVEL "
            "(reference bias, population_known)");
    }
    return v;
  }

  // --- Variant (reference-bias) rescue: highest priority -----------------------
  // A novel junction that is non-canonical on the reference but canonical on the
  // sample's haplotype is explained by reference bias, not new splicing.
  if (cfg_.use_variant && ev.variant_evaluable && ev.variant_rescue) {
    v.primary_mechanism = Mechanism::kVariant;
    if (ev.variant_circular) {
      // RNA-derived/unknown provenance: not independent evidence -> never promoted.
      v.confidence = ConfidenceClass::kAmbiguous;
      v.circularity_flag = true;
      trace("variant creates canonical motif on haplotype but provenance is circular-risk "
            "(RNA-derived/unknown) -> held AMBIGUOUS (not promoted)");
    } else {
      v.confidence = ConfidenceClass::kPanRefRescuedFalseNovel;
      trace("variant creates canonical splice motif on haplotype, non-canonical on reference "
            "-> PAN_REF_RESCUED_FALSE_NOVEL (reference bias)");
    }
    return v;
  }

  // --- Axis A: novelty support from short-read corroboration -------------------
  const bool sj_ok = ev.sj_evaluable && cfg_.use_short_read;
  NoveltySupport sup;
  if (!sj_ok || !ev.chain_available || ev.n_novel_junctions == 0) {
    sup = NoveltySupport::kUnknown;
    // Distinguish the unknowability mode for downstream filtering.
    if (!ev.chain_available) {
      trace("novelty-support UNKNOWN reason=caller_chain_absent");
    } else if (!ev.sj_evaluable) {
      trace("novelty-support UNKNOWN reason=short_read/catalog_axis_absent");
    } else {
      trace("novelty-support UNKNOWN reason=no_novel_junctions (e.g. NIC combinatorial novelty)");
    }
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
  // Priority: mapping > noncanonical > rt_switch > degradation. Rationale:
  // alignment reliability is the most upstream concern -- if the spanning reads
  // do not map confidently (low MAPQ / supplementary / multimapping), the junction
  // may be a mapping artifact and downstream motif/QC signals are moot. Only the
  // first-matching mechanism is reported as primary (see rule_trace for all flags).
  const bool mapping_flag = cfg_.use_mapping && ev.bam_evaluable && ev.bam_n_spanning_total > 0 &&
                            (ev.bam_max_frac_low_mapq > cfg_.bam_max_low_mapq_frac ||
                             ev.bam_max_frac_supplementary > cfg_.bam_max_supplementary_frac);
  Mechanism mech = Mechanism::kNone;
  if (mapping_flag) {
    mech = Mechanism::kMapping;
    trace("bam mapping artifact (low_mapq_frac=" + std::to_string(ev.bam_max_frac_low_mapq) +
          ", supplementary_frac=" + std::to_string(ev.bam_max_frac_supplementary) +
          ") -> mechanism=mapping_or_repeat");
  } else if (cfg_.use_noncanonical && ev.canon_evaluable && ev.noncanonical) {
    mech = Mechanism::kNoncanonical;
    trace("all_canonical=non_canonical -> mechanism=noncanonical");
  } else if (cfg_.use_rts && ev.rts_evaluable && ev.rts_stage) {
    mech = Mechanism::kRtSwitch;
    trace("RTS_stage=TRUE -> mechanism=rt_switch");
  } else if (cfg_.use_degradation && ev.percA_evaluable &&
             ev.perc_A_downstream_TTS >= cfg_.perc_A_degradation_threshold) {
    mech = Mechanism::kDegradation;
    trace("perc_A_downstream_TTS=" + std::to_string(ev.perc_A_downstream_TTS) + " >= " +
          std::to_string(cfg_.perc_A_degradation_threshold) + " -> mechanism=degradation");
  } else {
    trace("no artifact mechanism flagged -> mechanism=none");
  }
  v.primary_mechanism = mech;

  // --- Projection: (support x mechanism) -> confidence class -------------------
  ConfidenceClass cls;
  if (sup == NoveltySupport::kUnknown) {
    // A strong mapping artifact is decisive even when short-read support is not
    // evaluable; otherwise the call is held as AMBIGUOUS.
    cls = (mech == Mechanism::kMapping) ? ConfidenceClass::kArtifact : ConfidenceClass::kAmbiguous;
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

  // Reaching the grid means no reference-bias rescue fired (pangenome/variant
  // branches return early), so this verdict never rests on RNA-derived variants.
  v.circularity_flag = false;
  return v;
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
  return e;
}

}  // namespace panisoguard
