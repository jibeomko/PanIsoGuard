#include "catch2/catch.hpp"

#include "panisoguard/rules.hpp"

using namespace panisoguard;

namespace {
// A novel (NNC) isoform with `n` novel junctions, `k` of them SR-supported, and
// optional artifact-mechanism flags.
EvidenceVector novel_ev(int n, int k, bool noncanonical = false, bool rts = false,
                        double percA = -1.0) {
  EvidenceVector e;
  e.structural_category = "novel_not_in_catalog";
  e.is_novel = true;
  e.chain_available = true;
  e.sj_evaluable = true;
  e.n_novel_junctions = n;
  e.n_novel_jx_sr_supported = k;
  e.canon_evaluable = true;
  e.noncanonical = noncanonical;
  e.rts_evaluable = true;
  e.rts_stage = rts;
  if (percA >= 0) {
    e.percA_evaluable = true;
    e.perc_A_downstream_TTS = percA;
  }
  return e;
}
}  // namespace

TEST_CASE("rule projection over the 2-axis grid", "[rules]") {
  const RuleEngine eng;  // defaults

  CHECK(eng.evaluate(novel_ev(2, 2)).confidence == ConfidenceClass::kHighConfNovel);
  CHECK(eng.evaluate(novel_ev(2, 2, /*noncanonical=*/true)).confidence == ConfidenceClass::kMediumConfNovel);
  CHECK(eng.evaluate(novel_ev(2, 1)).confidence == ConfidenceClass::kMediumConfNovel);
  CHECK(eng.evaluate(novel_ev(2, 1, /*noncanonical=*/true)).confidence == ConfidenceClass::kLowConfPartial);
  CHECK(eng.evaluate(novel_ev(2, 0)).confidence == ConfidenceClass::kLowConfPartial);
  CHECK(eng.evaluate(novel_ev(2, 0, /*noncanonical=*/true)).confidence == ConfidenceClass::kArtifact);
  // UNSUPPORTED x rt_switch and x degradation also -> ARTIFACT (explicit grid coverage)
  CHECK(eng.evaluate(novel_ev(2, 0, false, /*rts=*/true)).confidence == ConfidenceClass::kArtifact);
  CHECK(eng.evaluate(novel_ev(2, 0, false, false, /*percA=*/70.0)).confidence == ConfidenceClass::kArtifact);

  // perc_A intra-priming threshold matches SQANTI3 default (flag >= 60)
  CHECK(eng.evaluate(novel_ev(2, 0, false, false, /*percA=*/60.0)).confidence == ConfidenceClass::kArtifact);
  CHECK(eng.evaluate(novel_ev(2, 0, false, false, /*percA=*/59.5)).confidence == ConfidenceClass::kLowConfPartial);
}

TEST_CASE("rules: unevaluable novelty support is AMBIGUOUS, not a guess", "[rules]") {
  EvidenceVector e = novel_ev(0, 0);
  e.sj_evaluable = false;  // no SR / catalog
  e.n_novel_junctions = 0;
  const Verdict v = RuleEngine().evaluate(e);
  CHECK(v.confidence == ConfidenceClass::kAmbiguous);
  CHECK_FALSE(v.rule_trace.empty());
  CHECK_FALSE(v.circularity_flag);  // Tier-0 never sets circularity
}

TEST_CASE("rules: BAM mapping artifact axis", "[rules]") {
  const RuleEngine eng;

  // SR-supported novel but reads map poorly -> demoted to MEDIUM (conflicting).
  EvidenceVector sup = novel_ev(1, 1);
  sup.bam_evaluable = true;
  sup.bam_n_spanning_total = 10;
  sup.bam_max_frac_low_mapq = 0.8;  // > default 0.5
  Verdict vs = eng.evaluate(sup);
  CHECK(vs.primary_mechanism == Mechanism::kMapping);
  CHECK(vs.confidence == ConfidenceClass::kMediumConfNovel);

  // Unsupported novel + mapping artifact -> ARTIFACT.
  EvidenceVector uns = novel_ev(1, 0);
  uns.bam_evaluable = true;
  uns.bam_n_spanning_total = 10;
  uns.bam_max_frac_supplementary = 0.9;  // > default 0.5
  CHECK(eng.evaluate(uns).confidence == ConfidenceClass::kArtifact);

  // SR axis absent (UNKNOWN support) but a strong mapping artifact is decisive.
  EvidenceVector unk = novel_ev(0, 0);
  unk.sj_evaluable = false;
  unk.n_novel_junctions = 0;
  unk.bam_evaluable = true;
  unk.bam_n_spanning_total = 5;
  unk.bam_max_frac_low_mapq = 0.9;
  CHECK(eng.evaluate(unk).confidence == ConfidenceClass::kArtifact);
}

TEST_CASE("rules: indel-near and soft-clip are mapping-artifact triggers", "[rules]") {
  const RuleEngine eng;  // defaults: indel_near / softclip thresholds both 0.5

  // An indel adjacent to the junction on most spanning reads (an alignment-ambiguous
  // indel mis-rendered as an intron) flags the mapping mechanism even when MAPQ and
  // supplementary fractions are clean. Validated on chr22: catches 41/124 false novel
  // junctions at 100% precision (benchmark/bam_axis).
  EvidenceVector indel = novel_ev(1, 0);
  indel.bam_evaluable = true;
  indel.bam_n_spanning_total = 20;
  indel.bam_max_frac_indel_near = 0.8;  // > default 0.5
  Verdict vi = eng.evaluate(indel);
  CHECK(vi.primary_mechanism == Mechanism::kMapping);
  CHECK(vi.confidence == ConfidenceClass::kArtifact);

  // Terminal soft-clipping on most spanning reads (reads that could not align through)
  // is the same mapping mechanism.
  EvidenceVector clip = novel_ev(1, 0);
  clip.bam_evaluable = true;
  clip.bam_n_spanning_total = 20;
  clip.bam_max_frac_softclip = 0.7;  // > default 0.5
  CHECK(eng.evaluate(clip).primary_mechanism == Mechanism::kMapping);

  // Below threshold (the genuine-junction regime, chr22 genuine max ~0.19): no mapping
  // flag, so a clean SR-supported novel stays HIGH (zero collateral damage).
  EvidenceVector ok = novel_ev(1, 1);
  ok.bam_evaluable = true;
  ok.bam_n_spanning_total = 20;
  ok.bam_max_frac_indel_near = 0.19;
  ok.bam_max_frac_softclip = 0.04;
  Verdict vo = eng.evaluate(ok);
  CHECK(vo.primary_mechanism == Mechanism::kNone);
  CHECK(vo.confidence == ConfidenceClass::kHighConfNovel);
}

TEST_CASE("rules: multi-caller consensus corroborates an UNKNOWN-support novel chain", "[rules]") {
  const RuleEngine eng;  // defaults: consensus_min_callers = 2

  // A novel isoform with NO short-read/catalog axis (UNKNOWN support) and no artifact
  // mechanism is held AMBIGUOUS by default...
  EvidenceVector base = novel_ev(0, 0);
  base.sj_evaluable = false;
  base.n_novel_junctions = 0;
  CHECK(eng.evaluate(base).confidence == ConfidenceClass::kAmbiguous);

  // ...but >= 2 independent callers recovering the chain promotes it to MEDIUM (never HIGH:
  // caller agreement is methodological, not experimental, corroboration).
  EvidenceVector consensus = base;
  consensus.consensus_evaluable = true;
  consensus.n_callers = 3;
  Verdict vc = eng.evaluate(consensus);
  CHECK(vc.confidence == ConfidenceClass::kMediumConfNovel);
  CHECK_FALSE(vc.rule_trace.empty());

  // A single caller does not reach the gate -> still AMBIGUOUS.
  EvidenceVector single = base;
  single.consensus_evaluable = true;
  single.n_callers = 1;
  CHECK(eng.evaluate(single).confidence == ConfidenceClass::kAmbiguous);

  // Consensus + an artifact mechanism (non-canonical) -> LOW_CONF_PARTIAL, not MEDIUM.
  EvidenceVector mech = consensus;
  mech.noncanonical = true;
  CHECK(eng.evaluate(mech).confidence == ConfidenceClass::kLowConfPartial);

  // A strong mapping artifact still dominates consensus -> ARTIFACT.
  EvidenceVector mapping = consensus;
  mapping.bam_evaluable = true;
  mapping.bam_n_spanning_total = 10;
  mapping.bam_max_frac_low_mapq = 0.9;
  CHECK(eng.evaluate(mapping).confidence == ConfidenceClass::kArtifact);

  // Ablating the consensus axis reverts the promotion (AMBIGUOUS again).
  CHECK(eng.with_axis_disabled("consensus").evaluate(consensus).confidence ==
        ConfidenceClass::kAmbiguous);
}

TEST_CASE("rules: known/partial categories pass through", "[rules]") {
  EvidenceVector fsm;
  fsm.structural_category = "full-splice_match";
  CHECK(RuleEngine().evaluate(fsm).confidence == ConfidenceClass::kHighConfKnown);

  EvidenceVector ism;
  ism.structural_category = "incomplete-splice_match";
  CHECK(RuleEngine().evaluate(ism).confidence == ConfidenceClass::kLowConfPartial);
}
