#include "panisoguard/adjudicator.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "panisoguard/fingerprint.hpp"

namespace panisoguard {

ReadPhase read_phase(const std::vector<JunctionRead>& reads, const std::vector<PhaseSite>& sites,
                     const std::vector<bool>& on, int margin) {
  ReadPhase p;
  p.n_sites = static_cast<int>(sites.size());
  std::vector<int> match(on.size());
  for (const JunctionRead& r : reads) {
    std::fill(match.begin(), match.end(), 0);
    for (const PhaseSite& s : sites) {   // in position order
      if (s.pos < r.start) continue;
      if (s.pos >= r.end) break;
      const char c = r.base_at(s.pos);
      if (c == 0) continue;
      for (std::size_t h = 0; h < on.size() && h < s.bases.size(); ++h) if (s.bases[h] == c) ++match[h];
    }
    int motif = 0, other = 0;   // the best agreement inside each group of haplotypes
    for (std::size_t h = 0; h < on.size(); ++h) {
      int& best = on[h] ? motif : other;
      best = std::max(best, match[h]);
    }
    if (motif - other >= margin) ++p.n_motif;
    else if (other - motif >= margin) ++p.n_other;
  }
  return p;
}

namespace {

constexpr std::size_t kJunctionReadCap = 200;   // junction reads looked at per junction

// What the variant axis found at one novel junction; the same for every isoform that uses it.
struct VariantAt {
  JunctionMotifs jm;
  std::vector<bool> on;            // per haplotype: it explains the junction (canonical motif, or it lacks the intron)
  std::string read_phase = "NA";
  ReadPhase reads;
};

VariantAt variant_at(const HaplotypeProvider& haplotype, const BamReader* bam, const std::string& chrom,
                     const Junction& intron, Strand strand, const RuleConfig& cfg) {
  VariantAt va;
  std::vector<JunctionRead> reads;
  bool loaded = false;
  auto load = [&]() {
    if (!loaded && bam != nullptr) reads = bam->junction_reads(chrom, intron, strand, kJunctionReadCap);
    loaded = true;
  };
  // An insertion exactly at an exon boundary: the reads carry it (with an insertion of that
  // length within a few bases on the exon side) or they do not. Read fractions need
  // bam_min_spanning reads, as on the mapping axis; between 20% and 80% the reads do not settle it.
  const BoundaryOracle oracle = [&](bool left, int len) {
    load();
    const int n = static_cast<int>(reads.size());
    if (len <= 0 || n < cfg.bam_min_spanning) return BoundaryReads::kUnknown;
    const int64_t b = left ? intron.start : intron.end;
    int k = 0;
    for (const JunctionRead& r : reads) {
      for (const auto& ins : r.insertions) {
        if (ins.second != len) continue;
        if (left ? (ins.first >= b - 3 - len && ins.first <= b) : (ins.first >= b && ins.first <= b + 3 + len)) { ++k; break; }
      }
    }
    return 5 * k <= n ? BoundaryReads::kWithout : 5 * k >= 4 * n ? BoundaryReads::kWith : BoundaryReads::kUnknown;
  };
  va.jm = haplotype.motifs(chrom, intron, strand, bam != nullptr ? oracle : BoundaryOracle());
  const JunctionMotifs& jm = va.jm;
  va.on.assign(jm.haps.size(), false);
  bool all_known = true, some = false, every = !jm.haps.empty();
  for (std::size_t h = 0; h < jm.haps.size(); ++h) {
    va.on[h] = jm.deleted[h] || (jm.haps[h].evaluable && jm.haps[h].canonical);
    if (!jm.deleted[h] && !jm.haps[h].evaluable) all_known = false;
    if (va.on[h]) some = true; else every = false;
  }
  // A junction explained on some haplotypes only is used by reads of those haplotypes: the
  // sample's other heterozygous sites say which haplotype each junction read comes from.
  const bool explained = jm.verdict == VariantVerdict::kCreated || jm.verdict == VariantVerdict::kDeleted;
  if (explained && bam != nullptr && some && !every && all_known) {
    load();
    va.read_phase = "uninformative";
    if (!reads.empty()) {
      int64_t lo = reads.front().start, hi = reads.front().end;
      for (const JunctionRead& r : reads) { lo = std::min(lo, r.start); hi = std::max(hi, r.end); }
      va.reads = read_phase(reads, haplotype.phase_sites(chrom, lo, hi, jm), va.on, cfg.variant_phase_margin);
      if (va.reads.n_motif + va.reads.n_other >= cfg.variant_phase_min_reads)
        va.read_phase = va.reads.n_other == 0 ? "consistent" : va.reads.n_motif == 0 ? "contradicted" : "mixed";
    }
  }
  return va;
}

// Assemble the evidence for one SQANTI record. `bam_cache` and `variant_cache` memoize the
// per-junction BAM features and variant findings across isoforms that share a junction.
EvidenceVector build_evidence(const SqantiRecord& r,
                              const std::map<std::string, IntronChain>& chains,
                              const Catalog* catalog, const SjTable* sj, const BamReader* bam,
                              const std::vector<const BamReader*>& validation_bams,
                              const HaplotypeProvider* haplotype, bool haplotype_circular,
                              const PangenomeJunctions* pangenome, bool pangenome_circular,
                              const CallerSupport* caller_support,
                              const RuleConfig& cfg,
                              std::unordered_map<std::string, JunctionBamFeatures>& bam_cache,
                              std::unordered_map<std::string, VariantAt>& variant_cache,
                              std::unordered_map<std::string, JunctionRow>* jrows) {
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
  // intron chain (from a `combine` matrix). Joined by the chain fingerprint (two callers may
  // reuse a native id for different chains); monoexonic isoforms have no chain, so by native id.
  if (caller_support != nullptr) {
    const int* n = nullptr;
    if (ev.chain_available && !it->second.monoexonic()) {
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

  // Whole-chain long-read support: reads carrying every intron of the caller's chain. From the
  // --validation-bam libraries when given (each counted on its own; the largest count is the
  // isoform's, and a library supports the chain when it reaches chain_min_full_reads alone),
  // else from the caller's own BAM.
  if (ev.chain_available && !it->second.monoexonic()) {
    if (!validation_bams.empty()) {
      ev.chain_independent = true;
      // libraries are named by file name; when two share one, by the path given
      std::vector<std::string> names;
      bool clash = false;
      for (const BamReader* vb : validation_bams) {
        const std::string& p = vb->path();
        const std::size_t cut = p.rfind('/');
        names.push_back(cut == std::string::npos ? p : p.substr(cut + 1));
        for (std::size_t k = 0; k + 1 < names.size(); ++k) if (names[k] == names.back()) clash = true;
      }
      for (std::size_t li = 0; li < validation_bams.size(); ++li) {
        const BamReader* vb = validation_bams[li];
        const ChainBamSupport cs = vb->chain_support(it->second);
        ev.chain_full_by_library.emplace_back(clash ? vb->path() : names[li], cs.evaluated ? cs.n_full : -1);
        if (!cs.evaluated) continue;
        ev.chain_evaluable = true;
        if (cs.n_full >= cfg.chain_min_full_reads) ++ev.chain_n_libraries_supporting;
        if (cs.n_full > ev.chain_n_full_reads || (cs.n_full == ev.chain_n_full_reads && cs.n_partial > ev.chain_n_partial_reads)) {
          ev.chain_n_full_reads = cs.n_full;
          ev.chain_n_partial_reads = cs.n_partial;
          ev.chain_n_spanning_reads = cs.n_span;
        }
      }
      if (bam != nullptr) {
        const ChainBamSupport own = bam->chain_support(it->second);
        if (own.evaluated) ev.chain_n_full_reads_own = own.n_full;
      }
    } else if (bam != nullptr) {
      const ChainBamSupport cs = bam->chain_support(it->second);
      ev.chain_evaluable = cs.evaluated;
      ev.chain_n_full_reads = cs.n_full;
      ev.chain_n_partial_reads = cs.n_partial;
      ev.chain_n_spanning_reads = cs.n_span;
    }
  }

  // Novel-junction axes (short-read corroboration + BAM mapping features) require
  // the caller chain AND a catalog to know which junctions are novel.
  if (ev.catalog_evaluable) {
    const IntronChain& chain = it->second;
    BamFeatureParams bp;
    bp.min_mapq = cfg.bam_min_mapq;
    bp.junction_window_bp = cfg.bam_junction_window_bp;

    int n_novel = 0, n_supported = 0, n_pan = 0, n_variant_created = 0, n_variant_deleted = 0, n_variant_disrupted = 0;
    int n_phase_uncertain = 0, n_phase_contradicted = 0;
    // explained junctions that some but not all haplotypes explain rest on a heterozygous
    // allele; a common haplotype across them needs their alleles in one phase set
    int n_phase_dependent = 0; bool phase_broken = false; long phase_set_common = kNoPhaseSet;
    int n_known = 0, n_known_supported = 0, n_bam_gated = 0;
    bool any_bam_eval = false;
    bool any_variant_eval = false;
    bool any_pan_eval = false;
    // Haplotype/path consistency: a transcript is one molecule from one haplotype, so a
    // rescue needs one haplotype (or graph path) that explains every novel junction, not
    // a different one per junction. Intersect the per-junction sets.
    std::vector<bool> hap_common(haplotype != nullptr ? haplotype->n_haplotypes() : 0, true);
    std::vector<bool> hap_at;
    std::vector<std::string> path_common;   // paths shared by every supporting junction that listed ids
    int n_pan_with_paths = 0;
    auto sr_supported = [&](const auto& intron, JunctionRow* jr) {
      const SjRecord* hit = sj->find_exact(chain.chrom, chain.strand, intron);
      const bool ok = hit != nullptr && hit->n_uniq >= cfg.sj_min_uniq_reads &&
                      (!cfg.sj_require_canonical_motif || hit->canonical());
      if (jr != nullptr) { jr->sj_evaluable = true; jr->sj_n_uniq = hit ? hit->n_uniq : 0; jr->sr_supported = ok; }
      return ok;
    };
    for (const auto& intron : chain.introns) {
      const std::string key = junction_key(chain.chrom, chain.strand, intron);
      JunctionRow* jr = nullptr;
      if (jrows != nullptr) {
        jr = &(*jrows)[key];
        if (jr->isoforms.empty()) { jr->chrom = chain.chrom; jr->intron = intron; jr->strand = chain.strand; }
        jr->isoforms.push_back(r.isoform);
      }
      if (catalog->has_intron(chain.chrom, chain.strand, intron)) {
        // Known junction: its short-read coverage tells whether the locus is covered at all
        // (the rule engine treats "0 reads at the novel junction" as absence only then).
        ++n_known;
        if (jr != nullptr) jr->known = true;
        if (sj != nullptr && sr_supported(intron, jr)) ++n_known_supported;
        continue;
      }
      ++n_novel;
      if (haplotype != nullptr) {
        auto vit = variant_cache.find(key);
        if (vit == variant_cache.end())
          vit = variant_cache.emplace(key, variant_at(*haplotype, bam, chain.chrom, intron, chain.strand, cfg)).first;
        const VariantAt& va = vit->second;
        const JunctionMotifs& jm = va.jm;
        const VariantVerdict vv = jm.verdict;
        hap_at = va.on;
        if (jr != nullptr) {
          jr->variant_evaluable = vv != VariantVerdict::kNotEvaluable;
          jr->ref_motif = jm.ref.evaluable ? jm.ref.motif : "NA";
          jr->hap_motifs.clear();
          jr->canonical_on.clear();
          for (std::size_t h = 0; h < jm.haps.size(); ++h) {
            jr->hap_motifs.push_back(jm.deleted[h] ? "DEL" : jm.haps[h].evaluable ? jm.haps[h].motif : "NA");
            if (jm.haps[h].evaluable && jm.haps[h].canonical) jr->canonical_on.push_back(static_cast<int>(h) + 1);
          }
          jr->variant = vv == VariantVerdict::kCreated ? "created" : vv == VariantVerdict::kDeleted ? "deleted"
                      : vv == VariantVerdict::kDisrupted ? "disrupted" : vv == VariantVerdict::kNone ? "none" : "NA";
          jr->variant_phase = jm.phase_unlinked ? "unlinked" : jm.phase == MotifPhase::kPhased ? "phased" : jm.phase == MotifPhase::kUnphased ? "unphased" : "NA";
          jr->read_phase = va.read_phase;
          jr->read_phase_n_motif = va.reads.n_motif;
          jr->read_phase_n_other = va.reads.n_other;
        }
        if (vv != VariantVerdict::kNotEvaluable) any_variant_eval = true;
        if (vv == VariantVerdict::kDisrupted) ++n_variant_disrupted;
        // kCreated (non-canonical on reference, canonical on a haplotype) and kDeleted (a
        // haplotype lacks the intron's bases) drive the reference-bias rescue; kDisrupted
        // (canonical on reference, non-canonical on every haplotype) is a non-canonical sign
        // for the mechanism axis (rules.cpp).
        if (vv == VariantVerdict::kCreated || vv == VariantVerdict::kDeleted) {
          if (va.read_phase == "contradicted") {
            ++n_phase_contradicted;   // the reads come from a haplotype that does not explain the junction
          } else {
            if (vv == VariantVerdict::kCreated) ++n_variant_created; else ++n_variant_deleted;
            if (jm.phase == MotifPhase::kUnphased) ++n_phase_uncertain;
            for (std::size_t h = 0; h < hap_common.size(); ++h) hap_common[h] = hap_common[h] && hap_at[h];
            const bool on_every_hap = std::find(hap_at.begin(), hap_at.end(), false) == hap_at.end();
            if (!on_every_hap && jm.phase != MotifPhase::kUnknown) {   // VCF input and a heterozygous allele
              ++n_phase_dependent;
              if (jm.phase == MotifPhase::kUnphased || jm.phase_mixed) phase_broken = true;
              else if (phase_set_common == kNoPhaseSet) phase_set_common = jm.phase_set;
              else if (phase_set_common != jm.phase_set) phase_broken = true;
            }
          }
        }
      }
      if (pangenome != nullptr) {
        any_pan_eval = true;
        const PangenomeJunctionRecord* prec = pangenome->find_exact(chain.chrom, chain.strand, intron);
        const bool pan_ok = prec != nullptr && prec->n_haplotypes >= cfg.pangenome_min_haplotypes;
        if (jr != nullptr) {
          jr->pangenome_evaluable = true;
          jr->pangenome_supported = pan_ok;
          jr->pangenome_n_haplotypes = prec ? prec->n_haplotypes : 0;
          jr->pangenome_paths = prec ? prec->haplotype_ids : std::vector<std::string>{};
        }
        if (pan_ok) {
          ++n_pan;
          const auto& ids = prec->haplotype_ids;
          if (!ids.empty()) {
            if (n_pan_with_paths++ == 0) {
              path_common = ids;
            } else {
              std::vector<std::string> keep;
              for (const auto& id : path_common)
                if (std::find(ids.begin(), ids.end(), id) != ids.end()) keep.push_back(id);
              path_common.swap(keep);
            }
          }
        }
      }
      if (sj != nullptr && sr_supported(intron, jr)) ++n_supported;
      if (bam != nullptr) {
        auto cit = bam_cache.find(key);
        if (cit == bam_cache.end()) {
          cit = bam_cache.emplace(key, bam->features_at_junction(chain.chrom, intron, chain.strand, bp)).first;
        }
        const JunctionBamFeatures& f = cit->second;
        if (jr != nullptr) { jr->bam_evaluable = f.evaluated; jr->bam = f; }
        if (f.evaluated) {
          any_bam_eval = true;
          ev.bam_n_spanning_total += f.n_spanning;
          // Fractions are only meaningful with enough reads: with 1-2 reads a single indel
          // is a fraction of 1.0. Junctions below bam_min_spanning do not feed the gate.
          if (f.n_spanning >= cfg.bam_min_spanning) {
            ++n_bam_gated;
            if (f.frac_low_mapq() > ev.bam_max_frac_low_mapq) ev.bam_max_frac_low_mapq = f.frac_low_mapq();
            if (f.frac_supplementary() > ev.bam_max_frac_supplementary) ev.bam_max_frac_supplementary = f.frac_supplementary();
            if (f.frac_indel_near() > ev.bam_max_frac_indel_near) ev.bam_max_frac_indel_near = f.frac_indel_near();
            if (f.frac_indel_near_excess() > ev.bam_max_frac_indel_near_excess) ev.bam_max_frac_indel_near_excess = f.frac_indel_near_excess();
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
      ev.sj_no_junction_in_span = !chain.introns.empty() &&
          !sj->any_in_span(chain.chrom, chain.strand, chain.introns.front().start, chain.introns.back().end);
    }
    ev.bam_evaluable = (bam != nullptr) && any_bam_eval;
    ev.bam_n_junctions_gated = n_bam_gated;
    ev.variant_evaluable = (haplotype != nullptr) && any_variant_eval;
    ev.n_novel_jx_phase_uncertain = n_phase_uncertain;
    ev.n_novel_jx_variant_created = n_variant_created;
    ev.n_novel_jx_variant_deleted = n_variant_deleted;
    ev.n_novel_jx_read_phase_contradicted = n_phase_contradicted;
    ev.n_novel_jx_variant_disrupted = n_variant_disrupted;
    ev.haplotype_consistent = n_variant_created + n_variant_deleted > 0 &&
                              std::find(hap_common.begin(), hap_common.end(), true) != hap_common.end();
    ev.haplotype_phase_unknown = n_phase_dependent >= 2 && phase_broken;
    // Rescue requires all of the isoform's novel junctions to be explained by the
    // reference-bias evidence -- a single explained junction does not make a
    // multi-novel-junction isoform "false novel" (the others may be new) --
    // and explained on one haplotype (the per-junction count alone would accept junction
    // 1 on haplotype A plus junction 2 on haplotype B).
    ev.variant_rescue = (n_novel > 0 && n_variant_created + n_variant_deleted == n_novel && ev.haplotype_consistent);
    // Keep the circular-risk flag consistent with evaluability (it is only read by the
    // rescue branch, which already requires *_evaluable, so this is a clarity guard).
    ev.variant_circular = ev.variant_evaluable && haplotype_circular;
    ev.pangenome_evaluable = (pangenome != nullptr) && any_pan_eval;
    ev.n_novel_jx_pangenome = n_pan;
    // Path identity is optional in the junction file. The rows that carry it must admit one
    // path shared by all of them; a row without ids cannot confirm a path, but it cannot hide
    // a contradiction between the rows that have them either.
    ev.pangenome_path_evaluable = n_pan_with_paths > 0;
    ev.pangenome_path_consistent = ev.pangenome_path_evaluable && !path_common.empty();
    ev.pangenome_path_complete = n_pan > 0 && n_pan_with_paths == n_pan;
    ev.pangenome_rescue = (n_novel > 0 && n_pan == n_novel &&
                           (!ev.pangenome_path_evaluable || ev.pangenome_path_consistent));
    ev.pangenome_circular = ev.pangenome_evaluable && pangenome_circular;
  }
  return ev;
}

}  // namespace

void check_same_coordinates(const SqantiTable& sqanti, const std::map<std::string, IntronChain>& chains) {
  std::size_t n = 0;
  std::string example;
  for (const SqantiRecord& r : sqanti.records) {
    const auto it = chains.find(r.isoform);
    if (it == chains.end()) continue;
    const IntronChain& c = it->second;
    const bool strand_differs = r.strand != Strand::kUnknown && c.strand != Strand::kUnknown && r.strand != c.strand;
    if ((r.chrom.empty() || r.chrom == c.chrom) && !strand_differs) continue;
    if (n++ == 0)
      example = r.isoform + ": " + r.chrom + "/" + strand_char(r.strand) + " in the classification, " + c.chrom + "/" +
                strand_char(c.strand) + " in the isoform file";
  }
  if (n > 0)
    throw std::runtime_error(std::to_string(n) + " isoform(s) differ in chromosome or strand between the classification "
                             "and the isoform file (e.g. " + example + "); the two files are not from the same run");
}

std::vector<AdjudicationResult> adjudicate(const AdjudicateInputs& in, const RuleEngine& engine,
                                           std::vector<JunctionRow>* junctions) {
  if (in.sqanti == nullptr || in.chains == nullptr) {
    throw std::invalid_argument("adjudicate: sqanti and chains inputs are required");
  }

  std::unordered_map<std::string, JunctionBamFeatures> bam_cache;
  std::unordered_map<std::string, VariantAt> variant_cache;
  std::unordered_map<std::string, JunctionRow> jrows;
  std::vector<AdjudicationResult> out;
  out.reserve(in.sqanti->records.size());
  for (const auto& r : in.sqanti->records) {
    AdjudicationResult res;
    res.isoform_id = r.isoform;
    res.chrom = r.chrom;
    res.strand = r.strand;
    res.structural_category = r.structural_category;
    res.evidence = build_evidence(r, *in.chains, in.catalog, in.sj, in.bam, in.validation_bams, in.haplotype,
                                  in.haplotype_circular, in.pangenome, in.pangenome_circular,
                                  in.caller_support, engine.config(), bam_cache, variant_cache,
                                  junctions != nullptr ? &jrows : nullptr);
    res.verdict = engine.evaluate(res.evidence);
    out.push_back(std::move(res));
  }
  if (junctions != nullptr) {
    junctions->clear();
    junctions->reserve(jrows.size());
    for (auto& kv : jrows) junctions->push_back(std::move(kv.second));
    std::sort(junctions->begin(), junctions->end(), [](const JunctionRow& a, const JunctionRow& b) {
      if (a.chrom != b.chrom) return a.chrom < b.chrom;
      if (a.intron.start != b.intron.start) return a.intron.start < b.intron.start;
      if (a.intron.end != b.intron.end) return a.intron.end < b.intron.end;
      return strand_char(a.strand) < strand_char(b.strand);
    });
  }
  return out;
}

}  // namespace panisoguard
