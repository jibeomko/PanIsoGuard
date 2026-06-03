#include "panisoguard/result_writer.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>

namespace panisoguard {
namespace {

std::string json_escape(const std::string& s) {
  std::string o;
  o.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '"':  o += "\\\""; break;
      case '\\': o += "\\\\"; break;
      case '\n': o += "\\n";  break;
      case '\r': o += "\\r";  break;
      case '\t': o += "\\t";  break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", c);
          o += buf;
        } else {
          o += c;
        }
    }
  }
  return o;
}

void write_tsv(const std::string& path, const std::vector<AdjudicationResult>& results) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << "isoform_id\tchrom\tstrand\tstructural_category\tnovelty_support\tprimary_mechanism\t"
         "confidence_class\tn_novel_junctions\tn_novel_jx_sr_supported\tcircularity_flag\n";
  for (const auto& r : results) {
    out << r.isoform_id << '\t' << r.chrom << '\t' << strand_char(r.strand) << '\t'
        << r.structural_category << '\t' << to_string(r.verdict.novelty_support) << '\t'
        << to_string(r.verdict.primary_mechanism) << '\t' << to_string(r.verdict.confidence)
        << '\t' << r.evidence.n_novel_junctions << '\t' << r.evidence.n_novel_jx_sr_supported
        << '\t' << (r.verdict.circularity_flag ? "true" : "false") << '\n';
  }
}

void write_jsonl(const std::string& path, const std::vector<AdjudicationResult>& results) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  for (const auto& r : results) {
    const EvidenceVector& e = r.evidence;
    out << '{'
        << "\"isoform\":\"" << json_escape(r.isoform_id) << "\","
        << "\"structural_category\":\"" << json_escape(r.structural_category) << "\","
        << "\"confidence_class\":\"" << to_string(r.verdict.confidence) << "\","
        << "\"novelty_support\":\"" << to_string(r.verdict.novelty_support) << "\","
        << "\"primary_mechanism\":\"" << to_string(r.verdict.primary_mechanism) << "\","
        << "\"circularity_flag\":" << (r.verdict.circularity_flag ? "true" : "false") << ","
        << "\"evidence\":{"
        << "\"chain_available\":" << (e.chain_available ? "true" : "false") << ","
        << "\"sj_evaluable\":" << (e.sj_evaluable ? "true" : "false") << ","
        << "\"n_novel_junctions\":" << e.n_novel_junctions << ","
        << "\"n_novel_jx_sr_supported\":" << e.n_novel_jx_sr_supported << ","
        << "\"rts_stage\":" << (e.rts_evaluable ? (e.rts_stage ? "true" : "false") : "null") << ","
        << "\"noncanonical\":" << (e.canon_evaluable ? (e.noncanonical ? "true" : "false") : "null") << ","
        << "\"bam_evaluable\":" << (e.bam_evaluable ? "true" : "false") << ","
        << "\"bam_n_spanning\":" << e.bam_n_spanning_total << ","
        << "\"bam_frac_low_mapq\":" << e.bam_max_frac_low_mapq << ","
        << "\"bam_frac_supplementary\":" << e.bam_max_frac_supplementary << ","
        << "\"bam_frac_softclip\":" << e.bam_max_frac_softclip << ","
        << "\"bam_frac_indel_near\":" << e.bam_max_frac_indel_near << ","
        << "\"variant_evaluable\":" << (e.variant_evaluable ? "true" : "false") << ","
        << "\"variant_rescue\":" << (e.variant_rescue ? "true" : "false") << ","
        << "\"variant_circular\":" << (e.variant_circular ? "true" : "false") << ","
        << "\"pangenome_evaluable\":" << (e.pangenome_evaluable ? "true" : "false") << ","
        << "\"pangenome_rescue\":" << (e.pangenome_rescue ? "true" : "false") << ","
        << "\"pangenome_circular\":" << (e.pangenome_circular ? "true" : "false") << ","
        << "\"n_novel_jx_pangenome\":" << e.n_novel_jx_pangenome << ","
        << "\"consensus_evaluable\":" << (e.consensus_evaluable ? "true" : "false") << ","
        << "\"n_callers\":" << (e.consensus_evaluable ? std::to_string(e.n_callers) : "null")
        << "},";
    // graph_trace: the same evidence re-expressed under the splice-graph framing
    // (docs/method_graph.md) -- a candidate isoform is a path through the gene-local
    // splice graph; novel junctions are edges absent from the known (reference) graph.
    // No new computation: this is a view of the EvidenceVector above.
    const char* site_type = "none";
    if (r.structural_category == "novel_in_catalog") site_type = "known_site_recombination";
    else if (r.structural_category == "novel_not_in_catalog") site_type = "novel_splice_site";
    out << "\"graph_trace\":{"
        << "\"novel_edges\":" << e.n_novel_junctions << ","
        << "\"path_distance_to_reference\":" << e.n_novel_junctions << ","
        << "\"edges_sr_supported\":" << e.n_novel_jx_sr_supported << ","
        << "\"novel_site_type\":\"" << site_type << "\","
        << "\"pangenome_edge_support\":" << (e.pangenome_rescue ? "true" : "false") << ","
        << "\"variant_canonicalized_edge\":" << (e.variant_rescue ? "true" : "false")
        << "},";
    // bio_flags: SQANTI3 QC descriptors passed through verbatim (consumed as priors,
    // NEVER recomputed) -- biological-plausibility annotations that do NOT affect the
    // verdict. PanIsoGuard does not reimplement SQANTI3's TSS/TTS/ORF/NMD/polyA QC; it
    // surfaces it alongside the independent path-level verdict (docs/relationship_to_sqanti3.md).
    out << "\"bio_flags\":{\"dist_to_CAGE_peak\":";
    if (sqanti_is_na(r.bio_dist_to_CAGE_peak)) out << "null"; else out << r.bio_dist_to_CAGE_peak;
    out << ",\"dist_to_polyA_site\":";
    if (sqanti_is_na(r.bio_dist_to_polyA_site)) out << "null"; else out << r.bio_dist_to_polyA_site;
    out << ",\"predicted_NMD\":"
        << (r.bio_predicted_NMD.empty() ? "null" : "\"" + json_escape(r.bio_predicted_NMD) + "\"")
        << ",\"polyA_motif_found\":"
        << (r.bio_polyA_motif_found.empty() ? "null" : "\"" + json_escape(r.bio_polyA_motif_found) + "\"")
        << "},";
    out << "\"rule_trace\":[";
    for (std::size_t i = 0; i < r.verdict.rule_trace.size(); ++i) {
      if (i) out << ',';
      out << '"' << json_escape(r.verdict.rule_trace[i]) << '"';
    }
    out << "]}" << '\n';
  }
}

void write_provenance(const std::string& path, const std::vector<AdjudicationResult>& results,
                      const RunProvenance& prov) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << "# PanIsoGuard adjudicate provenance\n";
  out << "tool_version\t" << prov.tool_version << '\n';
  out << "ruleset_version\t" << prov.ruleset_version << '\n';
  out << "sqanti3_version_target\t" << prov.sqanti3_version_target << '\n';
  out << "config\t" << (prov.config_path.empty() ? "<built-in defaults>" : prov.config_path) << '\n';
  out << "classification\t" << prov.classification_path << '\n';
  out << "isoforms\t" << prov.isoforms_path << '\n';
  out << "ref_gtf\t" << (prov.ref_gtf_path.empty() ? "<none>" : prov.ref_gtf_path) << '\n';
  out << "sj_tab\t" << (prov.sj_tab_path.empty() ? "<none>" : prov.sj_tab_path) << '\n';
  out << "bam\t" << (prov.bam_path.empty() ? "<none>" : prov.bam_path) << '\n';
  // Evidence-axis capabilities for this run.
  out << "axis.short_read\t" << (prov.sj_tab_path.empty() ? "not_evaluable" : "on") << '\n';
  out << "axis.catalog\t" << (prov.ref_gtf_path.empty() ? "not_evaluable" : "on") << '\n';
  out << "axis.bam\t" << (prov.bam_path.empty() ? "not_evaluable" : "on") << '\n';
  out << "axis.variant\t"
      << (prov.variant_axis_on ? (prov.variant_circular ? "on (circular-risk: held, not promoted)" : "on")
                               : "not_evaluable")
      << '\n';
  out << "axis.pangenome\t"
      << (prov.pangenome_axis_on
              ? (prov.pangenome_circular ? "on (circular-risk: held, not promoted)" : "on")
              : "not_evaluable")
      << '\n';

  std::map<std::string, std::size_t> class_counts;
  for (const auto& r : results) ++class_counts[to_string(r.verdict.confidence)];
  out << "# confidence-class counts\n";
  for (const auto& kv : class_counts) out << "class." << kv.first << '\t' << kv.second << '\n';
  out << "total\t" << results.size() << '\n';
}

}  // namespace

void write_adjudication_outputs(const std::string& out_prefix,
                                const std::vector<AdjudicationResult>& results,
                                const RunProvenance& prov) {
  write_tsv(out_prefix + ".adjudicated.tsv", results);
  write_jsonl(out_prefix + ".attribution.jsonl", results);
  write_provenance(out_prefix + ".provenance.log", results, prov);
}

}  // namespace panisoguard
