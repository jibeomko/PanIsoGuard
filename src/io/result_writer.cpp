#include "panisoguard/result_writer.hpp"

#include "panisoguard/gtf.hpp"

#include <cstdio>
#include <ctime>
#include <fstream>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <sys/stat.h>

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

// DNA-side summary, independent of the verdict: do the sample's haplotypes / the pangenome
// make the isoform's novel junctions canonical / realizable? "haplotype" and "pangenome"
// mean every novel junction, on one haplotype or path; "partial" some of them; "none" the
// axis was evaluated and explains nothing; "NA" no genome evidence was supplied.
std::string genomic_compatibility(const EvidenceVector& e) {
  if (!e.variant_evaluable && !e.pangenome_evaluable) return "NA";
  std::string s;
  if (e.variant_evaluable && e.variant_rescue) s += "haplotype";
  if (e.pangenome_evaluable && e.pangenome_rescue) s += std::string(s.empty() ? "" : ";") + "pangenome";
  if (!s.empty()) return s;
  if (e.n_novel_jx_variant_created > 0 || e.n_novel_jx_pangenome > 0) return "partial";
  return "none";
}

void write_tsv(const std::string& path, const std::vector<AdjudicationResult>& results) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << "isoform_id\tchrom\tstrand\tstructural_category\tnovelty_support\tprimary_mechanism\t"
         "confidence_class\tverdict\tverdict_reason\tn_novel_junctions\tn_novel_jx_sr_supported\t"
         "circularity_flag\tgenomic_compatibility\n";
  for (const auto& r : results) {
    out << r.isoform_id << '\t' << r.chrom << '\t' << strand_char(r.strand) << '\t'
        << r.structural_category << '\t' << to_string(r.verdict.novelty_support) << '\t'
        << to_string(r.verdict.primary_mechanism) << '\t' << to_string(r.verdict.confidence)
        << '\t' << r.verdict.verdict << '\t' << r.verdict.verdict_reason
        << '\t' << r.evidence.n_novel_junctions << '\t' << r.evidence.n_novel_jx_sr_supported
        << '\t' << (r.verdict.circularity_flag ? "true" : "false")
        << '\t' << genomic_compatibility(r.evidence) << '\n';
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
        << "\"verdict\":\"" << r.verdict.verdict << "\","
        << "\"verdict_reason\":\"" << json_escape(r.verdict.verdict_reason) << "\","
        << "\"novelty_support\":\"" << to_string(r.verdict.novelty_support) << "\","
        << "\"primary_mechanism\":\"" << to_string(r.verdict.primary_mechanism) << "\","
        << "\"circularity_flag\":" << (r.verdict.circularity_flag ? "true" : "false") << ","
        << "\"genomic_compatibility\":\"" << genomic_compatibility(e) << "\","
        << "\"evidence\":{"
        << "\"chain_available\":" << (e.chain_available ? "true" : "false") << ","
        << "\"sj_evaluable\":" << (e.sj_evaluable ? "true" : "false") << ","
        << "\"n_novel_junctions\":" << e.n_novel_junctions << ","
        << "\"n_novel_jx_sr_supported\":" << e.n_novel_jx_sr_supported << ","
        << "\"n_known_junctions\":" << e.n_known_junctions << ","
        << "\"n_known_jx_sr_supported\":" << e.n_known_jx_sr_supported << ","
        << "\"rts_stage\":" << (e.rts_evaluable ? (e.rts_stage ? "true" : "false") : "null") << ","
        << "\"noncanonical\":" << (e.canon_evaluable ? (e.noncanonical ? "true" : "false") : "null") << ","
        << "\"polyA_motif_found\":" << (e.polyA_motif_evaluable ? (e.polyA_motif_found ? "true" : "false") : "null") << ","
        << "\"bam_evaluable\":" << (e.bam_evaluable ? "true" : "false") << ","
        << "\"bam_n_spanning\":" << e.bam_n_spanning_total << ","
        << "\"bam_n_junctions_gated\":" << e.bam_n_junctions_gated << ","
        << "\"bam_frac_low_mapq\":" << e.bam_max_frac_low_mapq << ","
        << "\"bam_frac_supplementary\":" << e.bam_max_frac_supplementary << ","
        << "\"bam_frac_indel_near\":" << e.bam_max_frac_indel_near << ","
        << "\"variant_evaluable\":" << (e.variant_evaluable ? "true" : "false") << ","
        << "\"variant_rescue\":" << (e.variant_rescue ? "true" : "false") << ","
        << "\"n_novel_jx_variant_created\":" << e.n_novel_jx_variant_created << ","
        << "\"n_novel_jx_variant_disrupted\":" << e.n_novel_jx_variant_disrupted << ","
        << "\"haplotype_consistent\":" << (e.n_novel_jx_variant_created > 0 ? (e.haplotype_consistent ? "true" : "false") : "null") << ","
        << "\"variant_circular\":" << (e.variant_circular ? "true" : "false") << ","
        << "\"pangenome_evaluable\":" << (e.pangenome_evaluable ? "true" : "false") << ","
        << "\"pangenome_rescue\":" << (e.pangenome_rescue ? "true" : "false") << ","
        << "\"pangenome_circular\":" << (e.pangenome_circular ? "true" : "false") << ","
        << "\"n_novel_jx_pangenome\":" << e.n_novel_jx_pangenome << ","
        << "\"pangenome_path_consistent\":" << (e.pangenome_path_evaluable ? (e.pangenome_path_consistent ? "true" : "false") : "null") << ","
        << "\"consensus_evaluable\":" << (e.consensus_evaluable ? "true" : "false") << ","
        << "\"n_callers\":" << (e.consensus_evaluable ? std::to_string(e.n_callers) : "null") << ","
        << "\"chain_evaluable\":" << (e.chain_evaluable ? "true" : "false") << ","
        << "\"chain_independent\":" << (e.chain_independent ? "true" : "false") << ","
        << "\"chain_n_full_reads\":" << (e.chain_evaluable ? std::to_string(e.chain_n_full_reads) : "null") << ","
        << "\"chain_n_partial_reads\":" << (e.chain_evaluable ? std::to_string(e.chain_n_partial_reads) : "null")
        << "},";
    out << "\"rule_trace\":[";
    for (std::size_t i = 0; i < r.verdict.rule_trace.size(); ++i) {
      if (i) out << ',';
      out << '"' << json_escape(r.verdict.rule_trace[i]) << '"';
    }
    out << "]}" << '\n';
  }
}

// "size=<bytes>\tmtime=<UTC>" for one input file, so a re-run can tell whether it saw the
// same file.
std::string file_identity(const std::string& path) {
  struct stat st;
  if (stat(path.c_str(), &st) != 0) return "size=NA\tmtime=NA";
  char when[32] = "NA";
  struct tm tmv;
  if (gmtime_r(&st.st_mtime, &tmv) != nullptr) std::strftime(when, sizeof(when), "%Y-%m-%dT%H:%M:%SZ", &tmv);
  return "size=" + std::to_string(static_cast<long long>(st.st_size)) + "\tmtime=" + when;
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
  if (!prov.thresholds.empty()) out << "thresholds\t" << prov.thresholds << '\n';
  out << "classification\t" << prov.classification_path << '\n';
  out << "isoforms\t" << prov.isoforms_path << '\n';
  out << "ref_gtf\t" << (prov.ref_gtf_path.empty() ? "<none>" : prov.ref_gtf_path) << '\n';
  out << "sj_tab\t" << (prov.sj_tab_path.empty() ? "<none>" : prov.sj_tab_path) << '\n';
  out << "bam\t" << (prov.bam_path.empty() ? "<none>" : prov.bam_path) << '\n';
  out << "validation_bam\t" << (prov.validation_bam_path.empty() ? "<none>" : prov.validation_bam_path) << '\n';
  out << "caller_support\t"
      << (prov.caller_support_path.empty() ? "<none>" : prov.caller_support_path) << '\n';
  out << "reference\t" << (prov.reference_path.empty() ? "<none>" : prov.reference_path) << '\n';
  for (const auto& h : prov.reference_haplotype_paths) out << "reference_haplotype\t" << h << '\n';
  out << "pangenome_junctions\t" << (prov.pangenome_path.empty() ? "<none>" : prov.pangenome_path) << '\n';
  if (!prov.kept_path.empty()) {
    out << "kept_isoforms\t" << prov.kept_path << '\n';
    out << "kept_verdicts\t" << prov.kept_verdicts << '\n';
  }
  // Input identity (size, mtime): the paths above say what was named, these say what was
  // read. Lines start with "input." so an expected-output diff can skip them.
  out << "# inputs\n";
  auto ident = [&](const char* name, const std::string& p) {
    if (!p.empty()) out << "input." << name << '\t' << p << '\t' << file_identity(p) << '\n';
  };
  ident("classification", prov.classification_path);
  ident("isoforms", prov.isoforms_path);
  ident("ref_gtf", prov.ref_gtf_path);
  ident("sj_tab", prov.sj_tab_path);
  ident("bam", prov.bam_path);
  ident("validation_bam", prov.validation_bam_path);
  ident("caller_support", prov.caller_support_path);
  ident("reference", prov.reference_path);
  for (const auto& h : prov.reference_haplotype_paths) ident("reference_haplotype", h);
  ident("pangenome_junctions", prov.pangenome_path);
  ident("config", prov.config_path);
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
  out << "axis.consensus\t" << (prov.consensus_axis_on ? "on" : "not_evaluable") << '\n';
  out << "axis.chain\t"
      << (!prov.validation_bam_path.empty() ? "on (independent library: --validation-bam)"
          : prov.bam_path.empty() ? "not_evaluable" : "on (same reads as the caller)") << '\n';

  std::map<std::string, std::size_t> class_counts, verdict_counts;
  for (const auto& r : results) {
    ++class_counts[to_string(r.verdict.confidence)];
    ++verdict_counts[r.verdict.verdict];
  }
  out << "# confidence-class counts\n";
  for (const auto& kv : class_counts) out << "class." << kv.first << '\t' << kv.second << '\n';
  out << "# verdict counts\n";
  for (const auto& kv : verdict_counts) out << "verdict." << kv.first << '\t' << kv.second << '\n';
  out << "total\t" << results.size() << '\n';
}

}  // namespace

void write_junctions_tsv(const std::string& path, const std::vector<JunctionRow>& rows) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  auto join = [](const auto& v, const char* sep = ",") {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) s += sep; s += [](const auto& x) {
      if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::string>) return x; else return std::to_string(x); }(v[i]); }
    return s.empty() ? std::string("-") : s;
  };
  auto yn = [](bool evaluable, bool v) { return evaluable ? (v ? "true" : "false") : "NA"; };
  out << "chrom\tstart\tend\tstrand\tn_isoforms\tisoforms\tknown\tsr_n_uniq\tsr_supported\tref_motif\thap_motifs\t"
         "variant\tcanonical_on_haplotypes\tpangenome_n_haplotypes\tpangenome_paths\tpangenome_supported\t"
         "bam_n_spanning\tbam_frac_low_mapq\tbam_frac_supplementary\tbam_frac_indel_near\treference_bias\n";
  for (const auto& j : rows) {
    std::string rb;
    if (j.variant == "created") rb += "haplotype_variant";
    if (j.pangenome_supported) rb += std::string(rb.empty() ? "" : ";") + "pangenome_path";
    out << j.chrom << '\t' << j.intron.start + 1 << '\t' << j.intron.end << '\t' << strand_char(j.strand) << '\t'
        << j.isoforms.size() << '\t' << join(j.isoforms) << '\t' << (j.known ? "true" : "false") << '\t'
        << (j.sj_evaluable ? std::to_string(j.sj_n_uniq) : "NA") << '\t' << yn(j.sj_evaluable, j.sr_supported) << '\t'
        << (j.variant_evaluable ? j.ref_motif : "NA") << '\t' << (j.variant_evaluable ? join(j.hap_motifs) : "NA") << '\t'
        << j.variant << '\t' << (j.variant_evaluable ? join(j.canonical_on) : "NA") << '\t'
        << (j.pangenome_evaluable ? std::to_string(j.pangenome_n_haplotypes) : "NA") << '\t'
        << (j.pangenome_evaluable ? join(j.pangenome_paths) : "NA") << '\t' << yn(j.pangenome_evaluable, j.pangenome_supported) << '\t';
    if (j.bam_evaluable) {
      out << j.bam.n_spanning << '\t' << j.bam.frac_low_mapq() << '\t' << j.bam.frac_supplementary() << '\t'
          << j.bam.frac_indel_near();
    } else {
      out << "NA\tNA\tNA\tNA";
    }
    out << '\t' << (rb.empty() ? "-" : rb) << '\n';
  }
}

std::size_t write_kept_isoforms(const std::string& in_path, bool is_bed, const std::string& out_path,
                                const std::unordered_set<std::string>& keep) {
  std::ifstream in(in_path);
  if (!in) throw std::runtime_error("cannot open " + in_path);
  std::ofstream out(out_path);
  if (!out) throw std::runtime_error("cannot write " + out_path);
  std::unordered_set<std::string> written;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] == '#' || (is_bed && (line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0))) {
      if (!is_bed) out << line << '\n';
      continue;
    }
    std::string id;
    if (is_bed) {
      std::size_t a = 0;
      for (int col = 0; col < 3; ++col) { a = line.find('\t', a); if (a == std::string::npos) break; ++a; }
      if (a == std::string::npos) continue;
      id = line.substr(a, line.find('\t', a) - a);
    } else {
      const std::size_t tab = line.rfind('\t');
      if (tab == std::string::npos) continue;
      id = extract_gtf_attr(line.substr(tab + 1), "transcript_id");
    }
    if (id.empty() || keep.count(id) == 0) continue;
    out << line << '\n';
    written.insert(id);
  }
  return written.size();
}

void write_adjudication_outputs(const std::string& out_prefix,
                                const std::vector<AdjudicationResult>& results,
                                const RunProvenance& prov) {
  write_tsv(out_prefix + ".adjudicated.tsv", results);
  write_jsonl(out_prefix + ".attribution.jsonl", results);
  write_provenance(out_prefix + ".provenance.log", results, prov);
}

}  // namespace panisoguard
