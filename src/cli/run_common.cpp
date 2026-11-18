#include "run_common.hpp"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "panisoguard/bed12.hpp"
#include "panisoguard/classify.hpp"
#include "panisoguard/consensus.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/pangenome.hpp"
#include "panisoguard/sj_tab.hpp"
#include "panisoguard/tsv.hpp"
#include "panisoguard/variant_motif.hpp"
#include "panisoguard/version.hpp"

namespace panisoguard {

bool consume_common_arg(const std::string& a, int& i, int argc, char** argv, CommonArgs& c) {
  auto next = [&](const char* name) -> std::string {
    if (i + 1 >= argc) {
      std::fprintf(stderr, "panisoguard: %s requires an argument\n", name);
      std::exit(1);
    }
    return argv[++i];
  };
  if (a == "--classification")     { c.classification = next("--classification"); return true; }
  if (a == "--isoforms-bed")       { c.isoforms_bed = next("--isoforms-bed"); return true; }
  if (a == "--isoforms-gtf")       { c.isoforms_gtf = next("--isoforms-gtf"); return true; }
  if (a == "--ref-gtf")            { c.ref_gtf = next("--ref-gtf"); return true; }
  if (a == "--sj-tab")             { c.sj_tab = next("--sj-tab"); return true; }
  if (a == "--bam")                { c.bam_path = next("--bam"); return true; }
  if (a == "--validation-bam")     { c.validation_bams.push_back(next("--validation-bam")); return true; }
  if (a == "--reference")          { c.reference = next("--reference"); return true; }
  if (a == "--reference-haplotype"){ c.reference_haplotypes.push_back(next("--reference-haplotype")); return true; }
  if (a == "--haplotype-vcf")      { c.haplotype_vcfs.push_back(next("--haplotype-vcf")); return true; }
  if (a == "--vcf-sample")         { c.vcf_sample = next("--vcf-sample"); return true; }
  if (a == "--haplotype-provenance"){ c.haplotype_provenance = next("--haplotype-provenance"); return true; }
  if (a == "--pangenome-junctions"){ c.pangenome_junctions = next("--pangenome-junctions"); return true; }
  if (a == "--pangenome-provenance"){ c.pangenome_provenance = next("--pangenome-provenance"); return true; }
  if (a == "--caller-support")     { c.caller_support = next("--caller-support"); return true; }
  if (a == "--config")             { c.config = next("--config"); return true; }
  if (a == "--platform")           { c.platform = next("--platform"); return true; }
  if (a == "--no-sqanti3")         { c.no_sqanti3 = true; return true; }
  return false;
}

std::vector<std::string> input_paths(const CommonArgs& c) {
  std::vector<std::string> v = {c.classification, c.isoforms_bed, c.isoforms_gtf, c.ref_gtf, c.sj_tab, c.bam_path,
                                c.reference, c.pangenome_junctions, c.caller_support, c.config};
  v.insert(v.end(), c.validation_bams.begin(), c.validation_bams.end());
  v.insert(v.end(), c.reference_haplotypes.begin(), c.reference_haplotypes.end());
  v.insert(v.end(), c.haplotype_vcfs.begin(), c.haplotype_vcfs.end());
  return v;
}

bool common_args_ok(const CommonArgs& c, std::string& err) {
  if (c.no_sqanti3 && !c.classification.empty()) {
    err = "--no-sqanti3 computes the classification; drop --classification or drop the flag";
    return false;
  }
  if (c.no_sqanti3 && c.ref_gtf.empty()) { err = "--no-sqanti3 needs --ref-gtf to classify the isoforms"; return false; }
  if (!c.no_sqanti3 && c.classification.empty()) { err = "--classification is required (or --no-sqanti3)"; return false; }
  if (c.platform != "hifi" && c.platform != "ont") {
    err = "--platform must be hifi or ont (got '" + c.platform + "')";
    return false;
  }
  if (c.isoforms_bed.empty() && c.isoforms_gtf.empty()) {
    err = "one of --isoforms-bed / --isoforms-gtf is required";
    return false;
  }
  if (!c.isoforms_bed.empty() && !c.isoforms_gtf.empty()) {
    err = "--isoforms-bed and --isoforms-gtf are mutually exclusive; supply exactly one";
    return false;
  }
  if (!c.haplotype_vcfs.empty() && !c.reference_haplotypes.empty()) {
    err = "--haplotype-vcf and --reference-haplotype describe the same haplotypes two ways; supply one of them";
    return false;
  }
  if (!c.vcf_sample.empty() && c.haplotype_vcfs.empty()) {
    err = "--vcf-sample needs --haplotype-vcf";
    return false;
  }
  return true;
}

// Ids are joined by exact string match; a mismatch silently degrades verdicts
// (caller_chain_absent / no consensus), so say how many failed to join.
template <typename Map>
void warn_unjoined(const SqantiTable& sqanti, const Map& by_id, bool novel_only, const char* what) {
  std::size_t n = 0, missing = 0;
  std::string example;
  for (const auto& r : sqanti.records) {
    if (novel_only && !r.is_novel()) continue;
    ++n;
    if (by_id.count(r.isoform) == 0 && missing++ == 0) example = r.isoform;
  }
  if (missing > 0) {
    std::fprintf(stderr,
        "WARNING: %zu of %zu %sSQANTI3 isoform id(s) (e.g. \"%s\") are not in %s; "
        "check that the ids match exactly.\n",
        missing, n, novel_only ? "novel " : "", example.c_str(), what);
  }
}

LoadedRun load_and_adjudicate(const CommonArgs& c) {
  LoadedRun run;
  std::vector<CallerIsoform> isoforms;
  if (!c.isoforms_bed.empty()) {
    for (auto& r : read_bed12(c.isoforms_bed)) isoforms.push_back({r.name, std::move(r.chain), r.start, r.end});
  } else {
    for (auto& t : read_gtf_transcripts(c.isoforms_gtf)) isoforms.push_back({t.id, std::move(t.chain), t.start, t.end});
  }
  std::map<std::string, IntronChain> chains;
  for (const auto& iso : isoforms) chains.emplace(iso.id, iso.chain);
  std::fprintf(stderr, "read %zu caller isoform chains\n", chains.size());

  Catalog catalog;
  const Catalog* catalog_ptr = nullptr;
  if (!c.ref_gtf.empty()) {
    std::fprintf(stderr, "building reference catalog...\n");
    catalog = build_catalog_from_gtf(c.ref_gtf, /*with_transcripts=*/c.no_sqanti3);   // only the classifier reads the transcripts
    catalog_ptr = &catalog;
  }

  std::shared_ptr<FastaFetcher> ref_fa;
  if (c.no_sqanti3) {
    if (!c.reference.empty()) ref_fa = std::make_shared<FastaFetcher>(c.reference);
    run.sqanti = classify_isoforms(isoforms, catalog, ref_fa.get());
    std::map<std::string, std::size_t> by_cat;
    for (const auto& r : run.sqanti.records) ++by_cat[r.structural_category];
    std::string summary;
    for (const auto& kv : by_cat) summary += (summary.empty() ? "" : ", ") + kv.first + " " + std::to_string(kv.second);
    std::fprintf(stderr, "classified %zu isoforms without SQANTI3 (%s)%s\n", run.sqanti.records.size(), summary.c_str(),
                 ref_fa ? "" : "; no --reference: splice motifs and the 3' A content are not evaluable");
  } else {
    run.sqanti = read_sqanti_classification(c.classification);
    std::fprintf(stderr, "read %zu SQANTI records (%zu cols)\n", run.sqanti.records.size(), run.sqanti.n_columns);
    warn_unjoined(run.sqanti, chains, /*novel_only=*/false, "the caller isoform file");
    check_same_coordinates(run.sqanti, chains);
  }

  if (catalog_ptr == nullptr) {
    // Without a reference catalog there is no way to tell which junctions are novel,
    // so the whole novel-junction block is skipped: every novel isoform is held
    // UNKNOWN/UNRESOLVED and the short-read/BAM/variant/pangenome axes contribute
    // NOTHING even when their inputs are supplied. Warn loudly rather than fail soft.
    const bool novel_axes_supplied = !c.sj_tab.empty() || !c.bam_path.empty() || !c.validation_bams.empty() ||
                                     !c.reference_haplotypes.empty() || !c.haplotype_vcfs.empty() ||
                                     !c.pangenome_junctions.empty();
    std::fprintf(stderr,
        "WARNING: --ref-gtf not supplied: novel-junction classification is disabled; "
        "every novel isoform is held UNRESOLVED%s.\n",
        novel_axes_supplied ? " and the supplied SJ/BAM/haplotype/pangenome evidence "
                              "will NOT be used" : "");
  }

  SjTable sj;
  const SjTable* sj_ptr = nullptr;
  if (!c.sj_tab.empty()) {
    sj = read_sj_tab(c.sj_tab);
    std::fprintf(stderr, "read %zu short-read junctions\n", sj.size());
    if (sj.size() == 0) {
      // With no junctions at all, "no read at this junction" says nothing: the axis stays off
      // rather than marking every novel junction unsupported.
      std::fprintf(stderr, "WARNING: %s holds no junctions; the short-read axis is not evaluable in this run\n", c.sj_tab.c_str());
    } else {
      sj_ptr = &sj;
    }
  }

  std::unique_ptr<BamReader> bam;
  if (!c.bam_path.empty()) {
    bam = std::make_unique<BamReader>(c.bam_path, c.reference);
    std::fprintf(stderr, "opened BAM for mapping axis: %s\n", c.bam_path.c_str());
  }
  std::vector<std::unique_ptr<BamReader>> validation_bams;
  for (std::size_t i = 0; i < c.validation_bams.size(); ++i) {
    const std::string& vb = c.validation_bams[i];
    if (!c.bam_path.empty() && same_file(c.bam_path, vb))
      throw std::runtime_error("--validation-bam " + vb + " is the same file as --bam; the caller's own "
                               "reads are not an independent library (drop --validation-bam to count them as the caller's)");
    for (std::size_t j = 0; j < i; ++j)
      if (same_file(c.validation_bams[j], vb))
        throw std::runtime_error("--validation-bam " + vb + " and " + c.validation_bams[j] + " are the same file");
    validation_bams.push_back(std::make_unique<BamReader>(vb, c.reference));
    std::fprintf(stderr, "opened validation BAM %zu/%zu for whole-chain read support (independent library): %s\n",
                 i + 1, c.validation_bams.size(), vb.c_str());
  }

  std::vector<std::shared_ptr<FastaFetcher>> hap_fas;
  std::unique_ptr<HaplotypeProvider> haplo;
  if (!c.reference_haplotypes.empty()) {
    if (c.reference.empty()) {
      throw std::runtime_error("--reference (genome FASTA) is required with --reference-haplotype");
    }
    if (!ref_fa) ref_fa = std::make_shared<FastaFetcher>(c.reference);
    for (const auto& h : c.reference_haplotypes) hap_fas.push_back(std::make_shared<FastaFetcher>(h));
    haplo = std::make_unique<HaplotypeProvider>(ref_fa, hap_fas);
    std::fprintf(stderr, "variant axis: reference + %zu haplotype FASTA(s), provenance=%s\n",
                 hap_fas.size(), c.haplotype_provenance.c_str());
  } else if (!c.haplotype_vcfs.empty()) {
    if (c.reference.empty()) {
      throw std::runtime_error("--reference (genome FASTA) is required with --haplotype-vcf");
    }
    if (!ref_fa) ref_fa = std::make_shared<FastaFetcher>(c.reference);
    for (std::size_t i = 0; i < c.haplotype_vcfs.size(); ++i)
      for (std::size_t j = 0; j < i; ++j)
        if (same_file(c.haplotype_vcfs[j], c.haplotype_vcfs[i]))
          throw std::runtime_error("--haplotype-vcf " + c.haplotype_vcfs[i] + " and " + c.haplotype_vcfs[j] + " are the same file");
    for (const auto& h : c.haplotype_vcfs) run.haplotype_vcfs.push_back(std::make_shared<VcfHaplotypes>(h, c.vcf_sample));
    haplo = std::make_unique<HaplotypeProvider>(ref_fa, run.haplotype_vcfs);
    std::fprintf(stderr, "variant axis: reference + %zu haplotype(s) from %zu VCF(s)%s%s, provenance=%s\n",
                 haplo->n_haplotypes(), run.haplotype_vcfs.size(), c.vcf_sample.empty() ? "" : " sample ",
                 c.vcf_sample.c_str(), c.haplotype_provenance.c_str());
  }

  PangenomeJunctions pangenome;
  const PangenomeJunctions* pangenome_ptr = nullptr;
  if (!c.pangenome_junctions.empty()) {
    pangenome = read_pangenome_junctions(c.pangenome_junctions);
    pangenome_ptr = &pangenome;
    std::fprintf(stderr, "read %zu pangenome graph junctions, provenance=%s\n",
                 pangenome.size(), c.pangenome_provenance.c_str());
  }

  CallerSupport caller_support;
  const CallerSupport* caller_support_ptr = nullptr;
  if (!c.caller_support.empty()) {
    caller_support = read_caller_support(c.caller_support);
    caller_support_ptr = &caller_support;
    std::fprintf(stderr, "read caller-support for %zu chain(s) / %zu native isoform id(s) (consensus axis, "
                 "joined by chain fingerprint)\n", caller_support.by_chain.size(), caller_support.by_native_id.size());
    warn_unjoined(run.sqanti, caller_support.by_native_id, /*novel_only=*/true, "the --caller-support matrix");
  }

  run.engine = c.config.empty() ? RuleEngine() : RuleEngine::from_toml(c.config);
  run.engine.apply_platform_preset(c.platform);
  if (c.platform == "ont") {
    std::fprintf(stderr, "warning: --platform ont is deprecated and will be removed in 0.2: since ruleset 0.0.7 the indel-near "
                         "gate reads the excess over each read's own indel rate and the preset has not changed a verdict on ONT "
                         "data; it still turns that gate off in this run\n");
  }

  AdjudicateInputs in;
  in.sqanti = &run.sqanti;
  in.chains = &chains;
  in.catalog = catalog_ptr;
  in.sj = sj_ptr;
  in.bam = bam.get();
  for (const auto& vb : validation_bams) in.validation_bams.push_back(vb.get());
  in.haplotype = haplo.get();
  in.haplotype_circular = haplotype_provenance_is_circular(c.haplotype_provenance);
  in.pangenome = pangenome_ptr;
  in.pangenome_circular = pangenome_provenance_is_circular(c.pangenome_provenance);
  in.caller_support = caller_support_ptr;
  run.results = adjudicate(in, run.engine, &run.junctions);
  if (haplo) run.haplotype_fasta_length_mismatches = haplo->n_length_mismatches();
  return run;
}

RunProvenance provenance_for(const CommonArgs& c, const LoadedRun& run) {
  RunProvenance prov;
  prov.tool_version = kVersion;
  prov.ruleset_version = run.engine.config().ruleset_version;
  prov.sqanti3_version_target = run.engine.config().sqanti3_version_target;
  prov.classification_path = c.classification;
  prov.classification_computed = c.no_sqanti3;
  prov.isoforms_path = c.isoforms_bed.empty() ? c.isoforms_gtf : c.isoforms_bed;
  prov.ref_gtf_path = c.ref_gtf;
  prov.sj_tab_path = c.sj_tab;
  prov.bam_path = c.bam_path;
  prov.validation_bam_paths = c.validation_bams;
  prov.caller_support_path = c.caller_support;
  prov.reference_path = c.reference;
  prov.reference_haplotype_paths = c.reference_haplotypes;
  prov.haplotype_vcf_paths = c.haplotype_vcfs;
  prov.haplotype_vcf_sample = c.vcf_sample;
  for (const auto& v : run.haplotype_vcfs) {
    prov.haplotype_vcf_unphased_het += v->n_unphased_het_applied();
    prov.haplotype_vcf_indels += v->n_indels_applied();
    prov.haplotype_vcf_deletion_match += v->n_deletions_equal_to_junctions();
    prov.haplotype_vcf_unapplied += v->n_unapplied_alleles();
    prov.haplotype_vcf_missing += v->n_missing_alleles();
    prov.haplotype_vcf_conflicting += v->n_conflicting_records();
    prov.haplotype_vcf_ref_mismatch += v->n_ref_mismatches();
  }
  prov.haplotype_fasta_length_mismatch = run.haplotype_fasta_length_mismatches;
  prov.pangenome_path = c.pangenome_junctions;
  prov.variant_axis_on = !c.reference_haplotypes.empty() || !c.haplotype_vcfs.empty();
  prov.variant_circular = haplotype_provenance_is_circular(c.haplotype_provenance);
  prov.pangenome_axis_on = !c.pangenome_junctions.empty();
  prov.pangenome_circular = pangenome_provenance_is_circular(c.pangenome_provenance);
  prov.consensus_axis_on = !c.caller_support.empty();
  prov.config_path = c.config;
  prov.thresholds = run.engine.thresholds_summary();
  return prov;
}

}  // namespace panisoguard
