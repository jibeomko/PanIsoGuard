// `panisoguard combine` -- integrate multiple callers' isoform sets by intron-
// chain fingerprint, preserving each caller's native running IDs, and emit a
// caller-support matrix (optionally annotated known/novel vs a reference GTF).

#include <cstdio>
#include <string>
#include <vector>

#include "panisoguard/bed12.hpp"
#include "panisoguard/consensus.hpp"
#include "panisoguard/gtf.hpp"

namespace {

struct LabeledInput {
  std::string caller;
  std::string path;
};

// "label:path" -> {label, path}; if no ':' present, derive label from filename.
LabeledInput parse_labeled(const std::string& spec) {
  const std::size_t colon = spec.find(':');
  if (colon != std::string::npos) {
    return {spec.substr(0, colon), spec.substr(colon + 1)};
  }
  std::size_t slash = spec.find_last_of('/');
  std::string base = (slash == std::string::npos) ? spec : spec.substr(slash + 1);
  std::size_t dot = base.find('.');
  if (dot != std::string::npos) base = base.substr(0, dot);
  return {base, spec};
}

void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard combine [options] --out matrix.tsv\n"
      "\n"
      "  --gtf LABEL:PATH     caller isoform GTF (repeatable)\n"
      "  --bed LABEL:PATH     caller isoform BED12 (repeatable)\n"
      "  --ref-gtf PATH       reference GTF for known/novel annotation (optional)\n"
      "  --out PATH           output caller-support matrix TSV (required)\n"
      "\n"
      "LABEL identifies the caller (e.g. flair, isoquant, bambu); if omitted the\n"
      "filename stem is used. Isoforms are integrated by intron-chain fingerprint.\n");
}

}  // namespace

int cmd_combine(int argc, char** argv) {
  std::vector<LabeledInput> gtfs;
  std::vector<LabeledInput> beds;
  std::string ref_gtf;
  std::string out_path;

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "panisoguard combine: %s requires an argument\n", name);
        usage();
        std::exit(1);
      }
      return argv[++i];
    };
    if (a == "--gtf")            gtfs.push_back(parse_labeled(next("--gtf")));
    else if (a == "--bed")       beds.push_back(parse_labeled(next("--bed")));
    else if (a == "--ref-gtf")   ref_gtf = next("--ref-gtf");
    else if (a == "--out")       out_path = next("--out");
    else if (a == "-h" || a == "--help") { usage(); return 0; }
    else {
      std::fprintf(stderr, "panisoguard combine: unknown option '%s'\n", a.c_str());
      usage();
      return 1;
    }
  }

  if (gtfs.empty() && beds.empty()) {
    std::fprintf(stderr, "panisoguard combine: need at least one --gtf or --bed input\n");
    usage();
    return 1;
  }
  if (out_path.empty()) {
    std::fprintf(stderr, "panisoguard combine: --out is required\n");
    usage();
    return 1;
  }

  using namespace panisoguard;
  ConsensusBuilder builder;

  for (const auto& in : gtfs) {
    const auto tx = read_gtf_transcripts(in.path);
    for (const auto& t : tx) builder.add(in.caller, t.id, t.chain);
    std::fprintf(stderr, "  + gtf  %-12s %8zu isoforms  (%s)\n", in.caller.c_str(), tx.size(),
                 in.path.c_str());
  }
  for (const auto& in : beds) {
    const auto recs = read_bed12(in.path);
    for (const auto& r : recs) builder.add(in.caller, r.name, r.chain);
    std::fprintf(stderr, "  + bed  %-12s %8zu isoforms  (%s)\n", in.caller.c_str(), recs.size(),
                 in.path.c_str());
  }

  Catalog catalog;
  const Catalog* catalog_ptr = nullptr;
  if (!ref_gtf.empty()) {
    std::fprintf(stderr, "  reading reference catalog: %s\n", ref_gtf.c_str());
    catalog = build_catalog_from_gtf(ref_gtf);
    catalog_ptr = &catalog;
  }

  const auto isoforms = builder.build(catalog_ptr);
  write_caller_support_matrix(out_path, isoforms);

  // Summary to stderr.
  std::size_t multi = 0, known = 0, novel = 0, mono = 0;
  for (const auto& iso : isoforms) {
    if (iso.n_callers() > 1) ++multi;
    if (iso.monoexonic) ++mono;
    else if (iso.catalog_checked) (iso.known_in_catalog ? known : novel)++;
  }
  std::fprintf(stderr,
      "integrated %zu isoforms  (multi-caller: %zu",
      isoforms.size(), multi);
  if (catalog_ptr) std::fprintf(stderr, ", known: %zu, novel: %zu", known, novel);
  if (mono) std::fprintf(stderr, ", monoexonic: %zu", mono);
  std::fprintf(stderr, ")\nwrote %s\n", out_path.c_str());
  return 0;
}
