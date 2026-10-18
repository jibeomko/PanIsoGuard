// `panisoguard decoys` -- how often does whole-chain read support fire on a chain the
// annotation does not know? Recombine the annotated transcripts of each gene into
// novel-in-catalog decoy chains, count the reads that carry each decoy end to end in the
// BAM (the same count the chain axis uses), and report the fraction at or above the lift
// threshold next to the same fraction for the annotated parents.

#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "panisoguard/bam_features.hpp"
#include "panisoguard/decoys.hpp"
#include "panisoguard/gtf.hpp"
#include "panisoguard/version.hpp"

namespace {

void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard decoys --ref-gtf PATH --bam PATH --out PREFIX [options]\n"
      "\n"
      "  --ref-gtf PATH          reference annotation; the decoys recombine its transcripts\n"
      "  --bam PATH              indexed BAM/CRAM of the long reads (the same file adjudicate used)\n"
      "  --reference PATH        genome FASTA (only to decode a CRAM)\n"
      "  --min-full-reads N      lift threshold to evaluate (default 3 = chain_min_full_reads)\n"
      "  --max-per-gene N        at most N decoys per gene, sampled (default 10; 0 = all)\n"
      "  --seed N                sampling seed (default 1)\n"
      "  --out PREFIX            writes PREFIX.decoys.tsv (one row per decoy) and PREFIX.decoys.json\n"
      "\n"
      "A decoy is the 5' introns of one annotated transcript joined to the 3' introns of another\n"
      "transcript of the same gene at a shared intron, when the result is not an annotated chain\n"
      "or a piece of one. The fraction of decoys with >= N whole-chain reads is an upper bound\n"
      "on the chain axis's false-lift rate in this dataset (some decoys are real, unannotated).\n");
}

struct Tally {
  std::size_t n = 0, lifted = 0;
  void add(bool l) { ++n; lifted += l; }
  double rate() const { return n ? double(lifted) / double(n) : 0.0; }
};

std::string introns_1based(const panisoguard::IntronChain& c) {
  std::string s;
  for (const auto& j : c.introns) {
    if (!s.empty()) s += ',';
    s += std::to_string(j.start + 1) + '-' + std::to_string(j.end);
  }
  return s;
}

}  // namespace

int cmd_decoys(int argc, char** argv) {
  using namespace panisoguard;
  std::string ref_gtf, bam_path, reference, out_prefix;
  int min_full = 3, max_per_gene = 10;
  unsigned seed = 1;
  auto next = [&](int& i, const char* name) -> std::string {
    if (i + 1 >= argc) { std::fprintf(stderr, "decoys: %s requires an argument\n", name); std::exit(1); }
    return argv[++i];
  };
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--ref-gtf") ref_gtf = next(i, "--ref-gtf");
    else if (a == "--bam") bam_path = next(i, "--bam");
    else if (a == "--reference") reference = next(i, "--reference");
    else if (a == "--out") out_prefix = next(i, "--out");
    else if (a == "--min-full-reads") min_full = std::stoi(next(i, "--min-full-reads"));
    else if (a == "--max-per-gene") max_per_gene = std::stoi(next(i, "--max-per-gene"));
    else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(next(i, "--seed")));
    else if (a == "-h" || a == "--help") { usage(); return 0; }
    else { std::fprintf(stderr, "decoys: unknown option '%s'\n", a.c_str()); usage(); return 1; }
  }
  if (ref_gtf.empty() || bam_path.empty() || out_prefix.empty()) {
    std::fprintf(stderr, "decoys: --ref-gtf, --bam and --out are required\n");
    usage();
    return 1;
  }
  try {
    const std::vector<Transcript> ts = read_gtf_transcripts(ref_gtf);
    std::fprintf(stderr, "read %zu reference transcripts\n", ts.size());
    const std::vector<Decoy> decoys = make_decoys(ts, max_per_gene, seed);
    std::fprintf(stderr, "built %zu decoy chains\n", decoys.size());

    BamReader bam(bam_path, reference);
    // parents: the whole-chain read count of every annotated transcript a decoy was built from
    std::unordered_map<std::string, const Transcript*> by_id;
    for (const auto& t : ts) by_id.emplace(t.id, &t);
    std::unordered_map<std::string, int> parent_full;
    Tally parents;
    for (const auto& d : decoys) {
      for (const std::string* p : {&d.parent_a, &d.parent_b}) {
        if (parent_full.count(*p)) continue;
        const int n = bam.chain_support(by_id.at(*p)->chain).n_full;
        parent_full[*p] = n;
        parents.add(n >= min_full);
      }
    }
    std::ofstream tsv(out_prefix + ".decoys.tsv");
    if (!tsv) throw std::runtime_error("cannot write " + out_prefix + ".decoys.tsv");
    tsv << "decoy_id\tchrom\tstrand\tgene_id\tn_introns\tparent_a\tparent_b\tintrons_1based\t"
           "n_full_reads\tn_partial_reads\tparent_a_full_reads\tparent_b_full_reads\tlifted\n";
    Tally all, both_parents;
    std::map<std::size_t, Tally> by_n;
    std::size_t k = 0;
    for (const auto& d : decoys) {
      const ChainBamSupport s = bam.chain_support(d.chain);
      const bool lifted = s.evaluated && s.n_full >= min_full;
      const int pa = parent_full[d.parent_a], pb = parent_full[d.parent_b];
      all.add(lifted);
      by_n[d.chain.introns.size()].add(lifted);
      if (pa >= min_full && pb >= min_full) both_parents.add(lifted);
      tsv << "decoy_" << ++k << '\t' << d.chain.chrom << '\t' << strand_char(d.chain.strand) << '\t'
          << d.gene_id << '\t' << d.chain.introns.size() << '\t' << d.parent_a << '\t' << d.parent_b << '\t'
          << introns_1based(d.chain) << '\t' << s.n_full << '\t' << s.n_partial << '\t' << pa << '\t' << pb
          << '\t' << (lifted ? "true" : "false") << '\n';
      if (k % 20000 == 0) std::fprintf(stderr, "  %zu / %zu decoys\n", k, decoys.size());
    }
    std::ofstream js(out_prefix + ".decoys.json");
    if (!js) throw std::runtime_error("cannot write " + out_prefix + ".decoys.json");
    auto tally = [&](const char* name, const Tally& t, bool last = false) {
      js << "  \"" << name << "\": {\"n\": " << t.n << ", \"lifted\": " << t.lifted
         << ", \"rate\": " << t.rate() << "}" << (last ? "\n" : ",\n");
    };
    js << "{\n  \"tool_version\": \"" << kVersion << "\",\n  \"ref_gtf\": \"" << ref_gtf << "\",\n  \"bam\": \""
       << bam_path << "\",\n  \"min_full_reads\": " << min_full << ",\n  \"max_per_gene\": " << max_per_gene
       << ",\n  \"seed\": " << seed << ",\n";
    tally("decoys", all);
    tally("decoys_both_parents_lifted", both_parents);
    tally("parents", parents);
    js << "  \"decoys_by_n_introns\": {";
    bool first = true;
    for (const auto& kv : by_n) {
      js << (first ? "" : ", ") << "\"" << kv.first << "\": {\"n\": " << kv.second.n << ", \"lifted\": " << kv.second.lifted << "}";
      first = false;
    }
    js << "}\n}\n";
    std::fprintf(stderr,
        "decoys with >= %d whole-chain reads: %zu / %zu (%.4f)\n"
        "  among decoys whose two parents both have >= %d: %zu / %zu (%.4f)\n"
        "  annotated parents with >= %d: %zu / %zu (%.4f)\n"
        "wrote %s.{decoys.tsv,decoys.json}\n",
        min_full, all.lifted, all.n, all.rate(), min_full, both_parents.lifted, both_parents.n, both_parents.rate(),
        min_full, parents.lifted, parents.n, parents.rate(), out_prefix.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "decoys: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
