#include "panisoguard/consensus.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace panisoguard {

void ConsensusBuilder::add(const std::string& caller, const std::string& native_id,
                           const IntronChain& chain) {
  std::string key;
  if (chain.monoexonic()) {
    // Not chain-fingerprintable: keep each monoexonic isoform in its own group so
    // they never silently collapse onto each other.
    key = "MONO:" + std::to_string(mono_counter_++);
    Group g;
    g.monoexonic = true;
    g.chain = chain;
    g.support[caller].push_back(native_id);
    groups_.emplace(std::move(key), std::move(g));
    return;
  }

  const Fingerprint fp = fingerprint_intron_chain(chain);
  key = "FP:" + std::to_string(fp);
  auto it = groups_.find(key);
  if (it == groups_.end()) {
    Group g;
    g.fp = fp;
    g.monoexonic = false;
    g.chain = chain;
    it = groups_.emplace(std::move(key), std::move(g)).first;
  }
  it->second.support[caller].push_back(native_id);
}

std::vector<ConsensusIsoform> ConsensusBuilder::build(const Catalog* catalog) const {
  std::vector<ConsensusIsoform> out;
  out.reserve(groups_.size());

  for (const auto& kv : groups_) {
    const Group& g = kv.second;
    ConsensusIsoform iso;
    iso.chrom = g.chain.chrom;
    iso.strand = g.chain.strand;
    iso.chain = g.chain;
    iso.monoexonic = g.monoexonic;
    iso.support = g.support;
    if (catalog != nullptr && !g.monoexonic) {
      iso.catalog_checked = true;
      iso.known_in_catalog = catalog->has_chain(g.fp);
    }
    out.push_back(std::move(iso));
  }

  // Deterministic genomic ordering so PIG ids are stable across runs.
  auto first_start = [](const ConsensusIsoform& x) -> int64_t {
    return x.chain.introns.empty() ? -1 : x.chain.introns.front().start;
  };
  std::sort(out.begin(), out.end(), [&](const ConsensusIsoform& a, const ConsensusIsoform& b) {
    if (a.chrom != b.chrom) return a.chrom < b.chrom;
    if (a.strand != b.strand) return strand_char(a.strand) < strand_char(b.strand);
    const int64_t sa = first_start(a), sb = first_start(b);
    if (sa != sb) return sa < sb;
    // tie-break on full intron vector for total order
    return a.chain.introns < b.chain.introns;
  });

  char buf[32];
  for (std::size_t i = 0; i < out.size(); ++i) {
    std::snprintf(buf, sizeof(buf), "PIG.%06zu", i + 1);
    out[i].pig_id = buf;
  }
  return out;
}

namespace {

std::string join_native_ids(const std::map<std::string, std::vector<std::string>>& support) {
  std::string s;
  for (const auto& kv : support) {
    if (!s.empty()) s += ';';
    s += kv.first;
    s += '=';
    for (std::size_t i = 0; i < kv.second.size(); ++i) {
      if (i) s += '|';
      s += kv.second[i];
    }
  }
  return s;
}

std::string join_callers(const std::map<std::string, std::vector<std::string>>& support) {
  std::string s;
  for (const auto& kv : support) {
    if (!s.empty()) s += ',';
    s += kv.first;
  }
  return s;
}

}  // namespace

void write_caller_support_matrix(const std::string& path,
                                 const std::vector<ConsensusIsoform>& isoforms) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write caller support matrix: " + path);

  out << "pig_id\tchrom\tstrand\tn_introns\tn_callers\tn_isoforms\tcallers\tnovelty\tnative_ids\n";
  for (const auto& iso : isoforms) {
    std::string novelty = "NA";
    if (iso.monoexonic) {
      novelty = "monoexonic";
    } else if (iso.catalog_checked) {
      novelty = iso.known_in_catalog ? "known" : "novel";
    }
    out << iso.pig_id << '\t' << iso.chrom << '\t' << strand_char(iso.strand) << '\t'
        << iso.chain.introns.size() << '\t' << iso.n_callers() << '\t' << iso.n_isoforms()
        << '\t' << join_callers(iso.support) << '\t' << novelty << '\t'
        << join_native_ids(iso.support) << '\n';
  }
}

}  // namespace panisoguard
