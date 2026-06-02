#include "panisoguard/gtf.hpp"

#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "panisoguard/tsv.hpp"

namespace panisoguard {

std::string extract_gtf_attr(const std::string& attrs, const std::string& key) {
  // Find `key` as a token (preceded by start/space/;), then take the following
  // value, which may be quoted ("...") or bare up to ';'.
  std::size_t pos = 0;
  while ((pos = attrs.find(key, pos)) != std::string::npos) {
    const bool boundary_ok = (pos == 0) || attrs[pos - 1] == ' ' || attrs[pos - 1] == ';';
    std::size_t after = pos + key.size();
    if (boundary_ok && after < attrs.size() && (attrs[after] == ' ' || attrs[after] == '"' || attrs[after] == '=')) {
      // advance past spaces and an optional '=' to the value
      while (after < attrs.size() && (attrs[after] == ' ' || attrs[after] == '=')) ++after;
      if (after < attrs.size() && attrs[after] == '"') {
        const std::size_t b = after + 1;
        const std::size_t e = attrs.find('"', b);
        if (e == std::string::npos) return "";
        return attrs.substr(b, e - b);
      }
      const std::size_t e = attrs.find_first_of("; ", after);
      return attrs.substr(after, (e == std::string::npos ? attrs.size() : e) - after);
    }
    pos = after;
  }
  return "";
}

namespace {

struct ExonAcc {
  std::string chrom;
  Strand strand = Strand::kUnknown;
  std::vector<std::pair<int64_t, int64_t>> exons;  // 1-based inclusive [s, e]
};

IntronChain chain_from_exons(const ExonAcc& acc) {
  IntronChain chain;
  chain.chrom = acc.chrom;
  chain.strand = acc.strand;
  std::vector<std::pair<int64_t, int64_t>> exons = acc.exons;
  std::sort(exons.begin(), exons.end());
  for (std::size_t i = 0; i + 1 < exons.size(); ++i) {
    // GTF exons are 1-based inclusive [s, e]. In 0-based half-open coords the intron
    // between exon_i and exon_{i+1} is [e_i, s_{i+1}-1): the 1-based-inclusive exon end
    // e_i is the first intronic base (0-based), and s_{i+1}-1 is one past the last.
    const int64_t intron_start = exons[i].second;        // = e_i
    const int64_t intron_end = exons[i + 1].first - 1;   // = s_{i+1} - 1
    if (intron_end > intron_start) {
      chain.introns.push_back(Junction{intron_start, intron_end});
    }
  }
  chain.sort_introns();
  return chain;
}

}  // namespace

std::vector<Transcript> read_gtf_transcripts(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open GTF: " + path);

  // Preserve first-seen order of transcripts for deterministic output.
  std::unordered_map<std::string, ExonAcc> acc;
  std::vector<std::string> order;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    chomp(line);
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f = split_tsv(line);
    if (f.size() < 9) continue;
    if (f[2] != "exon") continue;
    const std::string tid = extract_gtf_attr(f[8], "transcript_id");
    if (tid.empty()) continue;
    auto it = acc.find(tid);
    if (it == acc.end()) {
      ExonAcc a;
      a.chrom = f[0];
      a.strand = parse_strand(f[6].empty() ? '.' : f[6][0]);
      it = acc.emplace(tid, std::move(a)).first;
      order.push_back(tid);
    }
    // Context-rich error on a malformed exon coordinate instead of a bare stoll terminate.
    const int64_t es = parse_int_field(f[3], "exon start", "GTF", lineno);
    const int64_t ee = parse_int_field(f[4], "exon end", "GTF", lineno);
    it->second.exons.emplace_back(es, ee);
  }

  std::vector<Transcript> out;
  out.reserve(order.size());
  for (const auto& tid : order) {
    out.push_back(Transcript{tid, chain_from_exons(acc[tid])});
  }
  return out;
}

std::string Catalog::intron_key(const std::string& chrom, Strand strand, const Junction& j) {
  return chrom + '|' + strand_char(strand) + '|' + std::to_string(j.start) + '-' +
         std::to_string(j.end);
}

std::string Catalog::site_key(const std::string& chrom, Strand strand, int64_t pos) {
  return chrom + '|' + strand_char(strand) + '|' + std::to_string(pos);
}

void Catalog::add_transcript(const Transcript& t) {
  if (!t.chain.introns.empty()) {
    chains_.insert(fingerprint_intron_chain(t.chain));
  }
  for (const auto& j : t.chain.introns) {
    introns_.insert(intron_key(t.chain.chrom, t.chain.strand, j));
    sites_.insert(site_key(t.chain.chrom, t.chain.strand, j.start));
    sites_.insert(site_key(t.chain.chrom, t.chain.strand, j.end));
  }
}

Catalog build_catalog_from_gtf(const std::string& path) {
  Catalog cat;
  for (const auto& t : read_gtf_transcripts(path)) {
    cat.add_transcript(t);
  }
  return cat;
}

}  // namespace panisoguard
