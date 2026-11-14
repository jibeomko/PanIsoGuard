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
  std::string gene_id;
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
    const Strand strand = parse_strand(f[6].empty() ? '.' : f[6][0]);
    auto it = acc.find(tid);
    if (it == acc.end()) {
      ExonAcc a;
      a.chrom = f[0];
      a.gene_id = extract_gtf_attr(f[8], "gene_id");
      a.strand = strand;
      it = acc.emplace(tid, std::move(a)).first;
      order.push_back(tid);
    } else if (it->second.chrom != f[0] || it->second.strand != strand) {
      throw std::runtime_error("GTF line " + std::to_string(lineno) + ": transcript_id \"" + tid +
                               "\" has exons on more than one chromosome or strand");
    }
    // Context-rich error on a malformed exon coordinate instead of a bare stoll terminate.
    const int64_t es = parse_int_field(f[3], "exon start", "GTF", lineno);
    const int64_t ee = parse_int_field(f[4], "exon end", "GTF", lineno);
    it->second.exons.emplace_back(es, ee);
  }

  std::vector<Transcript> out;
  out.reserve(order.size());
  for (const auto& tid : order) {
    Transcript t;
    t.id = tid;
    t.chain = chain_from_exons(acc[tid]);
    t.gene_id = acc[tid].gene_id;
    for (const auto& e : acc[tid].exons) {
      if (t.end == 0 || e.first - 1 < t.start) t.start = e.first - 1;   // 1-based inclusive exon -> 0-based half-open span
      if (e.second > t.end) t.end = e.second;
      t.exons.emplace_back(e.first - 1, e.second);
    }
    std::sort(t.exons.begin(), t.exons.end());
    out.push_back(std::move(t));
  }
  return out;
}

namespace {
std::string locus_key(const std::string& chrom, Strand strand) { return chrom + '|' + strand_char(strand); }
std::string site_key(const std::string& chrom, Strand strand, int64_t pos) { return locus_key(chrom, strand) + ':' + std::to_string(pos); }
}  // namespace

void Catalog::add_transcript(const Transcript& t) {
  txs_.push_back(t);
  const std::size_t idx = txs_.size() - 1;
  if (t.end > t.start) txs_by_locus_[locus_key(t.chain.chrom, t.chain.strand)].push_back(idx);
  if (t.chain.monoexonic()) return;
  chains_.insert(fingerprint_intron_chain(t.chain));
  for (const auto& j : t.chain.introns) {
    introns_.insert(junction_key(t.chain.chrom, t.chain.strand, j));
    intron_starts_.insert(site_key(t.chain.chrom, t.chain.strand, j.start));
    intron_ends_.insert(site_key(t.chain.chrom, t.chain.strand, j.end));
    txs_by_intron_[junction_key(t.chain.chrom, t.chain.strand, j)].push_back(idx);
  }
}

bool Catalog::has_intron_start(const std::string& chrom, Strand strand, int64_t pos) const {
  return intron_starts_.count(site_key(chrom, strand, pos)) != 0;
}
bool Catalog::has_intron_end(const std::string& chrom, Strand strand, int64_t pos) const {
  return intron_ends_.count(site_key(chrom, strand, pos)) != 0;
}

bool Catalog::has_subchain(const IntronChain& chain, int64_t start, int64_t end) const {
  if (chain.monoexonic()) return false;
  const auto it = txs_by_intron_.find(junction_key(chain.chrom, chain.strand, chain.introns.front()));
  if (it == txs_by_intron_.end()) return false;
  for (const std::size_t i : it->second) {
    const std::vector<Junction>& full = txs_[i].chain.introns;
    if (!contiguous_subchain(chain.introns, full)) continue;
    // the matched run is the reference introns equal to the chain; every other reference
    // intron must lie outside the isoform's span (otherwise the isoform retains it)
    bool retained = false;
    for (const auto& j : full) {
      bool in_chain = false;
      for (const auto& q : chain.introns) if (q == j) { in_chain = true; break; }
      if (!in_chain && j.start >= start && j.end <= end) { retained = true; break; }
    }
    if (!retained) return true;
  }
  return false;
}

std::string Catalog::mono_exon_category(const std::string& chrom, Strand strand, int64_t start, int64_t end) const {
  const auto it = txs_by_locus_.find(locus_key(chrom, strand));
  if (it == txs_by_locus_.end()) return "";
  bool in_exon_mono = false, in_exon_multi = false, touches_exon = false, in_span = false;
  for (const std::size_t i : it->second) {
    const Transcript& t = txs_[i];
    if (!(t.start < end && t.end > start)) continue;
    in_span = true;
    for (const auto& e : t.exons) {
      if (e.first <= start && end <= e.second) { (t.chain.monoexonic() ? in_exon_mono : in_exon_multi) = true; }
      else if (e.first < end && e.second > start) touches_exon = true;
    }
  }
  if (in_exon_mono) return "full-splice_match";
  if (in_exon_multi) return "incomplete-splice_match";
  if (touches_exon) return "genic";
  if (in_span) return "genic_intron";
  return "";
}

bool Catalog::overlaps(const std::string& chrom, Strand strand, int64_t start, int64_t end) const {
  const auto it = txs_by_locus_.find(locus_key(chrom, strand));
  if (it == txs_by_locus_.end()) return false;
  for (const std::size_t i : it->second)
    if (txs_[i].start < end && txs_[i].end > start) return true;
  return false;
}

bool Catalog::overlaps_any_strand(const std::string& chrom, int64_t start, int64_t end) const {
  return overlaps(chrom, Strand::kPlus, start, end) || overlaps(chrom, Strand::kMinus, start, end) ||
         overlaps(chrom, Strand::kUnknown, start, end);
}

Catalog build_catalog_from_gtf(const std::string& path) {
  Catalog cat;
  for (const auto& t : read_gtf_transcripts(path)) {
    cat.add_transcript(t);
  }
  return cat;
}

}  // namespace panisoguard
