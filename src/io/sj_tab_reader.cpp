#include "panisoguard/sj_tab.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {

void SjTable::add(SjRecord rec) {
  const std::string k = junction_key(rec.chrom, rec.strand, rec.intron);
  auto it = by_coord_.find(k);
  if (it != by_coord_.end()) {
    SjRecord& cur = records_[it->second];
    // Merged SJ.tab files can contain the same strand-aware junction more than once.
    // Keep one lookup record and aggregate read support instead of letting the first
    // row silently win.
    cur.annotated = cur.annotated || rec.annotated;
    cur.n_uniq += rec.n_uniq;
    cur.n_multi += rec.n_multi;
    cur.max_overhang = std::max(cur.max_overhang, rec.max_overhang);
    if (!cur.canonical() && rec.canonical()) cur.motif = rec.motif;
    return;
  }
  records_.push_back(std::move(rec));
  by_coord_.emplace(k, records_.size() - 1);
  span_indexed_ = false;
}

bool SjTable::any_in_span(const std::string& chrom, Strand strand, int64_t start, int64_t end) const {
  if (!span_indexed_) {
    by_chrom_.clear();
    longest_.clear();
    for (const SjRecord& r : records_) {
      by_chrom_[r.chrom].push_back({r.intron.start, r.intron.end, r.strand});
      longest_[r.chrom] = std::max(longest_[r.chrom], r.intron.end - r.intron.start);
    }
    for (auto& kv : by_chrom_)
      std::sort(kv.second.begin(), kv.second.end(), [](const Span& a, const Span& b) { return a.start < b.start; });
    span_indexed_ = true;
  }
  const auto it = by_chrom_.find(chrom);
  if (it == by_chrom_.end()) return false;
  const std::vector<Span>& v = it->second;
  // A junction overlapping [start, end) starts after start - longest intron and before end.
  auto s = std::lower_bound(v.begin(), v.end(), start - longest_.at(chrom),
                            [](const Span& e, int64_t x) { return e.start < x; });
  for (; s != v.end() && s->start < end; ++s) {
    if (s->end > start && (s->strand == Strand::kUnknown || strand == Strand::kUnknown || s->strand == strand)) return true;
  }
  return false;
}

const SjRecord* SjTable::find_exact(const std::string& chrom, Strand strand,
                                    const Junction& intron) const {
  auto it = by_coord_.find(junction_key(chrom, strand, intron));
  if (it == by_coord_.end()) return nullptr;
  return &records_[it->second];
}

SjTable read_sj_tab(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open SJ.tab: " + path);

  SjTable table;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    chomp(line);
    if (line.empty()) continue;
    std::vector<std::string> f = split_tsv(line);
    if (f.size() < 9) {
      throw std::runtime_error("SJ.tab line " + std::to_string(lineno) +
                               ": expected 9 columns, got " + std::to_string(f.size()));
    }
    SjRecord rec;
    rec.chrom = f[0];
    // 1-based inclusive intron [s, e] -> 0-based half-open [s-1, e].
    const int64_t s1 = parse_int_field(f[1], "intron_start", "SJ.tab", lineno);
    const int64_t e1 = parse_int_field(f[2], "intron_end", "SJ.tab", lineno);
    if (s1 > e1) {
      throw std::runtime_error("SJ.tab line " + std::to_string(lineno) +
                               ": intron_start (" + std::to_string(s1) + ") > intron_end (" +
                               std::to_string(e1) + ")");
    }
    rec.intron = Junction{s1 - 1, e1};
    // STAR strand: 0 = undefined (kept as Unknown), 1 = +, 2 = -.
    const int strand_code = static_cast<int>(parse_int_field(f[3], "strand", "SJ.tab", lineno));
    rec.strand = strand_code == 1 ? Strand::kPlus
               : strand_code == 2 ? Strand::kMinus
                                  : Strand::kUnknown;
    rec.motif = static_cast<int>(parse_int_field(f[4], "motif", "SJ.tab", lineno));
    rec.annotated = parse_int_field(f[5], "annotated", "SJ.tab", lineno) != 0;
    rec.n_uniq = static_cast<int>(parse_int_field(f[6], "n_uniq", "SJ.tab", lineno));
    rec.n_multi = static_cast<int>(parse_int_field(f[7], "n_multi", "SJ.tab", lineno));
    rec.max_overhang = static_cast<int>(parse_int_field(f[8], "max_overhang", "SJ.tab", lineno));
    table.add(std::move(rec));
  }
  return table;
}

}  // namespace panisoguard
