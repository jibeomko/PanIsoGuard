#include "panisoguard/sj_tab.hpp"

#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {

std::string SjTable::key(const std::string& chrom, Strand strand, const Junction& j) {
  std::string k;
  k.reserve(chrom.size() + 26);
  k += chrom;
  k += '|';
  k += strand_char(strand);
  k += ':';
  k += std::to_string(j.start);
  k += '-';
  k += std::to_string(j.end);
  return k;
}

void SjTable::add(SjRecord rec) {
  const std::string k = key(rec.chrom, rec.strand, rec.intron);
  records_.push_back(std::move(rec));
  by_coord_.emplace(k, records_.size() - 1);
}

const SjRecord* SjTable::find_exact(const std::string& chrom, Strand strand,
                                    const Junction& intron) const {
  auto it = by_coord_.find(key(chrom, strand, intron));
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
