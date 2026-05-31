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
    const int64_t s1 = std::stoll(f[1]);
    const int64_t e1 = std::stoll(f[2]);
    rec.intron = Junction{s1 - 1, e1};
    const int strand_code = std::stoi(f[3]);
    rec.strand = strand_code == 1 ? Strand::kPlus
               : strand_code == 2 ? Strand::kMinus
                                  : Strand::kUnknown;
    rec.motif = std::stoi(f[4]);
    rec.annotated = std::stoi(f[5]) != 0;
    rec.n_uniq = std::stoi(f[6]);
    rec.n_multi = std::stoi(f[7]);
    rec.max_overhang = std::stoi(f[8]);
    table.add(std::move(rec));
  }
  return table;
}

}  // namespace panisoguard
