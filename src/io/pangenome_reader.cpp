#include "panisoguard/pangenome.hpp"

#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {

std::string PangenomeJunctions::key(const std::string& chrom, Strand strand, const Junction& j) {
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

void PangenomeJunctions::add(PangenomeJunctionRecord rec) {
  const std::string k = key(rec.chrom, rec.strand, rec.intron);
  records_.push_back(std::move(rec));
  by_coord_.emplace(k, records_.size() - 1);
}

const PangenomeJunctionRecord* PangenomeJunctions::find_exact(const std::string& chrom, Strand strand,
                                                              const Junction& intron) const {
  auto it = by_coord_.find(key(chrom, strand, intron));
  if (it == by_coord_.end()) return nullptr;
  return &records_[it->second];
}

bool PangenomeJunctions::supports(const std::string& chrom, Strand strand, const Junction& intron,
                                  int min_haplotypes) const {
  const PangenomeJunctionRecord* r = find_exact(chrom, strand, intron);
  return r != nullptr && r->n_haplotypes >= min_haplotypes;
}

namespace {

Strand parse_strand(const std::string& s) {
  if (s == "+" || s == "1") return Strand::kPlus;
  if (s == "-" || s == "2") return Strand::kMinus;
  return Strand::kUnknown;
}

}  // namespace

PangenomeJunctions read_pangenome_junctions(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open pangenome junctions: " + path);

  PangenomeJunctions table;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    chomp(line);
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f = split_tsv(line);
    if (f.size() < 4) {
      throw std::runtime_error("pangenome junctions line " + std::to_string(lineno) +
                               ": expected >= 4 columns (chrom start end strand), got " +
                               std::to_string(f.size()));
    }
    PangenomeJunctionRecord rec;
    rec.chrom = f[0];
    // 1-based inclusive intron [s, e] -> 0-based half-open [s-1, e].
    const int64_t s1 = std::stoll(f[1]);
    const int64_t e1 = std::stoll(f[2]);
    rec.intron = Junction{s1 - 1, e1};
    rec.strand = parse_strand(f[3]);
    if (f.size() >= 5 && !f[4].empty()) {
      rec.n_haplotypes = std::stoi(f[4]);
      if (rec.n_haplotypes < 1) {
        throw std::runtime_error("pangenome junctions line " + std::to_string(lineno) +
                                 ": n_haplotypes must be >= 1, got " + f[4]);
      }
    }
    table.add(std::move(rec));
  }
  return table;
}

}  // namespace panisoguard
