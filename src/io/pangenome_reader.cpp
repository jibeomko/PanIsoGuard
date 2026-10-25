#include "panisoguard/pangenome.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {

void PangenomeJunctions::add(PangenomeJunctionRecord rec) {
  const std::string k = junction_key(rec.chrom, rec.strand, rec.intron);
  auto it = by_coord_.find(k);
  if (it != by_coord_.end()) {
    // Duplicate graph-junction rows should not create unreachable shadow records.
    // Treat n_haplotypes as an already-aggregated support count and preserve the
    // strongest support seen for this coordinate.
    PangenomeJunctionRecord& cur = records_[it->second];
    cur.n_haplotypes = std::max(cur.n_haplotypes, rec.n_haplotypes);
    for (const auto& id : rec.haplotype_ids) {
      if (std::find(cur.haplotype_ids.begin(), cur.haplotype_ids.end(), id) == cur.haplotype_ids.end())
        cur.haplotype_ids.push_back(id);
    }
    cur.n_haplotypes = std::max(cur.n_haplotypes, static_cast<int>(cur.haplotype_ids.size()));
    return;
  }
  records_.push_back(std::move(rec));
  by_coord_.emplace(k, records_.size() - 1);
}

const PangenomeJunctionRecord* PangenomeJunctions::find_exact(const std::string& chrom, Strand strand,
                                                              const Junction& intron) const {
  auto it = by_coord_.find(junction_key(chrom, strand, intron));
  if (it == by_coord_.end()) return nullptr;
  return &records_[it->second];
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
    const int64_t s1 = parse_int_field(f[1], "intron_start", "pangenome junctions", lineno);
    const int64_t e1 = parse_int_field(f[2], "intron_end", "pangenome junctions", lineno);
    if (s1 > e1) {
      throw std::runtime_error("pangenome junctions line " + std::to_string(lineno) +
                               ": intron_start (" + std::to_string(s1) + ") > intron_end (" +
                               std::to_string(e1) + ")");
    }
    rec.intron = Junction{s1 - 1, e1};
    rec.strand = parse_strand(f[3]);
    // A splice junction is strand-specific; reject a missing/typo'd strand instead of
    // silently storing never-matching kUnknown data.
    if (rec.strand == Strand::kUnknown) {
      throw std::runtime_error("pangenome junctions line " + std::to_string(lineno) +
                               ": unrecognized strand \"" + f[3] + "\" (expected +, -, 1, or 2)");
    }
    if (f.size() >= 5 && !f[4].empty()) {
      rec.n_haplotypes = static_cast<int>(parse_int_field(f[4], "n_haplotypes", "pangenome junctions", lineno));
      if (rec.n_haplotypes < 1) {
        throw std::runtime_error("pangenome junctions line " + std::to_string(lineno) +
                                 ": n_haplotypes must be >= 1, got " + f[4]);
      }
    }
    if (f.size() >= 6 && !f[5].empty() && f[5] != "NA" && f[5] != ".") {
      for (const std::string& id : split(f[5], ','))
        if (!id.empty()) rec.haplotype_ids.push_back(id);
      rec.n_haplotypes = std::max(rec.n_haplotypes, static_cast<int>(rec.haplotype_ids.size()));
    }
    table.add(std::move(rec));
  }
  return table;
}

}  // namespace panisoguard
