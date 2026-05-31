#include "panisoguard/bed12.hpp"

#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {
namespace {

// Parse a comma-separated list of integers (BED blockSizes / blockStarts;
// a trailing comma is allowed).
std::vector<int64_t> parse_int_list(const std::string& s) {
  std::vector<int64_t> out;
  std::size_t start = 0;
  while (start < s.size()) {
    std::size_t comma = s.find(',', start);
    const std::string tok = s.substr(start, (comma == std::string::npos ? s.size() : comma) - start);
    if (!tok.empty()) out.push_back(std::stoll(tok));
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  return out;
}

}  // namespace

std::vector<Bed12Record> read_bed12(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open BED12: " + path);

  std::vector<Bed12Record> out;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    chomp(line);
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
    std::vector<std::string> f = split_tsv(line);
    if (f.size() < 12) {
      throw std::runtime_error("BED12 line " + std::to_string(lineno) +
                               ": expected 12 columns, got " + std::to_string(f.size()));
    }

    Bed12Record rec;
    const std::string chrom = f[0];
    const int64_t chrom_start = std::stoll(f[1]);  // 0-based
    rec.name = f[3];
    rec.id_prefix = flair_id_prefix(rec.name);
    rec.score = f[4] == "." ? 0 : std::stoi(f[4]);

    rec.chain.chrom = chrom;
    rec.chain.strand = parse_strand(f[5].empty() ? '.' : f[5][0]);

    const int block_count = std::stoi(f[9]);
    const std::vector<int64_t> sizes = parse_int_list(f[10]);
    const std::vector<int64_t> starts = parse_int_list(f[11]);
    if (static_cast<int>(sizes.size()) != block_count ||
        static_cast<int>(starts.size()) != block_count) {
      throw std::runtime_error("BED12 line " + std::to_string(lineno) +
                               ": blockCount does not match blockSizes/blockStarts");
    }

    // Exons are 0-based half-open [chrom_start+start, chrom_start+start+size).
    // Intron between consecutive exons = [exon_i.end, exon_{i+1}.start).
    for (int i = 0; i + 1 < block_count; ++i) {
      const int64_t exon_i_end = chrom_start + starts[i] + sizes[i];
      const int64_t exon_next_start = chrom_start + starts[i + 1];
      if (exon_next_start > exon_i_end) {
        rec.chain.introns.push_back(Junction{exon_i_end, exon_next_start});
      }
    }
    rec.chain.sort_introns();
    out.push_back(std::move(rec));
  }
  return out;
}

}  // namespace panisoguard
