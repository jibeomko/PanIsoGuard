#include "panisoguard/sqanti.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

#include "panisoguard/tsv.hpp"

namespace panisoguard {
namespace {

const std::string& field_at(const std::vector<std::string>& f, int idx) {
  static const std::string kEmpty;
  if (idx < 0 || idx >= static_cast<int>(f.size())) return kEmpty;
  return f[idx];
}

bool is_na_token(const std::string& s) {
  return s.empty() || s == "NA" || s == "NaN" || s == "nan" || s == ".";
}

constexpr int kMaxParseWarnings = 5;

// A non-NA token that does not fully parse as a number is genuinely malformed (not an
// intentional NA). Fall back to the sentinel so a bad QC value degrades to
// not_evaluable rather than a wrong prior, but warn (bounded) so it is not SILENT.
void warn_unparseable(const char* what, std::size_t lineno, const std::string& s, int& warnings) {
  if (warnings < kMaxParseWarnings) {
    std::fprintf(stderr,
                 "WARNING: SQANTI3 classification line %zu: %s value \"%s\" is not numeric; "
                 "treated as NA (not_evaluable)\n",
                 lineno, what, s.c_str());
    if (++warnings == kMaxParseWarnings) {
      std::fprintf(stderr, "WARNING: further SQANTI3 numeric-parse warnings suppressed\n");
    }
  }
}

double parse_double_or_na(const std::string& s, const char* what, std::size_t lineno, int& warnings) {
  if (is_na_token(s)) return SqantiRecord::kNaN;
  try {
    std::size_t pos = 0;
    const double v = std::stod(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("trailing characters");
    return v;
  } catch (...) {
    warn_unparseable(what, lineno, s, warnings);
    return SqantiRecord::kNaN;
  }
}

int parse_int_or_na(const std::string& s, const char* what, std::size_t lineno, int& warnings) {
  if (is_na_token(s)) return SqantiRecord::kIntNA;
  try {
    std::size_t pos = 0;
    const long v = std::stol(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("trailing characters");
    return static_cast<int>(v);
  } catch (...) {
    warn_unparseable(what, lineno, s, warnings);
    return SqantiRecord::kIntNA;
  }
}

}  // namespace

SqantiTable read_sqanti_classification(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open SQANTI3 classification: " + path);

  std::string header_line;
  if (!std::getline(in, header_line)) {
    throw std::runtime_error("empty SQANTI3 classification: " + path);
  }
  chomp(header_line);
  const TsvHeader header(split_tsv(header_line));

  if (!header.has("isoform")) {
    throw std::runtime_error("SQANTI3 classification missing mandatory 'isoform' column");
  }

  // Resolve column indices once. Absent expected columns are recorded so the
  // engine can mark the corresponding evidence not_evaluable rather than crash.
  SqantiTable table;
  table.n_columns = header.ncol();

  const int c_isoform = header.col("isoform");
  const int c_chrom = header.col("chrom");
  const int c_strand = header.col("strand");
  const int c_cat = header.col("structural_category");
  const int c_subcat = header.col("subcategory");
  const int c_gene = header.col("associated_gene");
  const int c_tx = header.col("associated_transcript");
  const int c_rts = header.col("RTS_stage");
  const int c_canon = header.col("all_canonical");
  const int c_percA = header.col("perc_A_downstream_TTS");
  const int c_indels_junc = header.col("n_indels_junc");
  const int c_cage = header.col("dist_to_CAGE_peak");
  const int c_polya = header.col("dist_to_polyA_site");
  const int c_nmd = header.col("predicted_NMD");
  const int c_polya_motif = header.col("polyA_motif_found");
  // SQANTI3 filter column has varied names across versions.
  int c_filter = header.col("filter_result");
  if (c_filter < 0) c_filter = header.col("filter");

  for (const char* name : {"chrom", "structural_category", "RTS_stage", "all_canonical",
                           "perc_A_downstream_TTS", "n_indels_junc"}) {
    if (!header.has(name)) table.missing_columns.emplace_back(name);
  }

  std::string line;
  std::size_t lineno = 1;  // header was line 1
  int warnings = 0;
  while (std::getline(in, line)) {
    ++lineno;
    chomp(line);
    if (line.empty()) continue;
    std::vector<std::string> f = split_tsv(line);

    SqantiRecord r;
    r.isoform = field_at(f, c_isoform);
    r.chrom = field_at(f, c_chrom);
    const std::string& strand_s = field_at(f, c_strand);
    r.strand = parse_strand(strand_s.empty() ? '.' : strand_s[0]);
    r.structural_category = field_at(f, c_cat);
    r.subcategory = field_at(f, c_subcat);
    r.associated_gene = field_at(f, c_gene);
    r.associated_transcript = field_at(f, c_tx);
    r.RTS_stage = field_at(f, c_rts);
    r.all_canonical = field_at(f, c_canon);
    r.perc_A_downstream_TTS = parse_double_or_na(field_at(f, c_percA), "perc_A_downstream_TTS", lineno, warnings);
    r.n_indels_junc = parse_int_or_na(field_at(f, c_indels_junc), "n_indels_junc", lineno, warnings);
    r.dist_to_CAGE_peak = parse_double_or_na(field_at(f, c_cage), "dist_to_CAGE_peak", lineno, warnings);
    r.dist_to_polyA_site = parse_double_or_na(field_at(f, c_polya), "dist_to_polyA_site", lineno, warnings);
    r.predicted_NMD = field_at(f, c_nmd);
    r.polyA_motif_found = field_at(f, c_polya_motif);
    r.filter_result = field_at(f, c_filter);

    ++table.category_counts[r.structural_category];
    table.records.push_back(std::move(r));
  }
  return table;
}

}  // namespace panisoguard
