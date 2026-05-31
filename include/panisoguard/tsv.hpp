#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace panisoguard {

// Split an (already newline-stripped) line on tab. Bioinformatics TSVs are
// unquoted, so no quote handling is needed.
inline std::vector<std::string> split_tsv(const std::string& line) {
  std::vector<std::string> out;
  std::size_t start = 0;
  while (true) {
    const std::size_t tab = line.find('\t', start);
    if (tab == std::string::npos) {
      out.emplace_back(line.substr(start));
      break;
    }
    out.emplace_back(line.substr(start, tab - start));
    start = tab + 1;
  }
  return out;
}

// Strip a trailing CR and/or LF (CRLF-safe line handling).
inline void chomp(std::string& s) {
  while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
}

// Name-indexed header: tolerant column lookup that survives column reordering
// and appended columns (e.g. SQANTI3's appended filter_result columns).
class TsvHeader {
 public:
  explicit TsvHeader(const std::vector<std::string>& fields) {
    for (std::size_t i = 0; i < fields.size(); ++i) index_[fields[i]] = i;
    ncol_ = fields.size();
  }

  bool has(const std::string& name) const { return index_.count(name) != 0; }

  // Index of a column, or -1 if absent (caller decides not_evaluable handling).
  int col(const std::string& name) const {
    auto it = index_.find(name);
    return it == index_.end() ? -1 : static_cast<int>(it->second);
  }

  std::size_t ncol() const { return ncol_; }

 private:
  std::unordered_map<std::string, std::size_t> index_;
  std::size_t ncol_ = 0;
};

}  // namespace panisoguard
