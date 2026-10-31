#pragma once

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <sys/stat.h>

namespace panisoguard {

// Parse a base-10 integer from a TSV field, throwing a context-rich error
// (file kind + line number + offending field) rather than a bare std::stoll/stoi
// "stoll" message, so a malformed row is diagnosable.
inline long long parse_int_field(const std::string& s, const char* what,
                                 const char* file_kind, std::size_t lineno) {
  try {
    std::size_t pos = 0;
    const long long v = std::stoll(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("trailing characters");
    return v;
  } catch (const std::exception&) {
    throw std::runtime_error(std::string(file_kind) + " line " + std::to_string(lineno) +
                             ": " + what + " is not an integer: \"" + s + "\"");
  }
}

// Split on one separator character; empty fields are kept ("a,,b" -> three fields, "" ->
// one empty field).
inline std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::size_t start = 0;
  while (true) {
    const std::size_t pos = s.find(sep, start);
    if (pos == std::string::npos) {
      out.emplace_back(s.substr(start));
      return out;
    }
    out.emplace_back(s.substr(start, pos - start));
    start = pos + 1;
  }
}

// Split an (already newline-stripped) line on tab. Bioinformatics TSVs are
// unquoted, so no quote handling is needed.
inline std::vector<std::string> split_tsv(const std::string& line) { return split(line, '\t'); }

// Strip a trailing CR and/or LF (CRLF-safe line handling).
inline void chomp(std::string& s) {
  while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
}

// Finish writing a file and fail if any write did not reach it: a full disk or a closed
// pipe otherwise leaves a truncated file behind a successful exit code.
inline void close_checked(std::ofstream& out, const std::string& path) {
  out.flush();
  if (!out) throw std::runtime_error("write failed: " + path);
  out.close();
  if (out.fail()) throw std::runtime_error("write failed: " + path);
}

// True when the two paths name the same file (through symlinks and hard links).
inline bool same_file(const std::string& a, const std::string& b) {
  struct stat sa, sb;
  if (stat(a.c_str(), &sa) != 0 || stat(b.c_str(), &sb) != 0) return false;
  return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
}

// JSON string escaping: quotes, backslashes and control characters.
inline std::string json_escape(const std::string& s) {
  std::string o;
  o.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '"':  o += "\\\""; break;
      case '\\': o += "\\\\"; break;
      case '\n': o += "\\n";  break;
      case '\r': o += "\\r";  break;
      case '\t': o += "\\t";  break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", c);
          o += buf;
        } else {
          o += c;
        }
    }
  }
  return o;
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
