#pragma once

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <climits>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>

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

// The path a file will have once written: a final symlink followed to its target (which need
// not exist yet), then the directory resolved through symlinks (when it exists) plus the file
// name. Lets two paths be compared before either file exists.
inline std::string resolved_path(std::string p) {
  for (int hop = 0; hop < 40; ++hop) {   // bounded like the kernel's symlink-loop limit
    struct stat st;
    if (lstat(p.c_str(), &st) != 0 || !S_ISLNK(st.st_mode)) break;
    char target[PATH_MAX];
    const ssize_t n = readlink(p.c_str(), target, sizeof(target) - 1);
    if (n < 0) break;
    target[n] = '\0';
    const std::size_t at = p.rfind('/');
    p = target[0] == '/' || at == std::string::npos ? std::string(target) : p.substr(0, at + 1) + target;
  }
  const std::size_t slash = p.rfind('/');
  const std::string dir = slash == std::string::npos ? "." : slash == 0 ? "/" : p.substr(0, slash);
  const std::string name = slash == std::string::npos ? p : p.substr(slash + 1);
  char buf[PATH_MAX];
  const char* real = realpath(dir.c_str(), buf);
  return (real != nullptr ? std::string(real) : dir) + "/" + name;
}

// Refuse an output path that names an input file, or another output: writing it would erase
// the input, or an earlier output, behind a successful exit. Existing files are compared by
// inode (a symlink or hard link to an input is the input); files not written yet by their
// resolved path (a symlink to a file not written yet is that file).
inline void refuse_overwriting_inputs(const std::vector<std::string>& outputs, const std::vector<std::string>& inputs) {
  auto clash = [](const std::string& a, const std::string& b) { return same_file(a, b) || resolved_path(a) == resolved_path(b); };
  for (std::size_t i = 0; i < outputs.size(); ++i) {
    if (outputs[i].empty()) continue;
    for (const auto& in : inputs) {
      if (!in.empty() && clash(outputs[i], in))
        throw std::runtime_error("output " + outputs[i] + " is the input file " + in + "; choose another path");
    }
    for (std::size_t j = i + 1; j < outputs.size(); ++j) {
      if (!outputs[j].empty() && clash(outputs[i], outputs[j]))
        throw std::runtime_error("outputs " + outputs[i] + " and " + outputs[j] + " are the same file");
    }
  }
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
