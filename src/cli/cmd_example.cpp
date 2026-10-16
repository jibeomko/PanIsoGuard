// panisoguard example: copy a bundled example into a working folder.
//
// The examples ship with the binary (installed under <prefix>/share/panisoguard by
// `cmake --install`, so also in the bioconda package and the Docker image) and are
// located relative to the running executable; a build-tree binary falls back to the
// source checkout. They are copied rather than run in place because an install prefix
// is usually read-only and each example writes its output next to its run.sh.
//
// POSIX calls rather than <filesystem>: bioconda's osx-64 builds target macOS 10.9, where
// libc++'s std::filesystem is unavailable (it needs 10.15).

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace {

struct Example {
  const char* name;
  const char* rel;    // path under the share dir / source checkout
  const char* what;
  const char* next;   // what to run after copying (printed with the destination)
};

const Example kExamples[] = {
    {"tiny", "examples/tiny",
     "one caller: a known, a short-read-supported novel, and an artifact isoform", "run.sh"},
    {"multi_caller", "examples/multi_caller",
     "three callers: combine -> consensus verdict", "run.sh"},
    {"notes-data", "notes/data",
     "the study notes' toy gene (notes/ on GitHub walks through it)", nullptr},
};

bool is_dir(const std::string& p) {
  struct stat st;
  return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool exists(const std::string& p) {
  struct stat st;
  return lstat(p.c_str(), &st) == 0;
}

std::string canonical(const std::string& p) {
  char buf[PATH_MAX];
  return realpath(p.c_str(), buf) != nullptr ? std::string(buf) : std::string();
}

std::string parent(const std::string& p) {
  const std::size_t i = p.rfind('/');
  return i == std::string::npos ? std::string(".") : (i == 0 ? std::string("/") : p.substr(0, i));
}

// Directory of the running executable ("" if it cannot be found).
std::string exe_dir() {
#ifdef __APPLE__
  char buf[4096];
  uint32_t n = sizeof buf;
  if (_NSGetExecutablePath(buf, &n) == 0) {
    const std::string p = canonical(buf);
    if (!p.empty()) return parent(p);
  }
#else
  char buf[PATH_MAX];
  const ssize_t n = readlink("/proc/self/exe", buf, sizeof buf - 1);
  if (n > 0) { buf[n] = '\0'; return parent(buf); }
#endif
  return {};
}

// <prefix>/bin/panisoguard -> <prefix>/share/panisoguard/<rel>; otherwise the source tree.
std::string find_example(const Example& ex) {
  std::vector<std::string> roots;
  const std::string d = exe_dir();
  if (!d.empty()) roots.push_back(d + "/../share/panisoguard");
#ifdef PANISOGUARD_SOURCE_DIR
  roots.push_back(PANISOGUARD_SOURCE_DIR);
#endif
  for (const auto& r : roots) {
    const std::string p = r + "/" + ex.rel;
    if (is_dir(p)) return canonical(p);
  }
  return {};
}

bool copy_file(const std::string& from, const std::string& to, std::string& err) {
  std::ifstream in(from, std::ios::binary);
  std::ofstream out(to, std::ios::binary | std::ios::trunc);
  if (!in || !out) { err = to + ": cannot copy"; return false; }
  out << in.rdbuf();
  if (!out) { err = to + ": write failed"; return false; }
  struct stat st;
  if (stat(from.c_str(), &st) == 0) chmod(to.c_str(), st.st_mode & 0777);  // keeps run.sh executable
  return true;
}

// Recursive copy that skips an example's own `output/` from an earlier local run.
bool copy_tree(const std::string& src, const std::string& dst, std::string& err, bool top = true) {
  if (mkdir(dst.c_str(), 0755) != 0 && !is_dir(dst)) { err = dst + ": " + std::strerror(errno); return false; }
  DIR* d = opendir(src.c_str());
  if (d == nullptr) { err = src + ": " + std::strerror(errno); return false; }
  std::vector<std::string> names;
  while (dirent* e = readdir(d)) {
    const std::string n = e->d_name;
    if (n == "." || n == ".." || (top && n == "output")) continue;
    names.push_back(n);
  }
  closedir(d);
  for (const auto& n : names) {
    const std::string from = src + "/" + n, to = dst + "/" + n;
    if (is_dir(from)) { if (!copy_tree(from, to, err, false)) return false; }
    else if (!copy_file(from, to, err)) return false;
  }
  return true;
}

void usage() {
  std::fprintf(stderr,
      "Usage: panisoguard example                 list the bundled examples\n"
      "       panisoguard example <name> [DIR]    copy one into DIR (default ./<name>)\n\n"
      "Examples:\n");
  for (const auto& ex : kExamples) std::fprintf(stderr, "  %-14s %s\n", ex.name, ex.what);
}

}  // namespace

int cmd_example(int argc, char** argv) {
  if (argc >= 3 && (std::strcmp(argv[2], "-h") == 0 || std::strcmp(argv[2], "--help") == 0)) {
    usage();
    return 0;
  }
  if (argc < 3) {
    usage();
    return 0;
  }
  const Example* ex = nullptr;
  for (const auto& e : kExamples) if (std::strcmp(argv[2], e.name) == 0) ex = &e;
  if (ex == nullptr) {
    std::fprintf(stderr, "panisoguard example: unknown example '%s'\n\n", argv[2]);
    usage();
    return 1;
  }
  const std::string src = find_example(*ex);
  if (src.empty()) {
    std::fprintf(stderr, "panisoguard example: '%s' is not installed next to this binary "
                         "(expected <prefix>/share/panisoguard/%s)\n", ex->name, ex->rel);
    return 1;
  }
  const std::string dst = argc >= 4 ? argv[3] : ex->name;
  if (exists(dst)) {
    std::fprintf(stderr, "panisoguard example: '%s' already exists; choose another DIR\n", dst.c_str());
    return 1;
  }
  std::string err;
  if (!copy_tree(src, dst, err)) {
    std::fprintf(stderr, "panisoguard example: copy failed: %s\n", err.c_str());
    return 1;
  }
  std::printf("copied %s to %s\n", ex->name, dst.c_str());
  if (ex->next != nullptr) std::printf("next: %s/%s\n", dst.c_str(), ex->next);
  return 0;
}
