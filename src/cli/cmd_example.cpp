// panisoguard example: copy a bundled example into a working folder.
//
// The examples ship with the binary (installed under <prefix>/share/panisoguard by
// `cmake --install`, so also in the bioconda package and the Docker image) and are
// located relative to the running executable; a build-tree binary falls back to the
// source checkout. They are copied rather than run in place because an install prefix
// is usually read-only and each example writes its output next to its run.sh.

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace fs = std::filesystem;

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

fs::path exe_dir() {
  std::error_code ec;
#ifdef __APPLE__
  char buf[4096];
  uint32_t n = sizeof buf;
  if (_NSGetExecutablePath(buf, &n) == 0) {
    fs::path p = fs::weakly_canonical(fs::path(buf), ec);
    if (!ec) return p.parent_path();
  }
#else
  fs::path p = fs::read_symlink("/proc/self/exe", ec);
  if (!ec) return p.parent_path();
#endif
  return {};
}

// <prefix>/bin/panisoguard -> <prefix>/share/panisoguard/<rel>; otherwise the source tree.
fs::path find_example(const Example& ex) {
  std::vector<fs::path> roots;
  const fs::path d = exe_dir();
  if (!d.empty()) roots.push_back(d / ".." / "share" / "panisoguard");
#ifdef PANISOGUARD_SOURCE_DIR
  roots.push_back(PANISOGUARD_SOURCE_DIR);
#endif
  for (const auto& r : roots) {
    std::error_code ec;
    const fs::path p = r / ex.rel;
    if (fs::is_directory(p, ec)) return fs::weakly_canonical(p, ec);
  }
  return {};
}

// Recursive copy that skips an example's own `output/` from an earlier local run.
bool copy_tree(const fs::path& src, const fs::path& dst, std::string& err) {
  std::error_code ec;
  fs::create_directories(dst, ec);
  if (ec) { err = ec.message(); return false; }
  for (const auto& e : fs::recursive_directory_iterator(src, ec)) {
    const fs::path rel = fs::relative(e.path(), src, ec);
    if (rel.begin() != rel.end() && *rel.begin() == "output") continue;
    const fs::path to = dst / rel;
    if (e.is_directory()) {
      fs::create_directories(to, ec);
    } else {
      fs::copy_file(e.path(), to, fs::copy_options::none, ec);   // keeps run.sh executable
    }
    if (ec) { err = to.string() + ": " + ec.message(); return false; }
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
  const fs::path src = find_example(*ex);
  if (src.empty()) {
    std::fprintf(stderr, "panisoguard example: '%s' is not installed next to this binary "
                         "(expected <prefix>/share/panisoguard/%s)\n", ex->name, ex->rel);
    return 1;
  }
  const fs::path dst = argc >= 4 ? fs::path(argv[3]) : fs::path(ex->name);
  std::error_code ec;
  if (fs::exists(dst, ec)) {
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
