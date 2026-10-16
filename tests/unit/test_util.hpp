#pragma once

#include <cstdlib>
#include <string>
#include <unistd.h>

#ifndef PANISOGUARD_TEST_DATA_DIR
#define PANISOGUARD_TEST_DATA_DIR "."
#endif

// Absolute path to a bundled tiny fixture.
inline std::string tiny(const char* filename) {
  return std::string(PANISOGUARD_TEST_DATA_DIR) + "/" + filename;
}

// A path under $TMPDIR (or /tmp), unique per process, for a file a test writes; nothing a
// test writes lands in the working directory. No <filesystem>: macOS < 10.15 lacks it.
inline std::string scratch(const std::string& name) {
  const char* t = std::getenv("TMPDIR");
  return std::string(t != nullptr && *t ? t : "/tmp") + "/pig_test_" + std::to_string(getpid()) + "_" + name;
}
