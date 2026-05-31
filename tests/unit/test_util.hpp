#pragma once

#include <string>

#ifndef PANISOGUARD_TEST_DATA_DIR
#define PANISOGUARD_TEST_DATA_DIR "."
#endif

// Absolute path to a bundled tiny fixture.
inline std::string tiny(const char* filename) {
  return std::string(PANISOGUARD_TEST_DATA_DIR) + "/" + filename;
}
