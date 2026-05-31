#pragma once

namespace panisoguard {

// Defined by CMake (-DPANISOGUARD_VERSION); falls back for non-CMake builds.
#ifndef PANISOGUARD_VERSION
#define PANISOGUARD_VERSION "0.0.0-unknown"
#endif

inline constexpr const char* kVersion = PANISOGUARD_VERSION;

}  // namespace panisoguard
