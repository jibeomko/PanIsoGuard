#include "panisoguard/fingerprint.hpp"

#include <vector>

namespace panisoguard {
namespace {

constexpr uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr uint64_t kFnvPrime = 1099511628211ULL;

inline void fnv_bytes(uint64_t& h, const void* data, size_t n) {
  const auto* p = static_cast<const unsigned char*>(data);
  for (size_t i = 0; i < n; ++i) {
    h ^= p[i];
    h *= kFnvPrime;
  }
}

}  // namespace

Fingerprint fingerprint_intron_chain(const IntronChain& chain) {
  uint64_t h = kFnvOffset;
  fnv_bytes(h, chain.chrom.data(), chain.chrom.size());
  const char s = strand_char(chain.strand);
  fnv_bytes(h, &s, 1);

  // Hash introns in sorted order so input order does not affect the result.
  std::vector<Junction> introns = chain.introns;
  std::sort(introns.begin(), introns.end());
  for (const auto& j : introns) {
    fnv_bytes(h, &j.start, sizeof(j.start));
    fnv_bytes(h, &j.end, sizeof(j.end));
  }
  return h;
}

}  // namespace panisoguard
