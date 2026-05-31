#include "panisoguard/interval_index.hpp"

#include <stdexcept>

namespace panisoguard {

void IntervalIndex::add(int64_t start, int64_t end, uint32_t payload) {
  if (indexed_) {
    throw std::logic_error("IntervalIndex::add called after finalize()");
  }
  tree_.add(start, end, payload);
}

void IntervalIndex::finalize() {
  tree_.index();
  indexed_ = true;
}

std::vector<uint32_t> IntervalIndex::overlap(int64_t start, int64_t end) const {
  if (!indexed_) {
    throw std::logic_error("IntervalIndex::overlap called before finalize()");
  }
  std::vector<std::size_t> hits;
  tree_.overlap(start, end, hits);
  std::vector<uint32_t> out;
  out.reserve(hits.size());
  for (std::size_t i : hits) {
    out.push_back(tree_.data(i));
  }
  return out;
}

}  // namespace panisoguard
