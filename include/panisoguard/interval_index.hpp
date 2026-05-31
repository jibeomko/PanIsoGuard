#pragma once

#include <cstdint>
#include <vector>

#include "cgranges/IITree.h"

namespace panisoguard {

// RAII guard around cgranges' IITree.
//
// The bare IITree returns silently-empty overlaps if you query before calling
// index(); this wrapper enforces a build phase (add ...) followed by finalize()
// (which indexes exactly once) and makes overlap() throw if queried unindexed
// or add() throw if called after finalize(). Half-open intervals [start, end);
// payload is a uint32 (typically an index into a record vector).
class IntervalIndex {
 public:
  void add(int64_t start, int64_t end, uint32_t payload);
  void finalize();
  bool is_indexed() const { return indexed_; }
  std::size_t size() const { return tree_.size(); }

  // Payloads of all stored intervals overlapping [start, end). Throws if the
  // index has not been finalized.
  std::vector<uint32_t> overlap(int64_t start, int64_t end) const;

 private:
  IITree<int64_t, uint32_t> tree_;
  bool indexed_ = false;
};

}  // namespace panisoguard
