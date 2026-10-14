#pragma once

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "panisoguard/fingerprint.hpp"
#include "panisoguard/gtf.hpp"  // Catalog
#include "panisoguard/types.hpp"

namespace panisoguard {

// One integrated isoform: a single intron-chain identity onto which the native
// running IDs of one or more callers have been mapped.
struct ConsensusIsoform {
  std::string pig_id;             // stable PanIsoGuard id, e.g. "PIG.000001"
  std::string chrom;
  Strand strand = Strand::kUnknown;
  IntronChain chain;              // representative chain (introns identical within a group)
  bool monoexonic = false;        // monoexonic chains are NOT chain-fingerprintable -> never merged
  bool catalog_checked = false;   // true iff a catalog was supplied and the chain is multi-exon
  bool known_in_catalog = false;  // valid only when catalog_checked

  // caller label -> native running IDs that resolved to this chain
  std::map<std::string, std::vector<std::string>> support;

  int n_callers() const { return static_cast<int>(support.size()); }
  int n_isoforms() const {
    int n = 0;
    for (const auto& kv : support) n += static_cast<int>(kv.second.size());
    return n;
  }
};

// Accumulates caller isoforms and integrates them by caller-agnostic intron-chain
// fingerprint. Native running IDs are preserved; differently-named isoforms with
// the same splice chain collapse onto one ConsensusIsoform.
class ConsensusBuilder {
 public:
  // Register one isoform from a caller (e.g. caller="flair", native_id="1-1_...").
  void add(const std::string& caller, const std::string& native_id, const IntronChain& chain);

  // Materialize the integrated set with deterministic genomic-order PIG ids.
  // If catalog != nullptr, multi-exon chains are marked known/novel against it.
  std::vector<ConsensusIsoform> build(const Catalog* catalog = nullptr) const;

 private:
  struct Group {
    Fingerprint fp = 0;
    bool monoexonic = false;
    IntronChain chain;
    std::map<std::string, std::vector<std::string>> support;
  };
  std::map<std::string, Group> groups_;  // group key -> group
  long mono_counter_ = 0;                // keeps monoexonic isoforms distinct
};

// Write the caller-support matrix TSV (one row per ConsensusIsoform).
void write_caller_support_matrix(const std::string& path,
                                 const std::vector<ConsensusIsoform>& isoforms);

// A `combine` caller-support matrix as `adjudicate --caller-support` consumes it.
// The join key is the intron-chain fingerprint (column `fingerprint`, written since
// 0.0.6): two callers may reuse one native id ("transcript_1") for different chains, so
// a native-id join could hand one chain's support to another. The native-id map is kept
// for matrices without the column and for monoexonic isoforms (no chain fingerprint).
struct CallerSupport {
  std::unordered_map<Fingerprint, int> by_chain;       // fingerprint -> n_callers
  std::unordered_map<std::string, int> by_native_id;   // native running id -> n_callers (legacy join)
  int n_native_id_collisions = 0;  // native ids that appear under more than one chain
  bool has_fingerprints() const { return !by_chain.empty(); }
};
CallerSupport read_caller_support(const std::string& path);

}  // namespace panisoguard
