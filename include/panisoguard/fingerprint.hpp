#pragma once

#include "panisoguard/types.hpp"

namespace panisoguard {

// 64-bit FNV-1a fingerprint over (chrom, strand, sorted introns).
//
// Only INTRONS are hashed, so the 5' start of the first exon and the 3' end of
// the last exon are excluded by construction -- mirroring SQANTI's FSM/ISM
// end-ignoring and making two isoforms with the same splice chain but different
// TSS/TES collapse to the same fingerprint.
//
// A monoexonic chain (no introns) hashes to chrom+strand only, so it is not a usable
// identity: ConsensusBuilder keeps every monoexonic isoform in its own group and the
// caller-support matrix writes NA for its fingerprint.
Fingerprint fingerprint_intron_chain(const IntronChain& chain);

}  // namespace panisoguard
