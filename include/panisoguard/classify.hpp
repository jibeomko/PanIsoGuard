#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "panisoguard/gtf.hpp"
#include "panisoguard/sqanti.hpp"
#include "panisoguard/variant_motif.hpp"

namespace panisoguard {

// A caller isoform as the classifier needs it: id, intron chain and genomic span.
struct CallerIsoform {
  std::string id;
  IntronChain chain;
  int64_t start = 0;   // 0-based half-open span
  int64_t end = 0;
};

// The classification rows PanIsoGuard reads from SQANTI3, computed here instead
// (`adjudicate --no-sqanti3`), with SQANTI3's definitions. Multi-exonic: full-splice_match
// (the chain is annotated), incomplete-splice_match (a contiguous part of an annotated chain
// with no retained intron), novel_in_catalog (every splice site annotated, the combination
// not), novel_not_in_catalog (a novel splice site), each in a locus annotated on its strand;
// antisense (only the other strand is annotated there); intergenic (nothing is). Mono-exonic:
// full-splice_match inside an exon of a mono-exonic transcript, incomplete-splice_match inside
// an exon of a multi-exonic one, genic overlapping an exon, genic_intron inside a transcript
// span only, else antisense or intergenic. With the reference genome,
// all_canonical reads every intron's splice motif and perc_A_downstream_TTS the A content of
// the 20 bases after the 3' end on the transcript strand, as SQANTI3 defines them; without it
// both are NA. RTS_stage and polyA_motif_found are NA: this mode does not compute them.
SqantiTable classify_isoforms(const std::vector<CallerIsoform>& isoforms, const Catalog& catalog,
                              const FastaFetcher* reference);

}  // namespace panisoguard
