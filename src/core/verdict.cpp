#include "panisoguard/verdict.hpp"

namespace panisoguard {

const char* to_string(NoveltySupport s) {
  switch (s) {
    case NoveltySupport::kUnknown:     return "UNKNOWN";
    case NoveltySupport::kUnsupported: return "UNSUPPORTED";
    case NoveltySupport::kPartial:     return "PARTIAL";
    case NoveltySupport::kSupported:   return "SUPPORTED";
  }
  return "UNKNOWN";
}

const char* to_string(Mechanism m) {
  switch (m) {
    case Mechanism::kNone:           return "none";
    case Mechanism::kNoncanonical:   return "noncanonical";
    case Mechanism::kRtSwitch:       return "rt_switch";
    case Mechanism::kDegradation:    return "degradation";
    case Mechanism::kMapping:        return "mapping_or_repeat";
    case Mechanism::kVariant:        return "variant_created";
    case Mechanism::kPopulationKnown:return "population_known";
    case Mechanism::kUnknown:        return "unknown";
  }
  return "unknown";
}

const char* to_string(ConfidenceClass c) {
  switch (c) {
    case ConfidenceClass::kHighConfKnown:             return "HIGH_CONF_KNOWN";
    case ConfidenceClass::kHighConfNovel:             return "HIGH_CONF_NOVEL";
    case ConfidenceClass::kMediumConfNovel:           return "MEDIUM_CONF_NOVEL";
    case ConfidenceClass::kLowConfPartial:            return "LOW_CONF_PARTIAL";
    case ConfidenceClass::kPanRefRescuedFalseNovel:   return "PAN_REF_RESCUED_FALSE_NOVEL";
    case ConfidenceClass::kHaplotypeRescuedKnownLike: return "HAPLOTYPE_RESCUED_KNOWN_LIKE";
    case ConfidenceClass::kAmbiguous:                 return "AMBIGUOUS";
    case ConfidenceClass::kArtifact:                  return "ARTIFACT";
  }
  return "AMBIGUOUS";
}

}  // namespace panisoguard
