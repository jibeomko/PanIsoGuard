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
  }
  return "unknown";
}

const char* to_string(ConfidenceClass c) {
  switch (c) {
    case ConfidenceClass::kKnown:            return "KNOWN";
    case ConfidenceClass::kKnownFragment:    return "KNOWN_FRAGMENT";
    case ConfidenceClass::kNovelConfirmed:   return "NOVEL_CONFIRMED";
    case ConfidenceClass::kNovelProbable:    return "NOVEL_PROBABLE";
    case ConfidenceClass::kUnconfirmed:      return "UNCONFIRMED";
    case ConfidenceClass::kUnresolved:       return "UNRESOLVED";
    case ConfidenceClass::kReferenceBias:    return "REFERENCE_BIAS";
    case ConfidenceClass::kArtifact:         return "ARTIFACT";
  }
  return "UNRESOLVED";
}

}  // namespace panisoguard
