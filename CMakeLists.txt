cmake_minimum_required(VERSION 3.20)
project(panisoguard VERSION 0.0.3 LANGUAGES CXX C)

# --- Standards ---------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Options -----------------------------------------------------------------
# In-process gbwtgraph/GBZ support is EXPERIMENTAL and OFF by default. The MVP
# (v1.0) consumes pangenome evidence only via subprocess/file tiers (GFA, rpvg
# text output), so a plain `cmake .. && make` needs nothing but htslib.
option(WITH_PANGENOME_LIB "Link in-process gbwtgraph/GBZ support (experimental)" OFF)
option(PANISOGUARD_BUILD_TESTS "Build the unit-test suite" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Prefer the active conda prefix when locating htslib.
if(DEFINED ENV{CONDA_PREFIX})
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

# --- htslib (required, never vendored, dynamically linked) -------------------
find_package(HTSlib 1.18 REQUIRED)

# Hard-fail early if the linked htslib predates the symbols PanIsoGuard relies on
# (hts_set_thread_pool was added in 1.18; sam_itr_regarray is used for batched
# multi-locus BAM probing).
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_INCLUDES "${HTSLIB_INCLUDE_DIRS}")
set(CMAKE_REQUIRED_LIBRARIES "${HTSLIB_LIBRARIES}")
check_cxx_source_compiles("
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
int main() {
  (void)&sam_itr_regarray;
  (void)&hts_set_thread_pool;
  return 0;
}
" HTSLIB_HAS_REQUIRED_SYMBOLS)
unset(CMAKE_REQUIRED_INCLUDES)
unset(CMAKE_REQUIRED_LIBRARIES)
if(NOT HTSLIB_HAS_REQUIRED_SYMBOLS)
  message(FATAL_ERROR
    "Linked htslib is missing sam_itr_regarray / hts_set_thread_pool (need >= 1.18). "
    "Install a recent release, e.g.  conda install -c bioconda 'htslib>=1.18'")
endif()
message(STATUS "PanIsoGuard: using htslib ${HTSLIB_VERSION} (${HTSLIB_LIBRARY})")

# --- Core library (M1) -------------------------------------------------------
add_library(panisoguard_core STATIC
  src/core/fingerprint.cpp
  src/core/interval_index.cpp
  src/core/consensus.cpp
  src/core/verdict.cpp
  src/core/rules.cpp
  src/core/adjudicator.cpp
  src/core/analysis.cpp
  src/evidence/bam_features.cpp
  src/evidence/variant_motif.cpp
  src/io/sj_tab_reader.cpp
  src/io/pangenome_reader.cpp
  src/io/gtf_reader.cpp
  src/io/bed12_reader.cpp
  src/io/sqanti_reader.cpp
  src/io/result_writer.cpp
)
target_include_directories(panisoguard_core PUBLIC
  "${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/thirdparty"   # vendored: cgranges/IITree.h
)
target_link_libraries(panisoguard_core PUBLIC HTSlib::hts)

# --- CLI ---------------------------------------------------------------------
add_executable(panisoguard
  src/cli/main.cpp
  src/cli/run_common.cpp
  src/cli/cmd_combine.cpp
  src/cli/cmd_adjudicate.cpp
  src/cli/cmd_benchmark.cpp
  src/cli/cmd_ablate.cpp
)
target_link_libraries(panisoguard PRIVATE panisoguard_core)
target_compile_definitions(panisoguard PRIVATE PANISOGUARD_VERSION="${PROJECT_VERSION}")
if(WITH_PANGENOME_LIB)
  target_compile_definitions(panisoguard PRIVATE HAVE_GBWTGRAPH=1)
endif()

install(TARGETS panisoguard RUNTIME DESTINATION bin)

# --- Tests -------------------------------------------------------------------
if(PANISOGUARD_BUILD_TESTS)
  enable_testing()
  # M0 smoke test: the binary links and reports its version with exit code 0.
  add_test(NAME version_smoke COMMAND panisoguard --version)

  add_executable(panisoguard_tests
    tests/unit/test_main.cpp
    tests/unit/fingerprint_test.cpp
    tests/unit/interval_index_test.cpp
    tests/unit/sj_tab_test.cpp
    tests/unit/pangenome_test.cpp
    tests/unit/gtf_test.cpp
    tests/unit/bed12_test.cpp
    tests/unit/io_robustness_test.cpp
    tests/unit/sqanti_test.cpp
    tests/unit/consensus_test.cpp
    tests/unit/rules_test.cpp
    tests/unit/adjudicator_test.cpp
    tests/unit/bam_features_test.cpp
    tests/unit/analysis_test.cpp
    tests/unit/variant_motif_test.cpp
  )
  target_link_libraries(panisoguard_tests PRIVATE panisoguard_core)
  target_include_directories(panisoguard_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/thirdparty")  # vendored: catch2/catch.hpp
  target_compile_definitions(panisoguard_tests PRIVATE
    PANISOGUARD_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data/tiny")
  add_test(NAME unit COMMAND panisoguard_tests)

  # --- Truth-based integration tests (drive the built binary) ----------------
  # Self-contained Python drivers: synthesize a labelled truth set, adjudicate it,
  # and assert the verdicts are correct. controlled_truth is pure text; the
  # synthetic_axes BAM sub-check self-skips when samtools is absent (CI has only
  # htslib), so both stay green in a minimal environment.
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME integration_controlled_truth
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/benchmark/controlled_truth/ci_check.py"
              $<TARGET_FILE:panisoguard>)
    add_test(NAME integration_synthetic_axes
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/benchmark/synthetic_axes/ci_check.py"
              $<TARGET_FILE:panisoguard>)
    # GATE-1 pangenome rescue on a committed sample of real HPRC v1.1 chr22 junctions
    # (rescue fires under population provenance; circularity firewall holds otherwise).
    add_test(NAME integration_pangenome_gate1
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/benchmark/pangenome/gate1_check.py"
              $<TARGET_FILE:panisoguard>
              "${CMAKE_SOURCE_DIR}/benchmark/pangenome/testdata/sample_junctions.tsv")
  else()
    message(STATUS "PanIsoGuard: Python3 not found -- skipping truth-based integration tests")
  endif()
endif()
