cmake_minimum_required(VERSION 3.20)
project(panisoguard VERSION 0.0.7 LANGUAGES CXX C)

# --- Standards ---------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Options -----------------------------------------------------------------
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
# OBJECT, not STATIC: the core is internal (only the executable is installed), and an object
# library needs no archiver. conda-forge compilers >= 2.0 ship an unprefixed c++ but only a
# prefixed ar, so without system binutils CMake finds no `ar` and a STATIC build fails.
add_library(panisoguard_core OBJECT
  src/core/fingerprint.cpp
  src/core/consensus.cpp
  src/core/verdict.cpp
  src/core/rules.cpp
  src/core/adjudicator.cpp
  src/core/decoys.cpp
  src/evidence/bam_features.cpp
  src/evidence/variant_motif.cpp
  src/io/sj_tab_reader.cpp
  src/io/pangenome_reader.cpp
  src/io/gtf_reader.cpp
  src/io/bed12_reader.cpp
  src/io/sqanti_reader.cpp
  src/io/result_writer.cpp
)
target_include_directories(panisoguard_core PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(panisoguard_core SYSTEM PUBLIC "${CMAKE_SOURCE_DIR}/thirdparty")  # vendored: toml++
target_link_libraries(panisoguard_core PUBLIC HTSlib::hts)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(PANISOGUARD_WARNINGS -Wall -Wextra -Wshadow)
endif()
target_compile_options(panisoguard_core PRIVATE ${PANISOGUARD_WARNINGS})

# --- CLI ---------------------------------------------------------------------
add_executable(panisoguard
  src/cli/main.cpp
  src/cli/run_common.cpp
  src/cli/cmd_combine.cpp
  src/cli/cmd_adjudicate.cpp
  src/cli/cmd_example.cpp
  src/cli/cmd_decoys.cpp
)
target_link_libraries(panisoguard PRIVATE panisoguard_core)
target_compile_options(panisoguard PRIVATE ${PANISOGUARD_WARNINGS})
target_compile_definitions(panisoguard PRIVATE PANISOGUARD_VERSION="${PROJECT_VERSION}")
# `panisoguard example` looks for the bundled examples next to the installed binary
# (<prefix>/share/panisoguard) and falls back to the source checkout for a build-tree binary.
target_compile_definitions(panisoguard PRIVATE PANISOGUARD_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

install(TARGETS panisoguard RUNTIME DESTINATION bin)
# The bundled examples ship with the binary (so the bioconda package and the Docker image
# carry them); `panisoguard example <name>` copies one into a working folder.
install(DIRECTORY examples/tiny examples/multi_caller examples/worked_locus DESTINATION share/panisoguard/examples
        USE_SOURCE_PERMISSIONS PATTERN "output" EXCLUDE)

# --- Tests -------------------------------------------------------------------
if(PANISOGUARD_BUILD_TESTS)
  enable_testing()
  # M0 smoke test: the binary links and reports its version with exit code 0.
  add_test(NAME version_smoke COMMAND panisoguard --version)

  add_executable(panisoguard_tests
    tests/unit/test_main.cpp
    tests/unit/fingerprint_test.cpp
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
    tests/unit/variant_motif_test.cpp
    tests/unit/decoys_test.cpp
  )
  target_link_libraries(panisoguard_tests PRIVATE panisoguard_core)
  target_compile_options(panisoguard_tests PRIVATE ${PANISOGUARD_WARNINGS})
  target_include_directories(panisoguard_tests SYSTEM PRIVATE
    "${CMAKE_SOURCE_DIR}/thirdparty")  # vendored: catch2/catch.hpp
  target_compile_definitions(panisoguard_tests PRIVATE
    PANISOGUARD_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/tests/data/tiny")
  add_test(NAME unit COMMAND panisoguard_tests)
  # The README Quick start examples must keep matching their checked-in expected outputs.
  foreach(ex tiny multi_caller)
    add_test(NAME example_${ex}
      COMMAND bash "${CMAKE_SOURCE_DIR}/examples/${ex}/run.sh" "${CMAKE_BINARY_DIR}/example_${ex}_out")
    set_tests_properties(example_${ex} PROPERTIES ENVIRONMENT "PANISOGUARD=$<TARGET_FILE:panisoguard>")
  endforeach()
  # `panisoguard example` must find a bundled example, copy it, and the copy must still run.
  add_test(NAME example_copy
    COMMAND bash -c "rm -rf '${CMAKE_BINARY_DIR}/example_copy' && '$<TARGET_FILE:panisoguard>' example tiny '${CMAKE_BINARY_DIR}/example_copy' && bash '${CMAKE_BINARY_DIR}/example_copy/run.sh'")
  set_tests_properties(example_copy PROPERTIES ENVIRONMENT "PANISOGUARD=$<TARGET_FILE:panisoguard>")
  # --write-kept: the input GTF filtered to the kept verdicts (tiny: iso_known known,
  # iso_novel_supported confirmed, iso_novel_artifact artifact; the kept file must hold the
  # first two and not the third).
  add_test(NAME integration_write_kept
    COMMAND bash -c "mkdir -p '${CMAKE_BINARY_DIR}/write_kept' && cd '${CMAKE_SOURCE_DIR}/examples/tiny' && '$<TARGET_FILE:panisoguard>' adjudicate --classification data/classification.tsv --isoforms-gtf data/caller.gtf --ref-gtf data/reference.gtf --sj-tab data/short_read.SJ.tab --out-prefix '${CMAKE_BINARY_DIR}/write_kept/sample' --write-kept '${CMAKE_BINARY_DIR}/write_kept/kept.gtf' >/dev/null 2>&1 && grep -q 'transcript_id \"iso_novel_supported\"' '${CMAKE_BINARY_DIR}/write_kept/kept.gtf' && grep -q 'transcript_id \"iso_known\"' '${CMAKE_BINARY_DIR}/write_kept/kept.gtf' && ! grep -q 'iso_novel_artifact' '${CMAKE_BINARY_DIR}/write_kept/kept.gtf' && grep -q '^kept_verdicts' '${CMAKE_BINARY_DIR}/write_kept/sample.provenance.log'")

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
    # The shipped config/rules.default.toml must equal the built-in defaults (and the
    # consensus axis must fire) -- guards the silent-disable bug class.
    add_test(NAME integration_config_equivalence
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/benchmark/config_equivalence/ci_check.py"
              $<TARGET_FILE:panisoguard>)
    # The worked locus (examples/worked_locus): every verdict recomputed with the Python stdlib
    # and compared with the binary, and the 42-line rule replay against attribution.jsonl, so a
    # rule change cannot silently outdate the example or the replay.
    add_test(NAME integration_worked_locus
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/examples/worked_locus/check_worked_locus.py"
              $<TARGET_FILE:panisoguard>)
    add_test(NAME integration_rules_replay
      COMMAND ${Python3_EXECUTABLE}
              "${CMAKE_SOURCE_DIR}/examples/worked_locus/rules_replay.py"
              $<TARGET_FILE:panisoguard>)
  else()
    message(STATUS "PanIsoGuard: Python3 not found -- skipping truth-based integration tests")
  endif()
endif()
