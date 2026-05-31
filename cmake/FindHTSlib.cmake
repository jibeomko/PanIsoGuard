# FindHTSlib.cmake -- locate htslib and expose the HTSlib::hts imported target.
#
# Search order honours CMAKE_PREFIX_PATH (which the top-level CMakeLists seeds
# with $CONDA_PREFIX) before falling back to system locations, so the conda
# htslib is preferred. htslib is dynamically linked and never vendored.
#
# Result variables:
#   HTSlib_FOUND            - TRUE if htslib was located
#   HTSLIB_INCLUDE_DIRS     - include directory containing htslib/sam.h
#   HTSLIB_LIBRARIES        - full path to libhts
#   HTSLIB_VERSION          - parsed from htslib/hts.h (e.g. 1.21)
#   HTSlib::hts             - imported target

find_path(HTSLIB_INCLUDE_DIR
  NAMES htslib/sam.h
  HINTS ${CMAKE_PREFIX_PATH} $ENV{CONDA_PREFIX} $ENV{HTSLIB_ROOT}
  PATH_SUFFIXES include
)

find_library(HTSLIB_LIBRARY
  NAMES hts
  HINTS ${CMAKE_PREFIX_PATH} $ENV{CONDA_PREFIX} $ENV{HTSLIB_ROOT}
  PATH_SUFFIXES lib lib64
)

set(HTSLIB_VERSION "")
if(HTSLIB_INCLUDE_DIR AND EXISTS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h")
  file(STRINGS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h" _pg_hts_ver
       REGEX "^#define[ \t]+HTS_VERSION[ \t]+[0-9]+")
  if(_pg_hts_ver MATCHES "HTS_VERSION[ \t]+([0-9]+)")
    set(_pg_ver_int "${CMAKE_MATCH_1}")
    math(EXPR _pg_major "${_pg_ver_int} / 100000")
    math(EXPR _pg_minor "(${_pg_ver_int} / 100) % 1000")
    set(HTSLIB_VERSION "${_pg_major}.${_pg_minor}")
  endif()
  unset(_pg_hts_ver)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HTSlib
  REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR
  VERSION_VAR HTSLIB_VERSION
)

if(HTSlib_FOUND)
  set(HTSLIB_LIBRARIES "${HTSLIB_LIBRARY}")
  set(HTSLIB_INCLUDE_DIRS "${HTSLIB_INCLUDE_DIR}")
  if(NOT TARGET HTSlib::hts)
    add_library(HTSlib::hts UNKNOWN IMPORTED)
    set_target_properties(HTSlib::hts PROPERTIES
      IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)
