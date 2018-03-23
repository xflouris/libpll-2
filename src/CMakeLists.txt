cmake_minimum_required(VERSION 2.8)

include(CheckIncludeFile) 

set (LIBPLL_BASE_FLAGS "-Wall -Wsign-compare -D_GNU_SOURCE -O3 -fPIC ")

set (SSE_FLAGS "-msse3")
set (AVX_FLAGS "-mavx")
set (AVX2_FLAGS "-mfma -mavx2")

find_package(BISON)
find_package(FLEX)
set(LIBPLL_BISON_FLAGS "-y -d -p pll_utree_")
set(LIBPLL_FLEX_FLAGS "-P pll_utree_")
BISON_TARGET(parse_utree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/parse_utree.y ${CMAKE_CURRENT_BINARY_DIR}/parse_utree.c
  COMPILE_FLAGS ${LIBPLL_BISON_FLAGS})
FLEX_TARGET(lex_utree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/lex_utree.l ${CMAKE_CURRENT_BINARY_DIR}/lex_utree.c
  COMPILE_FLAGS ${LIBPLL_FLEX_FLAGS})
ADD_FLEX_BISON_DEPENDENCY(lex_utree_t parse_utree_t)
set(LIBPLL_BISON_FLAGS "-y -d -p pll_rtree_")
set(LIBPLL_FLEX_FLAGS "-P pll_rtree_")
BISON_TARGET(parse_rtree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/parse_rtree.y ${CMAKE_CURRENT_BINARY_DIR}/parse_rtree.c
  COMPILE_FLAGS ${LIBPLL_BISON_FLAGS})
FLEX_TARGET(lex_rtree_t 
  ${CMAKE_CURRENT_SOURCE_DIR}/lex_rtree.l ${CMAKE_CURRENT_BINARY_DIR}/lex_rtree.c
  COMPILE_FLAGS ${LIBPLL_FLEX_FLAGS})
ADD_FLEX_BISON_DEPENDENCY(lex_rtree_t parse_rtree_t)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#remove annoyting bison and flex warnings
SET_SOURCE_FILES_PROPERTIES( ${FLEX_lex_utree_t_OUTPUTS} PROPERTIES COMPILE_FLAGS -Wno-sign-compare )
SET_SOURCE_FILES_PROPERTIES( ${FLEX_lex_rtree_t_OUTPUTS} PROPERTIES COMPILE_FLAGS -Wno-sign-compare )

set(LIBPLL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/compress.c
  ${BISON_parse_utree_t_OUTPUTS}
  ${FLEX_lex_utree_t_OUTPUTS}
  ${BISON_parse_rtree_t_OUTPUTS}
  ${FLEX_lex_rtree_t_OUTPUTS}
  ${CMAKE_CURRENT_SOURCE_DIR}/core_derivatives.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_likelihood.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_partials.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_pmatrix.c
  ${CMAKE_CURRENT_SOURCE_DIR}/derivatives.c
  ${CMAKE_CURRENT_SOURCE_DIR}/fasta.c
  ${CMAKE_CURRENT_SOURCE_DIR}/fast_parsimony.c
  ${CMAKE_CURRENT_SOURCE_DIR}/gamma.c
  ${CMAKE_CURRENT_SOURCE_DIR}/hardware.c
  ${CMAKE_CURRENT_SOURCE_DIR}/likelihood.c
  ${CMAKE_CURRENT_SOURCE_DIR}/list.c
  ${CMAKE_CURRENT_SOURCE_DIR}/maps.c
  ${CMAKE_CURRENT_SOURCE_DIR}/models.c
  ${CMAKE_CURRENT_SOURCE_DIR}/output.c
  ${CMAKE_CURRENT_SOURCE_DIR}/parsimony.c
  ${CMAKE_CURRENT_SOURCE_DIR}/partials.c
  ${CMAKE_CURRENT_SOURCE_DIR}/phylip.c
  ${CMAKE_CURRENT_SOURCE_DIR}/pll.c
  ${CMAKE_CURRENT_SOURCE_DIR}/random.c
  ${CMAKE_CURRENT_SOURCE_DIR}/repeats.c
  ${CMAKE_CURRENT_SOURCE_DIR}/rtree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/stepwise.c
  ${CMAKE_CURRENT_SOURCE_DIR}/utree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/utree_moves.c
  ${CMAKE_CURRENT_SOURCE_DIR}/utree_svg.c
  )

file(GLOB LIBPLL_SSE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core_derivatives_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_likelihood_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_partials_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_pmatrix_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/fast_parsimony_sse.c
  )

file(GLOB LIBPLL_AVX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core_derivatives_avx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_likelihood_avx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_partials_avx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_pmatrix_avx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/fast_parsimony_avx.c
  )

file(GLOB LIBPLL_AVX2_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core_derivatives_avx2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_likelihood_avx2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_partials_avx2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/core_pmatrix_avx2.c
  ${CMAKE_CURRENT_SOURCE_DIR}/fast_parsimony_avx2.c
  )

# check that user did not disable simd
if (NOT DEFINED ENABLE_SSE)
  SET(ENABLE_SSE "True")
endif ()
if (NOT DEFINED ENABLE_AVX)
  SET(ENABLE_AVX "True")
endif ()
if (NOT DEFINED ENABLE_AVX2)
  SET(ENABLE_AVX2 "True")
endif ()

# check simd installed 
if (ENABLE_SSE)
  SET(_code " #include <immintrin.h>
  int main() {__m128d a = _mm_setzero_pd();  return 1;}")
  SET(_file ${CMAKE_CURRENT_BINARY_DIR}/testsse.c)
  FILE(WRITE "${_file}" "${_code}")
  TRY_COMPILE(SSE_COMPILED ${CMAKE_CURRENT_BINARY_DIR} ${_file}
    COMPILE_DEFINITIONS ${SSE_FLAGS})
  if (NOT SSE_COMPILED)
    message(STATUS "Disable sse simd, because not supported") 
    set(ENABLE_SSE "False")
  endif()
endif()
if (ENABLE_AVX)
  SET(_code " #include <immintrin.h>
  int main() {__m256d a = _mm256_setzero_pd ();  return 1;}")
  SET(_file ${CMAKE_CURRENT_BINARY_DIR}/testavx.c)
  FILE(WRITE "${_file}" "${_code}")
  TRY_COMPILE(AVX_COMPILED ${CMAKE_CURRENT_BINARY_DIR} ${_file}
    COMPILE_DEFINITIONS ${AVX_FLAGS})
  if (NOT AVX_COMPILED)
    message(STATUS "Disable avx simd, because not supported") 
    set(ENABLE_AVX "False")
  endif()
endif()
if (ENABLE_AVX2)
  SET(_code " #include <immintrin.h>
  int main() {__m256i a, b; b =  _mm256_abs_epi16(a); return 1;}")
  SET(_file ${CMAKE_CURRENT_BINARY_DIR}/testavx2.c)
  FILE(WRITE "${_file}" "${_code}")
  TRY_COMPILE(AVX2_COMPILED ${CMAKE_CURRENT_BINARY_DIR} ${_file}
    COMPILE_DEFINITIONS ${AVX2_FLAGS})
  if (NOT AVX2_COMPILED)
    message(STATUS "Disable avx2 simd, because not supported") 
    set(ENABLE_AVX2 "False")
  endif()
endif()


# set simd flags
if (ENABLE_SSE)
  add_definitions(-DHAVE_SSE)
  set(SIMD_FLAGS "${SIMD_FLAGS} ${SSE_FLAGS}")
  message(STATUS "SSE enabled. To disable it, run cmake with -DENABLE_SSE=false")
  set(LIBPLL_SOURCES ${LIBPLL_SOURCES} ${LIBPLL_SSE_SOURCES})
  SET_SOURCE_FILES_PROPERTIES( ${LIBPLL_SSE_SOURCES} PROPERTIES COMPILE_FLAGS ${SSE_FLAGS} )
endif ()
if (ENABLE_AVX)
  add_definitions(-DHAVE_AVX)
  set(SIMD_FLAGS "${SIMD_FLAGS} ${AVX_FLAGS}")
  message(STATUS "AVX enabled. To disable it, run cmake with -DENABLE_AVX=false")
  set(LIBPLL_SOURCES ${LIBPLL_SOURCES} ${LIBPLL_AVX_SOURCES})
  SET_SOURCE_FILES_PROPERTIES( ${LIBPLL_AVX_SOURCES} PROPERTIES COMPILE_FLAGS ${AVX_FLAGS} )
endif ()
if (ENABLE_AVX2)
  add_definitions(-DHAVE_AVX2)
  set(SIMD_FLAGS "${SIMD_FLAGS} ${AVX2_FLAGS}")
  message(STATUS "AVX2 enabled. To disable it, run cmake with -DENABLE_AVX2=false")
  set(LIBPLL_SOURCES ${LIBPLL_SOURCES} ${LIBPLL_AVX2_SOURCES})
  SET_SOURCE_FILES_PROPERTIES( ${LIBPLL_AVX2_SOURCES} PROPERTIES COMPILE_FLAGS ${AVX2_FLAGS} )
endif ()

add_definitions(-DHAVE_X86INTRIN_H)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LIBPLL_BASE_FLAGS}")

add_library(pll OBJECT ${LIBPLL_SOURCES})
set_property(TARGET pll PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(pll_shared SHARED $<TARGET_OBJECTS:pll>)
add_library(pll_static STATIC $<TARGET_OBJECTS:pll>)


