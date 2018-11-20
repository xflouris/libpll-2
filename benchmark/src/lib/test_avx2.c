#include "common.h"
#include <x86intrin.h>

double test_avx2(const double *A, const double *B, double *C, size_t elems, size_t max_repeat) {
  __m256d result = _mm256_setzero_pd();

  for (size_t repeat = 0; repeat < max_repeat; repeat++) {
    for (size_t i = 0; i < elems; i+=ELEM_PER_AVX2_REGISTER) {
      __m256d a = _mm256_load_pd(A + i);
      __m256d b = _mm256_load_pd(B + i);
      result = _mm256_fmadd_pd(a, b, result);
    }
    _mm256_store_pd(C, result);
  }

  double reg[4] __attribute__( ( aligned ( 32 ) ) ) ;
  _mm256_store_pd(reg, result);
  return reg[0] + reg[1] + reg[2] + reg[3];
}