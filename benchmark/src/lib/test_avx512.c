#include "common.h"
#include <x86intrin.h>

double test_avx512(const double *A, const double *B, double *C, size_t elems, size_t max_repeat) {
  __m512d result = _mm512_setzero_pd();

  for (size_t repeat = 0; repeat < max_repeat; repeat++) {
    for (size_t i = 0; i < elems; i+= ELEM_PER_AVX512_REGISTER) {
      __m512d a = _mm512_load_pd(A + i);
      __m512d b = _mm512_load_pd(B + i);
      result = _mm512_fmadd_pd(a, b, result);
    }
    _mm512_store_pd(C, result);
  }

  double reg[8] __attribute__( ( aligned ( 64 ) ) ) ;
  _mm512_store_pd(reg, result);
  return reg[0] + reg[1] + reg[2] + reg[3] + reg[4] + reg[5] + reg[6] + reg[7];
}