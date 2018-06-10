/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/* Include this file to enable SIMD instruction counting, Note that it will count every occurrence of AVX2 instructions */
#ifndef PLL_SIMD_COUNT_H
#define PLL_SIMD_COUNT_H

#ifdef SIMD_COUNTING_MODE
#include <stdlib.h>

/* Available counters counters */

extern size_t _mm256_setzero_pd_counter;
extern size_t _mm256_store_pd_counter;
extern size_t _mm256_fmadd_pd_counter;
extern size_t _mm256_set1_pd_counter;
extern size_t _mm256_load_pd_counter;
extern size_t _mm256_unpackhi_pd_counter;
extern size_t _mm256_unpacklo_pd_counter;
extern size_t _mm256_add_pd_counter;
extern size_t _mm256_mul_pd_counter;
extern size_t _mm256_cmp_pd_counter;
extern size_t _mm256_permute2f128_pd_counter;
extern size_t _mm256_movemask_pd_counter;
extern size_t _mm256_blend_pd_counter;

inline void pll_reset_simd_counter() {
  _mm256_setzero_pd_counter = 0;
  _mm256_store_pd_counter = 0;
  _mm256_fmadd_pd_counter = 0;
  _mm256_set1_pd_counter = 0;
  _mm256_load_pd_counter = 0;
  _mm256_unpackhi_pd_counter = 0;
  _mm256_unpacklo_pd_counter = 0;
  _mm256_add_pd_counter = 0;
  _mm256_mul_pd_counter = 0;
  _mm256_cmp_pd_counter = 0;
  _mm256_permute2f128_pd_counter = 0;
  _mm256_movemask_pd_counter = 0;
  _mm256_blend_pd_counter = 0;
}

/* check if appropriate compiler flags are present */
#if defined(__AVX__) || defined(__AVX2__)

inline __m256d _mm256_setzero_pd__inst__() {
  _mm256_setzero_pd_counter++;
  return _mm256_setzero_pd();
}
#undef _mm256_setzero_pd
#define _mm256_setzero_pd _mm256_setzero_pd__inst__

inline  __m256d _mm256_fmadd_pd__inst__(__m256d a, __m256d b, __m256d c) {
  _mm256_fmadd_pd_counter++;
  return _mm256_fmadd_pd(a, b, c);
}
#undef _mm256_fmadd_pd
#define _mm256_fmadd_pd _mm256_fmadd_pd__inst__

inline void _mm256_store_pd__inst__(double *a, __m256d b) {
  _mm256_store_pd_counter++;
  _mm256_store_pd(a, b);
}
#undef _mm256_store_pd
#define _mm256_store_pd _mm256_store_pd__inst__

inline __m256d _mm256_set1_pd__inst__(double v) {
  _mm256_set1_pd_counter++;
  return _mm256_set1_pd(v);
}
#undef _mm256_set1_pd
#define _mm256_set1_pd _mm256_set1_pd__inst__

inline __m256d _mm256_load_pd__inst__(double const *a) {
  _mm256_load_pd_counter++;
  return _mm256_load_pd(a);
}
#undef _mm256_load_pd
#define _mm256_load_pd _mm256_load_pd__inst__

inline __m256d _mm256_unpackhi_pd__inst__(__m256d a, __m256d b) {
  _mm256_unpackhi_pd_counter++;
  return _mm256_unpackhi_pd(a, b);
}
#undef _mm256_unpackhi_pd
#define _mm256_unpackhi_pd _mm256_unpackhi_pd__inst__

inline __m256d _mm256_unpacklo_pd__inst__(__m256d a, __m256d b) {
  _mm256_unpacklo_pd_counter++;
  return _mm256_unpacklo_pd(a, b);
}
#undef _mm256_unpacklo_pd
#define _mm256_unpacklo_pd _mm256_unpacklo_pd__inst__

inline __m256d _mm256_add_pd__inst__(__m256d a, __m256d b) {
  _mm256_add_pd_counter++;
  return _mm256_add_pd(a, b);
}
#undef _mm256_add_pd
#define _mm256_add_pd _mm256_add_pd__inst__

inline __m256d _mm256_mul_pd__inst__(__m256d a, __m256d b) {
  _mm256_mul_pd_counter++;
  return _mm256_mul_pd(a, b);
}
#undef _mm256_mul_pd
#define _mm256_mul_pd _mm256_mul_pd__inst__

inline __m256d _mm256_cmp_pd__inst__(__m256d a, __m256d b, const int predicate) {
  _mm256_cmp_pd_counter++;
  return _mm256_cmp_pd(a, b, predicate);
}
#undef _mm256_cmp_pd
#define _mm256_cmp_pd _mm256_cmp_pd__inst__

inline __m256d _mm256_permute2f128_pd__inst__(__m256d a, __m256d b, int control) {
  _mm256_permute2f128_pd_counter++;
  return _mm256_permute2f128_pd(a, b, control);
}
#undef _mm256_permute2f128_pd
#define _mm256_permute2f128_pd _mm256_permute2f128_pd__inst__

inline int _mm256_movemask_pd__inst__(__m256d a) {
  _mm256_movemask_pd_counter++;
  return _mm256_movemask_pd(a);
}
#undef _mm256_movemask_pd
#define _mm256_movemask_pd _mm256_movemask_pd__inst__

inline __m256d _mm256_blend_pd__inst__(__m256d a, __m256d b, const int mask) {
  _mm256_blend_pd_counter++;
  return _mm256_blend_pd(a, b, mask);
}
#undef _mm256_blend_pd
#define _mm256_blend_pd _mm256_blend_pd__inst__

#endif /* defined(__AVX__) || defined(__AVX2__) */

#endif /* SIMD_COUNTING_MODE */

#endif /* PLL_SIMD_COUNT_H */