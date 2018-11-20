#include "pll.h"

#ifdef SIMD_COUNTING_MODE

size_t _mm256_setzero_pd_counter = 0;
size_t _mm256_store_pd_counter = 0;
size_t _mm256_fmadd_pd_counter = 0;
size_t _mm256_set1_pd_counter = 0;
size_t _mm256_load_pd_counter = 0;
size_t _mm256_unpackhi_pd_counter = 0;
size_t _mm256_unpacklo_pd_counter = 0;
size_t _mm256_add_pd_counter = 0;
size_t _mm256_mul_pd_counter = 0;
size_t _mm256_cmp_pd_counter = 0;
size_t _mm256_permute2f128_pd_counter = 0;
size_t _mm256_movemask_pd_counter = 0;
size_t _mm256_blend_pd_counter = 0;

#endif /* SIMD_COUNTING_MODE */