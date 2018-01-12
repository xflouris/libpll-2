/*
    Copyright (C) 2016 Tomas Flouri, Diego Darriba, Alexey Kozlov

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
#include <limits.h>
#include "pll.h"

inline double reduce_add_pd(const __m512d zmm) {
  __m256d low = _mm512_castpd512_pd256(zmm);
  __m256d high = _mm512_extractf64x4_pd(zmm, 1);

  __m256d a = _mm256_add_pd(low, high);
  __m256d t1 = _mm256_hadd_pd(a, a);
  __m128d t2 = _mm256_extractf128_pd(t1, 1);
  __m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1), t2);
  return _mm_cvtsd_f64(t3);
}

#define COMPUTE_II_QCOL(q, offset) \
/* row 0 */ \
v_mat    = _mm512_load_pd(lm0 + (offset)); \
v_lterm0 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm0); \
v_mat    = _mm512_load_pd(rm0 + (offset)); \
v_rterm0 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm0); \
\
/* row 1 */ \
v_mat    = _mm512_load_pd(lm1 + (offset)); \
v_lterm1 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm1); \
v_mat    = _mm512_load_pd(rm1 + (offset)); \
v_rterm1 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm1); \
\
/* row 2 */ \
v_mat    = _mm512_load_pd(lm2 + (offset)); \
v_lterm2 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm2); \
v_mat    = _mm512_load_pd(rm2 + (offset)); \
v_rterm2 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm2); \
\
/* row 3 */ \
v_mat    = _mm512_load_pd(lm3 + (offset)); \
v_lterm3 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm3); \
v_mat    = _mm512_load_pd(rm3 + (offset)); \
v_rterm3 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm3); \
\
/* row 4 */ \
v_mat    = _mm512_load_pd(lm4 + (offset)); \
v_lterm4 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm4); \
v_mat    = _mm512_load_pd(rm4 + (offset)); \
v_rterm4 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm4); \
\
/* row 5 */ \
v_mat    = _mm512_load_pd(lm5 + (offset)); \
v_lterm5 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm5); \
v_mat    = _mm512_load_pd(rm5 + (offset)); \
v_rterm5 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm5); \
\
/* row 6 */ \
v_mat    = _mm512_load_pd(lm6 + (offset)); \
v_lterm6 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm6); \
v_mat    = _mm512_load_pd(rm6 + (offset)); \
v_rterm6 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm6); \
\
/* row 7 */ \
v_mat    = _mm512_load_pd(lm7 + (offset)); \
v_lterm7 = _mm512_fmadd_pd(v_mat, v_lclv[q], v_lterm7); \
v_mat    = _mm512_load_pd(rm7 + (offset)); \
v_rterm7 = _mm512_fmadd_pd(v_mat, v_rclv[q], v_rterm7);

#define COMPUTE_STATE_PART(i, j, k) \
v_inv_eigenvecs = _mm512_set1_pd(tt_inv_eigenvecs[(i) * states * states \
                                                          + (j) * states \
                                                          + (k)]); \
v_eigenvecs = _mm512_set1_pd(tt_eigenvecs[(i) * states * states \
                                                  + (j) * states \
                                                  + (k)]); \
v_lefterm = _mm512_fmadd_pd(v_clvp[k], v_inv_eigenvecs, v_lefterm); \
v_righterm = _mm512_fmadd_pd(v_clvc[k], v_eigenvecs, v_righterm); \

PLL_EXPORT int pll_core_update_sumtable_ii_20x20_avx512f(unsigned int sites,
                                                         unsigned int rate_cats,
                                                         const double *clvp,
                                                         const double *clvc,
                                                         const unsigned int *parent_scaler,
                                                         const unsigned int *child_scaler,
                                                         double *const *eigenvecs,
                                                         double *const *inv_eigenvecs,
                                                         double *const *freqs,
                                                         double *sumtable,
                                                         unsigned int attrib) {
  const double *t_clvp = clvp;
  const double *t_clvc = clvc;

  double *sum = sumtable;

  unsigned int states = 20;
  unsigned int states_padded = (states + 7) & (0xFFFFFFFF - 7);

/* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int *rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

/* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling) {
    rate_scalings = (unsigned int *) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (unsigned int i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i) {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

/* padded eigenvecs */
  double *tt_eigenvecs = (double *) pll_aligned_alloc(
          (states * states * rate_cats) * sizeof(double),
          PLL_ALIGNMENT_AVX512F);

  if (!tt_eigenvecs) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory for tt_eigenvecs");
    return PLL_FAILURE;
  }

/* transposed padded inv_eigenvecs */
  double *tt_inv_eigenvecs = (double *) pll_aligned_alloc(
          (states * states * rate_cats) * sizeof(double),
          PLL_ALIGNMENT_AVX512F);

  if (!tt_inv_eigenvecs) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

/* broadcast eigenvecs matrices and multiply with frequencies */
  for (unsigned int i = 0; i < rate_cats; ++i) {
    for (unsigned int j = 0; j < states; ++j)
      for (unsigned int k = 0; k < states; ++k) {
        tt_inv_eigenvecs[i * states * states +
                         j * states +
                         k] = inv_eigenvecs[i][k * states_padded + j] * freqs[i][k];
        tt_eigenvecs[i * states * states +
                     j * states +
                     k] = eigenvecs[i][j * states_padded + k];
      }
  }

  //TODO
  if (per_rate_scaling) {
    printf("Per rate scaling not supported in AVX512");
    exit(1);
  }

  __m512i v_index = _mm512_setr_epi64(0,
                                      1 * rate_cats * states_padded,
                                      2 * rate_cats * states_padded,
                                      3 * rate_cats * states_padded,
                                      4 * rate_cats * states_padded,
                                      5 * rate_cats * states_padded,
                                      6 * rate_cats * states_padded,
                                      7 * rate_cats * states_padded);
  __mmask8 gather_mask = 0x0F;

/* build sumtable */
  for (unsigned int n = 0; n < sites; n += ELEM_PER_AVX515_REGISTER) {
    for (unsigned int i = 0; i < rate_cats; ++i) {
      __m512d v_clvp[states];
      __m512d v_clvc[states];

      if (n < 16) {
        v_clvp[0] = _mm512_i64gather_pd(v_index, t_clvp, sizeof(double));
        v_clvc[0] = _mm512_i64gather_pd(v_index, t_clvc, sizeof(double));
        v_clvp[1] = _mm512_i64gather_pd(v_index, t_clvp + 1, sizeof(double));
        v_clvc[1] = _mm512_i64gather_pd(v_index, t_clvc + 1, sizeof(double));
        v_clvp[2] = _mm512_i64gather_pd(v_index, t_clvp + 2, sizeof(double));
        v_clvc[2] = _mm512_i64gather_pd(v_index, t_clvc + 2, sizeof(double));
        v_clvp[3] = _mm512_i64gather_pd(v_index, t_clvp + 3, sizeof(double));
        v_clvc[3] = _mm512_i64gather_pd(v_index, t_clvc + 3, sizeof(double));
        v_clvp[4] = _mm512_i64gather_pd(v_index, t_clvp + 4, sizeof(double));
        v_clvc[4] = _mm512_i64gather_pd(v_index, t_clvc + 4, sizeof(double));
        v_clvp[5] = _mm512_i64gather_pd(v_index, t_clvp + 5, sizeof(double));
        v_clvc[5] = _mm512_i64gather_pd(v_index, t_clvc + 5, sizeof(double));
        v_clvp[6] = _mm512_i64gather_pd(v_index, t_clvp + 6, sizeof(double));
        v_clvc[6] = _mm512_i64gather_pd(v_index, t_clvc + 6, sizeof(double));
        v_clvp[7] = _mm512_i64gather_pd(v_index, t_clvp + 7, sizeof(double));
        v_clvc[7] = _mm512_i64gather_pd(v_index, t_clvc + 7, sizeof(double));
        v_clvp[8] = _mm512_i64gather_pd(v_index, t_clvp + 8, sizeof(double));
        v_clvc[8] = _mm512_i64gather_pd(v_index, t_clvc + 8, sizeof(double));
        v_clvp[9] = _mm512_i64gather_pd(v_index, t_clvp + 9, sizeof(double));
        v_clvc[9] = _mm512_i64gather_pd(v_index, t_clvc + 9, sizeof(double));
        v_clvp[10] = _mm512_i64gather_pd(v_index, t_clvp + 10, sizeof(double));
        v_clvc[10] = _mm512_i64gather_pd(v_index, t_clvc + 10, sizeof(double));
        v_clvp[11] = _mm512_i64gather_pd(v_index, t_clvp + 11, sizeof(double));
        v_clvc[11] = _mm512_i64gather_pd(v_index, t_clvc + 11, sizeof(double));
        v_clvp[12] = _mm512_i64gather_pd(v_index, t_clvp + 12, sizeof(double));
        v_clvc[12] = _mm512_i64gather_pd(v_index, t_clvc + 12, sizeof(double));
        v_clvp[13] = _mm512_i64gather_pd(v_index, t_clvp + 13, sizeof(double));
        v_clvc[13] = _mm512_i64gather_pd(v_index, t_clvc + 13, sizeof(double));
        v_clvp[14] = _mm512_i64gather_pd(v_index, t_clvp + 14, sizeof(double));
        v_clvc[14] = _mm512_i64gather_pd(v_index, t_clvc + 14, sizeof(double));
        v_clvp[15] = _mm512_i64gather_pd(v_index, t_clvp + 15, sizeof(double));
        v_clvc[15] = _mm512_i64gather_pd(v_index, t_clvc + 15, sizeof(double));
        v_clvp[16] = _mm512_i64gather_pd(v_index, t_clvp + 16, sizeof(double));
        v_clvc[16] = _mm512_i64gather_pd(v_index, t_clvc + 16, sizeof(double));
        v_clvp[17] = _mm512_i64gather_pd(v_index, t_clvp + 17, sizeof(double));
        v_clvc[17] = _mm512_i64gather_pd(v_index, t_clvc + 17, sizeof(double));
        v_clvp[18] = _mm512_i64gather_pd(v_index, t_clvp + 18, sizeof(double));
        v_clvc[18] = _mm512_i64gather_pd(v_index, t_clvc + 18, sizeof(double));
        v_clvp[19] = _mm512_i64gather_pd(v_index, t_clvp + 19, sizeof(double));
        v_clvc[19] = _mm512_i64gather_pd(v_index, t_clvc + 19, sizeof(double));
      } else {
        v_clvp[0] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp, sizeof(double));
        v_clvc[0] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc, sizeof(double));
        v_clvp[1] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 1, sizeof(double));
        v_clvc[1] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 1, sizeof(double));
        v_clvp[2] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 2, sizeof(double));
        v_clvc[2] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 2, sizeof(double));
        v_clvp[3] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 3, sizeof(double));
        v_clvc[3] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 3, sizeof(double));
        v_clvp[4] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 4, sizeof(double));
        v_clvc[4] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 4, sizeof(double));
        v_clvp[5] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 5, sizeof(double));
        v_clvc[5] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 5, sizeof(double));
        v_clvp[6] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 6, sizeof(double));
        v_clvc[6] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 6, sizeof(double));
        v_clvp[7] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 7, sizeof(double));
        v_clvc[7] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 7, sizeof(double));
        v_clvp[8] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 8, sizeof(double));
        v_clvc[8] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 8, sizeof(double));
        v_clvp[9] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 9, sizeof(double));
        v_clvc[9] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 9, sizeof(double));
        v_clvp[10] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 10, sizeof(double));
        v_clvc[10] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 10, sizeof(double));
        v_clvp[11] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 11, sizeof(double));
        v_clvc[11] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 11, sizeof(double));
        v_clvp[12] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 12, sizeof(double));
        v_clvc[12] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 12, sizeof(double));
        v_clvp[13] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 13, sizeof(double));
        v_clvc[13] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 13, sizeof(double));
        v_clvp[14] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 14, sizeof(double));
        v_clvc[14] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 14, sizeof(double));
        v_clvp[15] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 15, sizeof(double));
        v_clvc[15] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 15, sizeof(double));
        v_clvp[16] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 16, sizeof(double));
        v_clvc[16] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 16, sizeof(double));
        v_clvc[17] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 17, sizeof(double));
        v_clvp[17] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 17, sizeof(double));
        v_clvp[18] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 18, sizeof(double));
        v_clvc[18] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 18, sizeof(double));
        v_clvp[19] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvp + 19, sizeof(double));
        v_clvc[19] = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), gather_mask, v_index, t_clvc + 19, sizeof(double));
      }

      for (unsigned int j = 0; j < states; ++j) {
        __m512d v_lefterm = _mm512_setzero_pd();
        __m512d v_righterm = _mm512_setzero_pd();

        __m512d v_inv_eigenvecs;
        __m512d v_eigenvecs;

        COMPUTE_STATE_PART(i, j, 0);
        COMPUTE_STATE_PART(i, j, 1);
        COMPUTE_STATE_PART(i, j, 2);
        COMPUTE_STATE_PART(i, j, 3);
        COMPUTE_STATE_PART(i, j, 4);
        COMPUTE_STATE_PART(i, j, 5);
        COMPUTE_STATE_PART(i, j, 6);
        COMPUTE_STATE_PART(i, j, 7);
        COMPUTE_STATE_PART(i, j, 8);
        COMPUTE_STATE_PART(i, j, 9);
        COMPUTE_STATE_PART(i, j, 10);
        COMPUTE_STATE_PART(i, j, 11);
        COMPUTE_STATE_PART(i, j, 12);
        COMPUTE_STATE_PART(i, j, 13);
        COMPUTE_STATE_PART(i, j, 14);
        COMPUTE_STATE_PART(i, j, 15);
        COMPUTE_STATE_PART(i, j, 16);
        COMPUTE_STATE_PART(i, j, 17);
        COMPUTE_STATE_PART(i, j, 18);
        COMPUTE_STATE_PART(i, j, 19);

        __m512d v_sum = _mm512_mul_pd(v_lefterm, v_righterm);

        _mm512_store_pd(sum, v_sum);
        sum += ELEM_PER_AVX515_REGISTER;
      }
      t_clvc += states_padded;
      t_clvp += states_padded;
    }
    //pointers already moved one site ahead, move another 7 sites forward,
    //so we start at the date of the 9th state
    t_clvc += rate_cats * states_padded * (ELEM_PER_AVX515_REGISTER - 1);
    t_clvp += rate_cats * states_padded * (ELEM_PER_AVX515_REGISTER - 1);
  }

  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii_avx512f(unsigned int states,
                                                   unsigned int sites,
                                                   unsigned int rate_cats,
                                                   const double *clvp,
                                                   const double *clvc,
                                                   const unsigned int *parent_scaler,
                                                   const unsigned int *child_scaler,
                                                   double *const *eigenvecs,
                                                   double *const *inv_eigenvecs,
                                                   double *const *freqs,
                                                   double *sumtable,
                                                   unsigned int attrib) {
  unsigned int i, j, k, n;

  /* build sumtable */
  double *sum = sumtable;

  const double *t_clvp = clvp;
  const double *t_clvc = clvc;
  double *t_freqs;

  /* dedicated functions for 4x4 and 20x20 matrices */
  if (states == 4) {
    /* call AVX variant */
    return pll_core_update_sumtable_ii_avx(states,
                                           sites,
                                           rate_cats,
                                           clvp,
                                           clvc,
                                           parent_scaler,
                                           child_scaler,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs,
                                           sumtable,
                                           attrib);
  } else if (states == 20) {
    return pll_core_update_sumtable_ii_20x20_avx512f(sites,
                                                     rate_cats,
                                                     clvp,
                                                     clvc,
                                                     parent_scaler,
                                                     child_scaler,
                                                     eigenvecs,
                                                     inv_eigenvecs,
                                                     freqs,
                                                     sumtable,
                                                     attrib);
  }

  unsigned int states_padded = (states + 7) & (0xFFFFFFFF - 7);

  __m512i permute_mask = _mm512_setr_epi64(0 | 4,
                                           0 | 5,
                                           0 | 6,
                                           0 | 7,
                                           8 | 0,
                                           8 | 1,
                                           8 | 2,
                                           8 | 3);
  __m512i permute_mask_final_stage = _mm512_setr_epi64(0 | 2,
                                                       0 | 3,
                                                       8 | 0,
                                                       8 | 1,
                                                       0 | 6,
                                                       0 | 7,
                                                       8 | 4,
                                                       8 | 5);

  /* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int *rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m512d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling) {
    rate_scalings = (unsigned int *) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings) {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i) {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm512_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  double *tt_eigenvecs = (double *) pll_aligned_alloc(
          (states_padded * states_padded * rate_cats) * sizeof(double),
          PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory for tt_eigenvecs");
    return PLL_FAILURE;
  }

  /* transposed padded inv_eigenvecs */
  double *tt_inv_eigenvecs = (double *) pll_aligned_alloc(
          (states_padded * states_padded * rate_cats) * sizeof(double),
          PLL_ALIGNMENT_AVX);

  if (!tt_inv_eigenvecs) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i) {
    t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k) {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
                         + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
                     + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++) {

    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling) {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i) {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n * rate_cats + i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n * rate_cats + i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i) {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    const double *c_eigenvecs = tt_eigenvecs;
    const double *ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i) {
      for (j = 0; j < states_padded; j += ELEM_PER_AVX515_REGISTER) {
        /* point to the eight rows of the eigenvecs matrix */
        const double *em0 = c_eigenvecs;
        const double *em1 = em0 + states_padded;
        const double *em2 = em1 + states_padded;
        const double *em3 = em2 + states_padded;
        const double *em4 = em3 + states_padded;
        const double *em5 = em4 + states_padded;
        const double *em6 = em5 + states_padded;
        const double *em7 = em6 + states_padded;
        c_eigenvecs += ELEM_PER_AVX515_REGISTER * states_padded;

        /* point to the eight rows of the inv_eigenvecs matrix */
        const double *im0 = ct_inv_eigenvecs;
        const double *im1 = im0 + states_padded;
        const double *im2 = im1 + states_padded;
        const double *im3 = im2 + states_padded;
        const double *im4 = im3 + states_padded;
        const double *im5 = im4 + states_padded;
        const double *im6 = im5 + states_padded;
        const double *im7 = im6 + states_padded;
        ct_inv_eigenvecs += ELEM_PER_AVX515_REGISTER * states_padded;

        __m512d v_lefterm0 = _mm512_setzero_pd();
        __m512d v_righterm0 = _mm512_setzero_pd();
        __m512d v_lefterm1 = _mm512_setzero_pd();
        __m512d v_righterm1 = _mm512_setzero_pd();
        __m512d v_lefterm2 = _mm512_setzero_pd();
        __m512d v_righterm2 = _mm512_setzero_pd();
        __m512d v_lefterm3 = _mm512_setzero_pd();
        __m512d v_righterm3 = _mm512_setzero_pd();
        __m512d v_lefterm4 = _mm512_setzero_pd();
        __m512d v_righterm4 = _mm512_setzero_pd();
        __m512d v_lefterm5 = _mm512_setzero_pd();
        __m512d v_righterm5 = _mm512_setzero_pd();
        __m512d v_lefterm6 = _mm512_setzero_pd();
        __m512d v_righterm6 = _mm512_setzero_pd();
        __m512d v_lefterm7 = _mm512_setzero_pd();
        __m512d v_righterm7 = _mm512_setzero_pd();

        __m512d v_eigen;
        __m512d v_clvp;
        __m512d v_clvc;

        for (k = 0; k < states_padded; k += ELEM_PER_AVX515_REGISTER) {
          v_clvp = _mm512_load_pd(t_clvp + k);
          v_clvc = _mm512_load_pd(t_clvc + k);

          /* row 0 */
          v_eigen = _mm512_load_pd(im0 + k);
          v_lefterm0 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm0);

          v_eigen = _mm512_load_pd(em0 + k);
          v_righterm0 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm0);

          /* row 1 */
          v_eigen = _mm512_load_pd(im1 + k);
          v_lefterm1 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm1);

          v_eigen = _mm512_load_pd(em1 + k);
          v_righterm1 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm1);

          /* row 2 */
          v_eigen = _mm512_load_pd(im2 + k);
          v_lefterm2 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm2);

          v_eigen = _mm512_load_pd(em2 + k);
          v_righterm2 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm2);

          /* row 3 */
          v_eigen = _mm512_load_pd(im3 + k);
          v_lefterm3 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm3);

          v_eigen = _mm512_load_pd(em3 + k);
          v_righterm3 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm3);

          /* row 4 */
          v_eigen = _mm512_load_pd(im4 + k);
          v_lefterm4 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm4);

          v_eigen = _mm512_load_pd(em4 + k);
          v_righterm4 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm4);

          /* row 5 */
          v_eigen = _mm512_load_pd(im5 + k);
          v_lefterm5 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm5);

          v_eigen = _mm512_load_pd(em5 + k);
          v_righterm5 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm5);

          /* row 6 */
          v_eigen = _mm512_load_pd(im6 + k);
          v_lefterm6 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm6);

          v_eigen = _mm512_load_pd(em6 + k);
          v_righterm6 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm6);

          /* row 7 */
          v_eigen = _mm512_load_pd(im7 + k);
          v_lefterm7 = _mm512_fmadd_pd(v_eigen, v_clvp, v_lefterm7);

          v_eigen = _mm512_load_pd(em7 + k);
          v_righterm7 = _mm512_fmadd_pd(v_eigen, v_clvc, v_righterm7);
        }

        /* compute lefterm */
        __m512d xmm0 = _mm512_add_pd(_mm512_unpackhi_pd(v_lefterm0, v_lefterm1),
                                     _mm512_unpacklo_pd(v_lefterm0, v_lefterm1));
        __m512d xmm1 = _mm512_add_pd(_mm512_unpackhi_pd(v_lefterm2, v_lefterm3),
                                     _mm512_unpacklo_pd(v_lefterm2, v_lefterm3));
        __m512d xmm2 = _mm512_add_pd(_mm512_unpackhi_pd(v_lefterm4, v_lefterm5),
                                     _mm512_unpacklo_pd(v_lefterm4, v_lefterm5));
        __m512d xmm3 = _mm512_add_pd(_mm512_unpackhi_pd(v_lefterm6, v_lefterm7),
                                     _mm512_unpacklo_pd(v_lefterm6, v_lefterm7));

        __m512d ymm0 = _mm512_add_pd(_mm512_permutex2var_pd(xmm0, permute_mask, xmm2),
                                     _mm512_mask_blend_pd(0xF0, xmm0, xmm2));

        __m512d ymm1 = _mm512_add_pd(_mm512_permutex2var_pd(xmm1, permute_mask, xmm3),
                                     _mm512_mask_blend_pd(0xF0, xmm1, xmm3));

        __m512d v_lefterm_sum = _mm512_add_pd(_mm512_permutex2var_pd(ymm0,
                                                                     permute_mask_final_stage,
                                                                     ymm1),
                                              _mm512_mask_blend_pd(0xCC, ymm0, ymm1));

        /* compute righterm */
        xmm0 = _mm512_add_pd(_mm512_unpackhi_pd(v_righterm0, v_righterm1),
                             _mm512_unpacklo_pd(v_righterm0, v_righterm1));
        xmm1 = _mm512_add_pd(_mm512_unpackhi_pd(v_righterm2, v_righterm3),
                             _mm512_unpacklo_pd(v_righterm2, v_righterm3));
        xmm2 = _mm512_add_pd(_mm512_unpackhi_pd(v_righterm4, v_righterm5),
                             _mm512_unpacklo_pd(v_righterm4, v_righterm5));
        xmm3 = _mm512_add_pd(_mm512_unpackhi_pd(v_righterm6, v_righterm7),
                             _mm512_unpacklo_pd(v_righterm6, v_righterm7));

        ymm0 = _mm512_add_pd(_mm512_permutex2var_pd(xmm0, permute_mask, xmm2),
                             _mm512_mask_blend_pd(0xF0, xmm0, xmm2));

        ymm1 = _mm512_add_pd(_mm512_permutex2var_pd(xmm1, permute_mask, xmm3),
                             _mm512_mask_blend_pd(0xF0, xmm1, xmm3));

        __m512d v_righterm_sum = _mm512_add_pd(_mm512_permutex2var_pd(ymm0,
                                                                      permute_mask_final_stage,
                                                                      ymm1),
                                               _mm512_mask_blend_pd(0xCC, ymm0, ymm1));

        /* update sum */
        __m512d v_prod = _mm512_mul_pd(v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0) {
          v_prod = _mm512_mul_pd(v_prod, v_scale_minlh[rate_scalings[i] - 1]);
        }

        _mm512_store_pd(sum + j, v_prod);
      }

      t_clvc += states_padded;
      t_clvp += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free(tt_inv_eigenvecs);
  pll_aligned_free(tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti_avx512f(unsigned int states,
                                                   unsigned int sites,
                                                   unsigned int rate_cats,
                                                   const double *parent_clv,
                                                   const unsigned char *left_tipchars,
                                                   const unsigned int *parent_scaler,
                                                   double *const *eigenvecs,
                                                   double *const *inv_eigenvecs,
                                                   double *const *freqs,
                                                   const pll_state_t *tipmap,
                                                   unsigned int tipmap_size,
                                                   double *sumtable,
                                                   unsigned int attrib) {
  assert(0);
  //TODO: Not implemented!
  return PLL_FAILURE;
}

PLL_EXPORT
int pll_core_likelihood_derivatives_avx512f(unsigned int states,
                                            unsigned int states_padded,
                                            unsigned int rate_cats,
                                            unsigned int ef_sites,
                                            const unsigned int *pattern_weights,
                                            const double *rate_weights,
                                            const int *invariant,
                                            const double *prop_invar,
                                            double *const *freqs,
                                            const double *sumtable,
                                            const double *diagptable,
                                            double *d_f,
                                            double *dd_f) {
  /* vectors for accumulating LH, 1st and 2nd derivatives */
  __m512d v_df = _mm512_setzero_pd();
  __m512d v_ddf = _mm512_setzero_pd();
  __m512d v_all1 = _mm512_set1_pd(1.);

  __m512d site_lk[3];

  const double *sum = sumtable;
  const int *invariant_ptr = invariant;

  double *t_diagp = (double *) pll_aligned_alloc(
          ELEM_PER_AVX515_REGISTER * 3 * rate_cats * states * sizeof(double), PLL_ALIGNMENT_AVX512F);

  if (!t_diagp) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* transpose diagptable */
  for (unsigned int i = 0; i < rate_cats; ++i) {
    for (unsigned int j = 0; j < states; ++j) {
      for (unsigned int k = 0; k < 3; ++k) {
        __m512d v_diagp = _mm512_set1_pd(diagptable[i * states * 4 + j * 4 + k]);
        _mm512_store_pd(t_diagp +
                        i * ELEM_PER_AVX515_REGISTER * 3 * states +
                        j * ELEM_PER_AVX515_REGISTER * 3 +
                        k * ELEM_PER_AVX515_REGISTER,
                        v_diagp);
      }
    }
  }

  for (unsigned int n = 0;
       n < ef_sites;
       n += ELEM_PER_AVX515_REGISTER, invariant_ptr += ELEM_PER_AVX515_REGISTER) {
    site_lk[0] = _mm512_setzero_pd();
    site_lk[1] = _mm512_setzero_pd();
    site_lk[2] = _mm512_setzero_pd();

    const double *diagp = t_diagp;

    for (unsigned int i = 0; i < rate_cats; ++i) {

      __m512d v_cat_sitelk[3];
      v_cat_sitelk[0] = _mm512_setzero_pd();
      v_cat_sitelk[1] = _mm512_setzero_pd();
      v_cat_sitelk[2] = _mm512_setzero_pd();

      for (unsigned int j = 0;
           j < states; j++, diagp += 3 * ELEM_PER_AVX515_REGISTER, sum += ELEM_PER_AVX515_REGISTER) {
        __m512d v_sum = _mm512_load_pd(sum);
        __m512d v_diagp;

        v_diagp = _mm512_load_pd(diagp);
        //v_diagp = _mm512_set1_pd(diagp[0]);
        v_cat_sitelk[0] = _mm512_fmadd_pd(v_sum, v_diagp, v_cat_sitelk[0]);

        v_diagp = _mm512_load_pd(diagp + ELEM_PER_AVX515_REGISTER);
        v_cat_sitelk[1] = _mm512_fmadd_pd(v_sum, v_diagp, v_cat_sitelk[1]);

        v_diagp = _mm512_load_pd(diagp + 2 * ELEM_PER_AVX515_REGISTER);
        v_cat_sitelk[2] = _mm512_fmadd_pd(v_sum, v_diagp, v_cat_sitelk[2]);
      }

      /* account for invariant sites */
      double t_prop_invar = prop_invar[i];
      if (t_prop_invar > 0) {

        //TODO Vectorize?
        double inv_site_lk_0 =
                (n + 0 >= ef_sites || invariant_ptr[0] == -1) ? 0 : freqs[i][invariant_ptr[0]] * t_prop_invar;
        double inv_site_lk_1 =
                (n + 1 >= ef_sites || invariant_ptr[1] == -1) ? 0 : freqs[i][invariant_ptr[1]] * t_prop_invar;
        double inv_site_lk_2 =
                (n + 2 >= ef_sites || invariant_ptr[2] == -1) ? 0 : freqs[i][invariant_ptr[2]] * t_prop_invar;
        double inv_site_lk_3 =
                (n + 3 >= ef_sites || invariant_ptr[3] == -1) ? 0 : freqs[i][invariant_ptr[3]] * t_prop_invar;
        double inv_site_lk_4 =
                (n + 4 >= ef_sites || invariant_ptr[4] == -1) ? 0 : freqs[i][invariant_ptr[4]] * t_prop_invar;
        double inv_site_lk_5 =
                (n + 5 >= ef_sites || invariant_ptr[5] == -1) ? 0 : freqs[i][invariant_ptr[5]] * t_prop_invar;
        double inv_site_lk_6 =
                (n + 6 >= ef_sites || invariant_ptr[6] == -1) ? 0 : freqs[i][invariant_ptr[6]] * t_prop_invar;
        double inv_site_lk_7 =
                (n + 7 >= ef_sites || invariant_ptr[7] == -1) ? 0 : freqs[i][invariant_ptr[7]] * t_prop_invar;

        __m512d v_inv_site_lk = _mm512_setr_pd(inv_site_lk_0,
                                               inv_site_lk_1,
                                               inv_site_lk_2,
                                               inv_site_lk_3,
                                               inv_site_lk_4,
                                               inv_site_lk_5,
                                               inv_site_lk_6,
                                               inv_site_lk_7);

        __m512d v_prop_invar = _mm512_set1_pd(1. - t_prop_invar);

        v_cat_sitelk[0] = _mm512_add_pd(_mm512_mul_pd(v_cat_sitelk[0], v_prop_invar), v_inv_site_lk);
        v_cat_sitelk[1] = _mm512_mul_pd(v_cat_sitelk[1], v_prop_invar);
        v_cat_sitelk[2] = _mm512_mul_pd(v_cat_sitelk[2], v_prop_invar);
      }

      /* apply rate category weights */
      __m512d v_weight = _mm512_set1_pd(rate_weights[i]);
      site_lk[0] = _mm512_fmadd_pd(v_cat_sitelk[0], v_weight, site_lk[0]);
      site_lk[1] = _mm512_fmadd_pd(v_cat_sitelk[1], v_weight, site_lk[1]);
      site_lk[2] = _mm512_fmadd_pd(v_cat_sitelk[2], v_weight, site_lk[2]);
    }

    /* build derivatives */
    __m512d v_recip0 = _mm512_div_pd(v_all1, site_lk[0]);
    __m512d v_deriv1 = _mm512_mul_pd(site_lk[1], v_recip0);
    __m512d v_deriv2 = _mm512_sub_pd(_mm512_mul_pd(v_deriv1, v_deriv1),
                                     _mm512_mul_pd(site_lk[2], v_recip0));

    /* eliminates nan values on padded states */
    if (n + ELEM_PER_AVX515_REGISTER > ef_sites) {
      __mmask8 mask = _mm512_cmp_pd_mask(site_lk[0], _mm512_setzero_pd(), _CMP_NEQ_UQ);

      v_deriv1 = _mm512_maskz_expand_pd(mask, v_deriv1);
      v_deriv2 = _mm512_maskz_expand_pd(mask, v_deriv2);
    }

    v_df = _mm512_fnmadd_pd(v_deriv1, _mm512_set1_pd(pattern_weights[n]), v_df);
    v_ddf = _mm512_fmadd_pd(v_deriv2, _mm512_set1_pd(pattern_weights[n]), v_ddf);
  }

  *d_f = reduce_add_pd(v_df);
  *dd_f = reduce_add_pd(v_ddf);

  return PLL_SUCCESS;
}
