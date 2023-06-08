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

#define COMPUTE_II_QCOL(q, offset) \
/* row 0 */ \
v_mat    = _mm256_load_pd(lm0 + offset); \
v_lterm0 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_lterm0); \
v_mat    = _mm256_load_pd(rm0 + offset); \
v_rterm0 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm0); \
 \
/* row 1 */ \
v_mat    = _mm256_load_pd(lm1 + offset); \
v_lterm1 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_lterm1); \
v_mat    = _mm256_load_pd(rm1 + offset); \
v_rterm1 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm1); \
\
/* row 2 */ \
v_mat    = _mm256_load_pd(lm2 + offset); \
v_lterm2 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_lterm2); \
v_mat    = _mm256_load_pd(rm2 + offset); \
v_rterm2 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm2); \
\
/* row 3 */ \
v_mat    = _mm256_load_pd(lm3 + offset); \
v_lterm3 = _mm256_fmadd_pd(v_mat, v_lclv[q], v_lterm3); \
v_mat    = _mm256_load_pd(rm3 + offset); \
v_rterm3 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm3);


PLL_EXPORT int pll_core_update_sumtable_ii_20x20_avx2(unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * clvp,
                                                const double * clvc,
                                                const unsigned int * parent_scaler,
                                                const unsigned int * child_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                double * sumtable,
                                                unsigned int attrib)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  const double * t_lclv = clvp;
  const double * t_rclv = clvc;
  double * t_freqs;

  unsigned int states = 20;
  unsigned int states_padded = states;

  double * tt_inv_eigenvecs = NULL;
  double * tt_eigenvecs = NULL;

  /* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_AVX);

  /* transposed padded inv_eigenvecs */
  tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    if (tt_eigenvecs)
      pll_aligned_free(tt_eigenvecs);
    if (tt_inv_eigenvecs)
      pll_aligned_free(tt_inv_eigenvecs);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    const double * c_eigenvecs      = tt_eigenvecs;
    const double * ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {
      __m256d v_lclv[5];
      __m256d v_rclv[5];
      for (j = 0; j < 5; ++j)
      {
        v_lclv[j]    = _mm256_load_pd(t_lclv + j*4);
        v_rclv[j]    = _mm256_load_pd(t_rclv + j*4);
      }

      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvecs matrix */
        const double * lm0 = ct_inv_eigenvecs;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;
        ct_inv_eigenvecs += 4*states_padded;

        /* point to the four rows of the inv_eigenvecs matrix */
        const double * rm0 = c_eigenvecs;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;
        c_eigenvecs += 4*states_padded;

        __m256d v_lterm0 = _mm256_setzero_pd ();
        __m256d v_rterm0 = _mm256_setzero_pd ();
        __m256d v_lterm1 = _mm256_setzero_pd ();
        __m256d v_rterm1 = _mm256_setzero_pd ();
        __m256d v_lterm2 = _mm256_setzero_pd ();
        __m256d v_rterm2 = _mm256_setzero_pd ();
        __m256d v_lterm3 = _mm256_setzero_pd ();
        __m256d v_rterm3 = _mm256_setzero_pd ();

        __m256d v_mat;

        /* iterate over quadruples of columns */
        COMPUTE_II_QCOL(0, 0);
        COMPUTE_II_QCOL(1, 4);
        COMPUTE_II_QCOL(2, 8);
        COMPUTE_II_QCOL(3, 12);
        COMPUTE_II_QCOL(4, 16);

        /* compute lefterm */
        __m256d xmm0 = _mm256_unpackhi_pd (v_lterm0, v_lterm1);
        __m256d xmm1 = _mm256_unpacklo_pd (v_lterm0, v_lterm1);
        __m256d xmm2 = _mm256_unpackhi_pd (v_lterm2, v_lterm3);
        __m256d xmm3 = _mm256_unpacklo_pd (v_lterm2, v_lterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* compute righterm */
        xmm0 = _mm256_unpackhi_pd (v_rterm0, v_rterm1);
        xmm1 = _mm256_unpacklo_pd (v_rterm0, v_rterm1);
        xmm2 = _mm256_unpackhi_pd (v_rterm2, v_rterm3);
        xmm3 = _mm256_unpacklo_pd (v_rterm2, v_rterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* update sum */
        __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm256_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd (sum + j, v_prod);
      }

      t_lclv += states_padded;
      t_rclv += states_padded;
      sum    += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * clvp,
                                                const double * clvc,
                                                const unsigned int * parent_scaler,
                                                const unsigned int * child_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                double * sumtable,
                                                unsigned int attrib)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  const double * t_clvp = clvp;
  const double * t_clvc = clvc;
  double * t_freqs;

  /* dedicated functions for 4x4 and 20x20 matrices */
  if (states == 4)
  {
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
  }
  else if (states == 20)
  {
    /* call AVX variant */
    return pll_core_update_sumtable_ii_20x20_avx2(sites,
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


  unsigned int states_padded = (states+3) & 0xFFFFFFFC;

  double * tt_inv_eigenvecs = NULL;
  double * tt_eigenvecs = NULL;

  /* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_AVX);

  /* transposed padded inv_eigenvecs */
  tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    if (tt_eigenvecs)
      pll_aligned_free(tt_eigenvecs);
    if (tt_inv_eigenvecs)
      pll_aligned_free(tt_inv_eigenvecs);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    const double * c_eigenvecs      = tt_eigenvecs;
    const double * ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {
      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvecs matrix */
        const double * em0 = c_eigenvecs;
        const double * em1 = em0 + states_padded;
        const double * em2 = em1 + states_padded;
        const double * em3 = em2 + states_padded;
        c_eigenvecs += 4*states_padded;

        /* point to the four rows of the inv_eigenvecs matrix */
        const double * im0 = ct_inv_eigenvecs;
        const double * im1 = im0 + states_padded;
        const double * im2 = im1 + states_padded;
        const double * im3 = im2 + states_padded;
        ct_inv_eigenvecs += 4*states_padded;

        __m256d v_lefterm0 = _mm256_setzero_pd ();
        __m256d v_righterm0 = _mm256_setzero_pd ();
        __m256d v_lefterm1 = _mm256_setzero_pd ();
        __m256d v_righterm1 = _mm256_setzero_pd ();
        __m256d v_lefterm2 = _mm256_setzero_pd ();
        __m256d v_righterm2 = _mm256_setzero_pd ();
        __m256d v_lefterm3 = _mm256_setzero_pd ();
        __m256d v_righterm3 = _mm256_setzero_pd ();

        __m256d v_eigen;
        __m256d v_clvp;
        __m256d v_clvc;

        for (k = 0; k < states_padded; k += 4)
        {
          v_clvp = _mm256_load_pd (t_clvp + k);
          v_clvc = _mm256_load_pd (t_clvc + k);

          /* row 0 */
          v_eigen = _mm256_load_pd (im0 + k);
          v_lefterm0 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm0);


          v_eigen = _mm256_load_pd (em0 + k);
          v_righterm0 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm0);

          /* row 1 */
          v_eigen = _mm256_load_pd (im1 + k);
          v_lefterm1 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm1);

          v_eigen = _mm256_load_pd (em1 + k);
          v_righterm1 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm1);


          /* row 2 */
          v_eigen = _mm256_load_pd (im2 + k);
          v_lefterm2 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm2);

          v_eigen = _mm256_load_pd (em2 + k);
          v_righterm2 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm2);

          /* row 3 */
          v_eigen = _mm256_load_pd (im3 + k);
          v_lefterm3 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm3);

          v_eigen = _mm256_load_pd (em3 + k);
          v_righterm3 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm3);

        }

        /* compute lefterm */
        __m256d xmm0 = _mm256_unpackhi_pd (v_lefterm0, v_lefterm1);
        __m256d xmm1 = _mm256_unpacklo_pd (v_lefterm0, v_lefterm1);
        __m256d xmm2 = _mm256_unpackhi_pd (v_lefterm2, v_lefterm3);
        __m256d xmm3 = _mm256_unpacklo_pd (v_lefterm2, v_lefterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* compute righterm */
        xmm0 = _mm256_unpackhi_pd (v_righterm0, v_righterm1);
        xmm1 = _mm256_unpacklo_pd (v_righterm0, v_righterm1);
        xmm2 = _mm256_unpackhi_pd (v_righterm2, v_righterm3);
        xmm3 = _mm256_unpacklo_pd (v_righterm2, v_righterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* update sum */
        __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm256_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd (sum + j, v_prod);
      }

      t_clvc += states_padded;
      t_clvp += states_padded;
      sum    += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_repeats_20x20_avx2(unsigned int sites,
                                                           unsigned int parent_sites,
                                                           unsigned int rate_cats,
                                                           const double * clvp,
                                                           const double * clvc,
                                                           const unsigned int * parent_scaler,
                                                           const unsigned int * child_scaler,
                                                           double * const * eigenvecs,
                                                           double * const * inv_eigenvecs,
                                                           double * const * freqs,
                                                           double *sumtable,
                                                           const unsigned int * parent_site_id,
                                                           const unsigned int * child_site_id,
                                                           double * bclv_buffer,
                                                           unsigned int inv,
                                                           unsigned int attrib)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  double * t_freqs;

  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int span_padded = rate_cats * states_padded;

  double * tt_inv_eigenvecs = NULL;
  double * tt_eigenvecs = NULL;

  /* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_AVX);

  /* transposed padded inv_eigenvecs */
  tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    if (tt_eigenvecs)
      pll_aligned_free(tt_eigenvecs);
    if (tt_inv_eigenvecs)
      pll_aligned_free(tt_inv_eigenvecs);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    unsigned int pid = PLL_GET_ID(parent_site_id, n);
    unsigned int cid = PLL_GET_ID(child_site_id, n);
    const double * t_lclv = &clvp[pid * span_padded];
    const double * t_rclv = &clvc[cid * span_padded];
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[pid*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[cid*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    const double * c_eigenvecs      = tt_eigenvecs;
    const double * ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {

      __m256d v_lclv[5];
      __m256d v_rclv[5];
      for (j = 0; j < 5; ++j)
      {
        v_lclv[j]    = _mm256_load_pd(t_lclv + j*4);
        v_rclv[j]    = _mm256_load_pd(t_rclv + j*4);
      }

      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvecs matrix */
        const double * lm0 = ct_inv_eigenvecs;
        const double * lm1 = lm0 + states_padded;
        const double * lm2 = lm1 + states_padded;
        const double * lm3 = lm2 + states_padded;
        ct_inv_eigenvecs += 4*states_padded;

        /* point to the four rows of the inv_eigenvecs matrix */
        const double * rm0 = c_eigenvecs;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;
        c_eigenvecs += 4*states_padded;

        __m256d v_lterm0 = _mm256_setzero_pd ();
        __m256d v_rterm0 = _mm256_setzero_pd ();
        __m256d v_lterm1 = _mm256_setzero_pd ();
        __m256d v_rterm1 = _mm256_setzero_pd ();
        __m256d v_lterm2 = _mm256_setzero_pd ();
        __m256d v_rterm2 = _mm256_setzero_pd ();
        __m256d v_lterm3 = _mm256_setzero_pd ();
        __m256d v_rterm3 = _mm256_setzero_pd ();

        __m256d v_mat;

        /* iterate over quadruples of columns */
        COMPUTE_II_QCOL(0, 0);
        COMPUTE_II_QCOL(1, 4);
        COMPUTE_II_QCOL(2, 8);
        COMPUTE_II_QCOL(3, 12);
        COMPUTE_II_QCOL(4, 16);

        /* compute lefterm */
        __m256d xmm0 = _mm256_unpackhi_pd (v_lterm0, v_lterm1);
        __m256d xmm1 = _mm256_unpacklo_pd (v_lterm0, v_lterm1);
        __m256d xmm2 = _mm256_unpackhi_pd (v_lterm2, v_lterm3);
        __m256d xmm3 = _mm256_unpacklo_pd (v_lterm2, v_lterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* compute righterm */
        xmm0 = _mm256_unpackhi_pd (v_rterm0, v_rterm1);
        xmm1 = _mm256_unpacklo_pd (v_rterm0, v_rterm1);
        xmm2 = _mm256_unpackhi_pd (v_rterm2, v_rterm3);
        xmm3 = _mm256_unpacklo_pd (v_rterm2, v_rterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* update sum */
        __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm256_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd (sum + j, v_prod);
      }

      t_lclv += states_padded;
      t_rclv += states_padded;
      sum    += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}


PLL_EXPORT int pll_core_update_sumtable_repeats_generic_avx2(unsigned int states,
                                                             unsigned int sites,
                                                             unsigned int parent_sites,
                                                             unsigned int rate_cats,
                                                             const double * clvp,
                                                             const double * clvc,
                                                             const unsigned int * parent_scaler,
                                                             const unsigned int * child_scaler,
                                                             double * const * eigenvecs,
                                                             double * const * inv_eigenvecs,
                                                             double * const * freqs,
                                                             double *sumtable,
                                                             const unsigned int * parent_site_id,
                                                             const unsigned int * child_site_id,
                                                             double * bclv_buffer,
                                                             unsigned int inv,
                                                             unsigned int attrib)
{
  unsigned int i, j, k, n;

  /* build sumtable */
  double * sum = sumtable;

  double * t_freqs;

  if (states == 20)
  {
    return pll_core_update_sumtable_repeats_20x20_avx2(sites,
                                                       parent_sites,
                                                       rate_cats,
                                                       clvp,
                                                       clvc,
                                                       parent_scaler,
                                                       child_scaler,
                                                       eigenvecs,
                                                       inv_eigenvecs,
                                                       freqs,
                                                       sumtable,
                                                       parent_site_id,
                                                       child_site_id,
                                                       bclv_buffer,
                                                       inv,
                                                       attrib);
  }

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span_padded = rate_cats * states_padded;

  double * tt_eigenvecs = NULL;
  double * tt_inv_eigenvecs = NULL;

  /* scaling stuff */
  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_AVX);

  /* transposed padded inv_eigenvecs */
  tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    if (tt_eigenvecs)
      pll_aligned_free(tt_eigenvecs);
    if (tt_inv_eigenvecs)
      pll_aligned_free(tt_inv_eigenvecs);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* vectorized loop from update_sumtable() */
  for (n = 0; n < sites; n++)
  {
    unsigned int pid = PLL_GET_ID(parent_site_id, n);
    unsigned int cid = PLL_GET_ID(child_site_id, n);
    const double * t_clvp = &clvp[pid * span_padded];
    const double * t_clvc = &clvc[cid * span_padded];
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[pid*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[cid*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    const double * c_eigenvecs      = tt_eigenvecs;
    const double * ct_inv_eigenvecs = tt_inv_eigenvecs;
    for (i = 0; i < rate_cats; ++i)
    {
      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvecs matrix */
        const double * em0 = c_eigenvecs;
        const double * em1 = em0 + states_padded;
        const double * em2 = em1 + states_padded;
        const double * em3 = em2 + states_padded;
        c_eigenvecs += 4*states_padded;

        /* point to the four rows of the inv_eigenvecs matrix */
        const double * im0 = ct_inv_eigenvecs;
        const double * im1 = im0 + states_padded;
        const double * im2 = im1 + states_padded;
        const double * im3 = im2 + states_padded;
        ct_inv_eigenvecs += 4*states_padded;

        __m256d v_lefterm0 = _mm256_setzero_pd ();
        __m256d v_righterm0 = _mm256_setzero_pd ();
        __m256d v_lefterm1 = _mm256_setzero_pd ();
        __m256d v_righterm1 = _mm256_setzero_pd ();
        __m256d v_lefterm2 = _mm256_setzero_pd ();
        __m256d v_righterm2 = _mm256_setzero_pd ();
        __m256d v_lefterm3 = _mm256_setzero_pd ();
        __m256d v_righterm3 = _mm256_setzero_pd ();

        __m256d v_eigen;
        __m256d v_clvp;
        __m256d v_clvc;

        for (k = 0; k < states_padded; k += 4)
        {
          v_clvp = _mm256_load_pd (t_clvp + k);
          v_clvc = _mm256_load_pd (t_clvc + k);

          /* row 0 */
          v_eigen = _mm256_load_pd (im0 + k);
          v_lefterm0 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm0);


          v_eigen = _mm256_load_pd (em0 + k);
          v_righterm0 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm0);

          /* row 1 */
          v_eigen = _mm256_load_pd (im1 + k);
          v_lefterm1 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm1);

          v_eigen = _mm256_load_pd (em1 + k);
          v_righterm1 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm1);


          /* row 2 */
          v_eigen = _mm256_load_pd (im2 + k);
          v_lefterm2 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm2);

          v_eigen = _mm256_load_pd (em2 + k);
          v_righterm2 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm2);

          /* row 3 */
          v_eigen = _mm256_load_pd (im3 + k);
          v_lefterm3 = _mm256_fmadd_pd(v_eigen, v_clvp, v_lefterm3);

          v_eigen = _mm256_load_pd (em3 + k);
          v_righterm3 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm3);

        }

        /* compute lefterm */
        __m256d xmm0 = _mm256_unpackhi_pd (v_lefterm0, v_lefterm1);
        __m256d xmm1 = _mm256_unpacklo_pd (v_lefterm0, v_lefterm1);
        __m256d xmm2 = _mm256_unpackhi_pd (v_lefterm2, v_lefterm3);
        __m256d xmm3 = _mm256_unpacklo_pd (v_lefterm2, v_lefterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_lefterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* compute righterm */
        xmm0 = _mm256_unpackhi_pd (v_righterm0, v_righterm1);
        xmm1 = _mm256_unpacklo_pd (v_righterm0, v_righterm1);
        xmm2 = _mm256_unpackhi_pd (v_righterm2, v_righterm3);
        xmm3 = _mm256_unpacklo_pd (v_righterm2, v_righterm3);
        xmm0 = _mm256_add_pd (xmm0, xmm1);
        xmm1 = _mm256_add_pd (xmm2, xmm3);
        xmm2 = _mm256_permute2f128_pd (xmm0, xmm1, _MM_SHUFFLE(0, 2, 0, 1));
        xmm3 = _mm256_blend_pd (xmm0, xmm1, 12);
        __m256d v_righterm_sum = _mm256_add_pd (xmm2, xmm3);

        /* update sum */
        __m256d v_prod = _mm256_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm256_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd (sum + j, v_prod);
      }

      t_clvc += states_padded;
      t_clvp += states_padded;
      sum    += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

#define COMPUTE_TI_QCOL(q, offset) \
/* row 0 */ \
v_mat    = _mm256_load_pd(rm0 + offset); \
v_rterm0 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm0); \
 \
/* row 1 */ \
v_mat    = _mm256_load_pd(rm1 + offset); \
v_rterm1 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm1); \
\
/* row 2 */ \
v_mat    = _mm256_load_pd(rm2 + offset); \
v_rterm2 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm2); \
\
/* row 3 */ \
v_mat    = _mm256_load_pd(rm3 + offset); \
v_rterm3 = _mm256_fmadd_pd(v_mat, v_rclv[q], v_rterm3);

PLL_EXPORT int pll_core_update_sumtable_ti_20x20_avx2(unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * parent_clv,
                                                const unsigned char * left_tipchars,
                                                const unsigned int * parent_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                const pll_state_t * tipmap,
                                                unsigned int tipmap_size,
                                                double *sumtable,
                                                unsigned int attrib)
{
  unsigned int states = 20;
  unsigned int states_padded = states;
  unsigned int span = states_padded * rate_cats;
  unsigned int maxstates = tipmap_size;

  unsigned int i, j, k, n;
  unsigned int tipstate;

  double * eigenvecs_padded = NULL;
  double * precomp_left = NULL;

  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  double * sum = sumtable;
  const double * t_rclv = parent_clv;
  const double * t_eigenvecs_padded;

  eigenvecs_padded = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  precomp_left = (double *) pll_aligned_alloc (
      (maxstates * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!eigenvecs_padded || !precomp_left)
  {
    if (eigenvecs_padded)
      pll_aligned_free(eigenvecs_padded);
    if (precomp_left)
      pll_aligned_free(precomp_left);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  /* add padding to eigenvecs matrix -> for efficient vectorization */
  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states_padded; ++j)
      for (k = 0; k < states_padded; ++k)
      {
        eigenvecs_padded[i*states_padded*states_padded + j*states_padded + k] =
            (j < states && k < states) ? eigenvecs[i][j*states_padded + k] : 0.;
      }
  }

  /* precompute left terms since they are the same for every site */
  double * t_precomp = precomp_left;
  for (n = 0; n < maxstates; ++n)
  {
    unsigned int state = tipmap ? (unsigned int) tipmap[n] : n;

    int ss = PLL_POPCNT32(state) == 1 ? PLL_CTZ32(state) : -1;

    for (i = 0; i < rate_cats; ++i)
    {
      for (j = 0; j < states_padded; j += 4)
      {
        __m256d v_lefterm;

        if (ss != -1)
        {
          /* special case for non-ambiguous state */
          __m256d v_freqs = _mm256_set1_pd(freqs[i][ss]);
          __m256d v_eigen = _mm256_load_pd(inv_eigenvecs[i] +
                                                       ss*states_padded + j);
          v_lefterm =  _mm256_mul_pd(v_eigen, v_freqs);
        }
        else
        {
          v_lefterm = _mm256_setzero_pd();
          for (k = 0; k < states; ++k)
          {
            if ((state>>k) & 1)
            {
              __m256d v_freqs = _mm256_set1_pd(freqs[i][k]);
              __m256d v_eigen = _mm256_load_pd(inv_eigenvecs[i] +
                                                           k*states_padded + j);

              v_lefterm = _mm256_fmadd_pd(v_eigen, v_freqs, v_lefterm);
            }
          }
        }

        _mm256_store_pd(t_precomp, v_lefterm);
        t_precomp += 4;
      }
    }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    tipstate = (unsigned int) left_tipchars[n];

    unsigned int loffset = tipstate * span;

    t_eigenvecs_padded = eigenvecs_padded;
    t_precomp = precomp_left + loffset;

    for (i = 0; i < rate_cats; ++i)
    {
      __m256d v_rclv[5];
      for (j = 0; j < 5; ++j)
      {
        v_rclv[j]    = _mm256_load_pd(t_rclv + j*4);
      }

      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvec matrix */
        const double * rm0 = t_eigenvecs_padded;
        const double * rm1 = rm0 + states_padded;
        const double * rm2 = rm1 + states_padded;
        const double * rm3 = rm2 + states_padded;
        t_eigenvecs_padded += 4*states_padded;

        __m256d v_rterm0 = _mm256_setzero_pd();
        __m256d v_rterm1 = _mm256_setzero_pd();
        __m256d v_rterm2 = _mm256_setzero_pd();
        __m256d v_rterm3 = _mm256_setzero_pd();

        __m256d v_mat;

        COMPUTE_TI_QCOL(0, 0);
        COMPUTE_TI_QCOL(1, 4);
        COMPUTE_TI_QCOL(2, 8);
        COMPUTE_TI_QCOL(3, 12);
        COMPUTE_TI_QCOL(4, 16);

        /* reduce righterm */
        __m256d xmm0, xmm1, xmm2, xmm3;
        xmm0 = _mm256_unpackhi_pd(v_rterm0,v_rterm1);
        xmm1 = _mm256_unpacklo_pd(v_rterm0,v_rterm1);

        xmm2 = _mm256_unpackhi_pd(v_rterm2,v_rterm3);
        xmm3 = _mm256_unpacklo_pd(v_rterm2,v_rterm3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_righterm = _mm256_add_pd(xmm2,xmm3);

        __m256d v_lefterm = _mm256_load_pd(t_precomp + j);

        __m256d v_sum = _mm256_mul_pd(v_lefterm, v_righterm);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_sum = _mm256_mul_pd(v_sum, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd(sum + j, v_sum);
      }

      t_rclv += states_padded;
      t_precomp += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free(eigenvecs_padded);
  pll_aligned_free(precomp_left);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti_avx2(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                const double * parent_clv,
                                                const unsigned char * left_tipchars,
                                                const unsigned int * parent_scaler,
                                                double * const * eigenvecs,
                                                double * const * inv_eigenvecs,
                                                double * const * freqs,
                                                const pll_state_t * tipmap,
                                                unsigned int tipmap_size,
                                                double *sumtable,
                                                unsigned int attrib)
{
  if (states == 4)
  {
    /* call AVX version for the 4x4 case */
    return pll_core_update_sumtable_ti_avx(states,
                                           sites,
                                           rate_cats,
                                           parent_clv,
                                           left_tipchars,
                                           parent_scaler,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs,
                                           tipmap,
                                           tipmap_size,
                                           sumtable,
                                           attrib);
  }
  else if (states == 20)
  {
    return pll_core_update_sumtable_ti_20x20_avx2(sites,
                                           rate_cats,
                                           parent_clv,
                                           left_tipchars,
                                           parent_scaler,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs,
                                           tipmap,
                                           tipmap_size,
                                           sumtable,
                                           attrib);
  }

  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int span = states_padded * rate_cats;
  unsigned int maxstates = tipmap_size;

  unsigned int i, j, k, n;
  pll_state_t tipstate;

  double * eigenvecs_padded = NULL;
  double * precomp_left = NULL;

  unsigned int min_scaler = 0;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m256d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200, "Cannot allocate memory for rate scalers");
      return PLL_FAILURE;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm256_set1_pd(scale_factor);
    }
  }

  double * sum = sumtable;
  const double * t_clvc = parent_clv;
  const double * t_eigenvecs_padded;

  eigenvecs_padded = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  precomp_left = (double *) pll_aligned_alloc (
      (maxstates * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_AVX);

  if (!eigenvecs_padded || !precomp_left)
  {
    if (eigenvecs_padded)
      pll_aligned_free(eigenvecs_padded);
    if (precomp_left)
      pll_aligned_free(precomp_left);
    if (rate_scalings)
      free(rate_scalings);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  /* add padding to eigenvecs matrix -> for efficient vectorization */
  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states_padded; ++j)
      for (k = 0; k < states_padded; ++k)
      {
        eigenvecs_padded[i*states_padded*states_padded + j*states_padded + k] =
            (j < states && k < states) ? eigenvecs[i][j*states_padded + k] : 0.;
      }
  }

  /* precompute left terms since they are the same for every site */
  double * t_precomp = precomp_left;
  for (n = 0; n < maxstates; ++n)
  {
    pll_state_t state = tipmap ? tipmap[n] : n;

    int ss = PLL_STATE_POPCNT(state) == 1 ? PLL_STATE_CTZ(state) : -1;

    for (i = 0; i < rate_cats; ++i)
    {
      for (j = 0; j < states_padded; j += 4)
      {
        __m256d v_lefterm;

        if (ss != -1)
        {
          /* special case for non-ambiguous state */
          __m256d v_freqs = _mm256_set1_pd(freqs[i][ss]);
          __m256d v_eigen = _mm256_load_pd(inv_eigenvecs[i] +
                                                       ss*states_padded + j);
          v_lefterm =  _mm256_mul_pd(v_eigen, v_freqs);
        }
        else
        {
          v_lefterm = _mm256_setzero_pd();
          for (k = 0; k < states; ++k)
          {
            if ((state>>k) & 1)
            {
              __m256d v_freqs = _mm256_set1_pd(freqs[i][k]);
              __m256d v_eigen = _mm256_load_pd(inv_eigenvecs[i] +
                                                           k*states_padded + j);

              v_lefterm = _mm256_fmadd_pd(v_eigen, v_freqs, v_lefterm);
            }
          }
        }

        _mm256_store_pd(t_precomp, v_lefterm);
        t_precomp += 4;
      }
    }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    /* compute per-rate scalers and obtain minimum value (within site) */
    if (per_rate_scaling)
    {
      min_scaler = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < min_scaler)
          min_scaler = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - min_scaler,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }

    tipstate = (unsigned int) left_tipchars[n];

    unsigned int loffset = tipstate * span;

    t_eigenvecs_padded = eigenvecs_padded;
    t_precomp = precomp_left + loffset;

    for (i = 0; i < rate_cats; ++i)
    {
      for (j = 0; j < states_padded; j += 4)
      {
        /* point to the four rows of the eigenvec matrix */
        const double * em0 = t_eigenvecs_padded;
        const double * em1 = em0 + states_padded;
        const double * em2 = em1 + states_padded;
        const double * em3 = em2 + states_padded;
        t_eigenvecs_padded += 4*states_padded;

        __m256d v_righterm0 = _mm256_setzero_pd();
        __m256d v_righterm1 = _mm256_setzero_pd();
        __m256d v_righterm2 = _mm256_setzero_pd();
        __m256d v_righterm3 = _mm256_setzero_pd();

        for (k = 0; k < states_padded; k += 4)
        {
          /* load 4 entries of CLV */
          __m256d v_clvc = _mm256_load_pd(t_clvc + k);

          /* row 0 */
          __m256d v_eigen = _mm256_load_pd(em0 + k);
          v_righterm0 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm0);

          /* row 1 */
          v_eigen = _mm256_load_pd(em1 + k);
          v_righterm1 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm1);

          /* row 2 */
          v_eigen = _mm256_load_pd(em2 + k);
          v_righterm2 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm2);

          /* row 3 */
          v_eigen = _mm256_load_pd(em3 + k);
          v_righterm3 = _mm256_fmadd_pd(v_eigen, v_clvc, v_righterm3);
        }

        /* reduce righterm */
        __m256d xmm0, xmm1, xmm2, xmm3;
        xmm0 = _mm256_unpackhi_pd(v_righterm0,v_righterm1);
        xmm1 = _mm256_unpacklo_pd(v_righterm0,v_righterm1);

        xmm2 = _mm256_unpackhi_pd(v_righterm2,v_righterm3);
        xmm3 = _mm256_unpacklo_pd(v_righterm2,v_righterm3);

        xmm0 = _mm256_add_pd(xmm0,xmm1);
        xmm1 = _mm256_add_pd(xmm2,xmm3);

        xmm2 = _mm256_permute2f128_pd(xmm0,xmm1, _MM_SHUFFLE(0,2,0,1));

        xmm3 = _mm256_blend_pd(xmm0,xmm1,12);

        __m256d v_righterm = _mm256_add_pd(xmm2,xmm3);

        __m256d v_lefterm = _mm256_load_pd(t_precomp + j);

        __m256d v_sum = _mm256_mul_pd(v_lefterm, v_righterm);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_sum = _mm256_mul_pd(v_sum, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm256_store_pd(sum + j, v_sum);
      }

      t_clvc += states_padded;
      t_precomp += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free(eigenvecs_padded);
  pll_aligned_free(precomp_left);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT
int pll_core_likelihood_derivatives_avx2(unsigned int states,
                                         unsigned int states_padded,
                                         unsigned int rate_cats,
                                         unsigned int ef_sites,
                                         const double * pattern_weights,
                                         const double * rate_weights,
                                         const int * invariant,
                                         const double * prop_invar,
                                         double * const * freqs,
                                         const double * sumtable,
                                         const double * diagptable,
                                         double * d_f,
                                         double * dd_f)
{
  unsigned int i,j,k,n;
  unsigned int span_padded = rate_cats * states_padded;

  double * t_diagp = NULL;
  const double * diagp_start = NULL;
  double * invar_lk = NULL;

  /* check for special cases in which we can save some computation later on */
  int use_pinv = 0;
  int eq_weights = 1;
  for (i = 0; i < rate_cats; ++i)
  {
    /* check if proportion of invariant site is used */
    use_pinv |= (prop_invar[i] > 0);

    /* check if rate weights are all equal (e.g. GAMMA) */
    eq_weights &= (rate_weights[i] == rate_weights[0]);
  }

  if (use_pinv)
  {
    invar_lk = (double *) pll_aligned_alloc(rate_cats * states * sizeof(double),
                                            PLL_ALIGNMENT_AVX);

    if (!invar_lk)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return PLL_FAILURE;
    }

    /* pre-compute invariant site likelihoods*/
    for(i = 0; i < states; ++i)
    {
      for(j = 0; j < rate_cats; ++j)
      {
        invar_lk[i * rate_cats + j] = freqs[j][i] * prop_invar[j];
      }
    }
  }

  if (states == 4)
  {
    diagp_start = diagptable;
  }
  else
  {
    t_diagp = (double *) pll_aligned_alloc(
                                      3 * span_padded * sizeof(double),
                                      PLL_ALIGNMENT_AVX);

    if (!t_diagp)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return PLL_FAILURE;
    }

    memset(t_diagp, 0, 3 * span_padded * sizeof(double));

    /* transpose diagptable */
    for(i = 0; i < rate_cats; ++i)
    {
      for(j = 0; j < states; ++j)
      {
        for(k = 0; k < 3; ++k)
        {
          t_diagp[i * states_padded * 3 + k * states_padded  + j] =
              diagptable[i * states * 4 + j * 4 + k];
        }
      }
    }

    diagp_start = t_diagp;
  }

  /* here we will temporary store per-site LH, 1st and 2nd derivatives */
  double site_lk[16] __attribute__( ( aligned ( PLL_ALIGNMENT_AVX ) ) ) ;

  /* vectors for accumulating 1st and 2nd derivatives */
  __m256d v_df = _mm256_setzero_pd ();
  __m256d v_ddf = _mm256_setzero_pd ();
  __m256d v_all1 = _mm256_set1_pd(1.);

  const double *sum = sumtable;
  const int * invariant_ptr = invariant;
  unsigned int offset = 0;
  for (n = 0; n < ef_sites; ++n)
  {
    const double * diagp = diagp_start;

    __m256d v_sitelk = _mm256_setzero_pd ();
    for (i = 0; i < rate_cats; ++i)
    {
      __m256d v_cat_sitelk = _mm256_setzero_pd ();

      if (states == 4)
      {
        /* use unrolled loop */
        __m256d v_diagp = _mm256_load_pd(diagp);
        __m256d v_sum = _mm256_set1_pd(sum[0]);
        v_cat_sitelk = _mm256_fmadd_pd(v_sum, v_diagp, v_cat_sitelk);

        v_diagp = _mm256_load_pd(diagp + 4);
        v_sum = _mm256_set1_pd(sum[1]);
        v_cat_sitelk = _mm256_fmadd_pd(v_sum, v_diagp, v_cat_sitelk);

        v_diagp = _mm256_load_pd(diagp + 8);
        v_sum = _mm256_set1_pd(sum[2]);
        v_cat_sitelk = _mm256_fmadd_pd(v_sum, v_diagp, v_cat_sitelk);

        v_diagp = _mm256_load_pd(diagp + 12);
        v_sum = _mm256_set1_pd(sum[3]);
        v_cat_sitelk = _mm256_fmadd_pd(v_sum, v_diagp, v_cat_sitelk);

        diagp += 16;
        sum += 4;
      }
      else
      {
        /* pointer to 3 "rows" of diagp with values for lk0, lk1 and lk2 */
        const double * r0 = diagp;
        const double * r1 = r0 + states_padded;
        const double * r2 = r1 + states_padded;

        /* unroll 1st iteration to save a couple of adds */
        __m256d v_sum = _mm256_load_pd(sum);
        __m256d v_diagp = _mm256_load_pd(r0);
        __m256d v_lk0 = _mm256_mul_pd(v_sum, v_diagp);

        v_diagp = _mm256_load_pd(r1);
        __m256d v_lk1 = _mm256_mul_pd(v_sum, v_diagp);

        v_diagp = _mm256_load_pd(r2);
        __m256d v_lk2 = _mm256_mul_pd(v_sum, v_diagp);

        /* iterate over remaining states (if any) */
        for (j = 4; j < states_padded; j+= 4)
        {
          v_sum = _mm256_load_pd(sum + j);
          v_diagp = _mm256_load_pd(r0 + j);
          v_lk0 =  _mm256_fmadd_pd (v_sum, v_diagp, v_lk0);

          v_diagp = _mm256_load_pd(r1 + j);
          v_lk1 =  _mm256_fmadd_pd (v_sum, v_diagp, v_lk1);

          v_diagp = _mm256_load_pd(r2 + j);
          v_lk2 =  _mm256_fmadd_pd (v_sum, v_diagp, v_lk2);
        }

        /* reduce lk0 (=LH), lk1 (=1st deriv) and v_lk2 (=2nd deriv) */
        v_lk0 = _mm256_hadd_pd(v_lk0, v_lk0);
        double lk0 = ((double *)&v_lk0)[0] + ((double *)&v_lk0)[2];

        v_lk1 = _mm256_hadd_pd(v_lk1, v_lk1);
        double lk1 = ((double *)&v_lk1)[0] + ((double *)&v_lk1)[2];

        v_lk2 = _mm256_hadd_pd(v_lk2, v_lk2);
        double lk2 = ((double *)&v_lk2)[0] + ((double *)&v_lk2)[2];

        v_cat_sitelk = _mm256_setr_pd (lk0, lk1, lk2, 0.);

        sum += states_padded;
        diagp += 3 * states_padded;
      }

      /* account for invariant sites */
      if (use_pinv && prop_invar[i] > 0)
      {
        __m256d v_inv_prop = _mm256_set1_pd(1. - prop_invar[i]);
        v_cat_sitelk = _mm256_mul_pd(v_cat_sitelk, v_inv_prop);

        if (invariant && *invariant_ptr != -1)
        {
          double site_invar_lk = invar_lk[(*invariant_ptr) * rate_cats + i];
          __m256d v_inv_lk = _mm256_setr_pd(site_invar_lk, 0., 0., 0.);
          v_cat_sitelk = _mm256_add_pd(v_cat_sitelk, v_inv_lk);
        }
      }

      /* apply rate category weights */
      if (eq_weights)
      {
        /* all rate weights are equal -> no multiplication needed */
        v_sitelk = _mm256_add_pd (v_sitelk, v_cat_sitelk);
      }
      else
      {
        __m256d v_weight = _mm256_set1_pd(rate_weights[i]);
        v_sitelk = _mm256_fmadd_pd(v_cat_sitelk, v_weight, v_sitelk);
      }
    }

    _mm256_store_pd(&site_lk[offset], v_sitelk);
    offset += 4;

    invariant_ptr++;

    /* build derivatives for 4 adjacent sites at once */
    if (offset == 16)
    {
      __m256d v_term0 = _mm256_setr_pd(site_lk[0], site_lk[4],
                                       site_lk[8], site_lk[12]);
      __m256d v_term1 = _mm256_setr_pd(site_lk[1], site_lk[5],
                                       site_lk[9], site_lk[13]);
      __m256d v_term2 = _mm256_setr_pd(site_lk[2], site_lk[6],
                                       site_lk[10], site_lk[14]);

      __m256d v_recip0 = _mm256_div_pd(v_all1, v_term0);
      __m256d v_deriv1 = _mm256_mul_pd(v_term1, v_recip0);
      __m256d v_deriv2 = _mm256_sub_pd(_mm256_mul_pd(v_deriv1, v_deriv1),
                                       _mm256_mul_pd(v_term2, v_recip0));

      /* assumption: no zero weights */
      //if ((pattern_weights[n-3] | pattern_weights[n-2] |
      //     pattern_weights[n-1] | pattern_weights[n]) == 1)
      if (pattern_weights[n-3] == 1.0 && pattern_weights[n-2] == 1.0 &&
          pattern_weights[n-1] == 1.0 && pattern_weights[n] == 1.0)
      {
        /* all 4 weights are 1 -> no multiplication needed */
        v_df = _mm256_sub_pd (v_df, v_deriv1);
        v_ddf = _mm256_add_pd (v_ddf, v_deriv2);
      }
      else
      {
        __m256d v_patw = _mm256_setr_pd(pattern_weights[n-3], pattern_weights[n-2],
                                        pattern_weights[n-1], pattern_weights[n]);

        v_df = _mm256_fnmadd_pd (v_deriv1, v_patw, v_df);
        v_ddf = _mm256_fmadd_pd (v_deriv2, v_patw, v_ddf);
      }
      offset = 0;
    }
  }

  *d_f = *dd_f = 0.;

  /* remainder loop */
  while (offset > 0)
  {
    offset -= 4;
    n--;
    double deriv1 = (-site_lk[offset+1] / site_lk[offset]);
    double deriv2 = (deriv1 * deriv1 - (site_lk[offset+2] / site_lk[offset]));
    *d_f += pattern_weights[n] * deriv1;
    *dd_f += pattern_weights[n] * deriv2;
  }

  assert(offset == 0 && n == ef_sites / 4 * 4);

  /* reduce 1st derivative */
  _mm256_store_pd(site_lk, v_df);
  *d_f += site_lk[0] + site_lk[1] + site_lk[2] + site_lk[3];

  /* reduce 2nd derivative */
  _mm256_store_pd(site_lk, v_ddf);
  *dd_f += site_lk[0] + site_lk[1] + site_lk[2] + site_lk[3];

  if (t_diagp)
    pll_aligned_free(t_diagp);
  if (invar_lk)
    pll_aligned_free(invar_lk);

  return PLL_SUCCESS;
}
