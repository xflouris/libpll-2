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

#include <limits.h>
#include "pll.h"

/*
    TODO: Vectorize these functions which were created to make derivatives
          work with SSE (states_padded is set to the corresponding value)
*/

static int core_update_sumtable_ti_4x4_sse(unsigned int sites,
                                           unsigned int rate_cats,
                                           const double * parent_clv,
                                           const unsigned char * left_tipchars,
                                           const unsigned int * parent_scaler,
                                           double * const * eigenvecs,
                                           double * const * inv_eigenvecs,
                                           double * const * freqs_indices,
                                           double * sumtable,
                                           unsigned int attrib)
{
  unsigned int i,j,k,n;
  unsigned int tipstate;
  unsigned int states = 4;
  double lterm = 0;
  double rterm = 0;

  const double * clvc = parent_clv;
  const double * ev;
  const double * invei;
  const double * freqs;

  double * sum = sumtable;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
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

    for (i = 0; i < rate_cats; ++i)
    {
      ev    = eigenvecs[i];
      invei = inv_eigenvecs[i];
      freqs = freqs_indices[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = (unsigned int)left_tipchars[n];
        lterm = 0;
        rterm = 0;

        for (k = 0; k < states; ++k)
        {
          lterm += (tipstate & 1) * freqs[k] * invei[k*states+j];
          rterm += ev[j*states+k] * clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lterm*rterm;

        if (rate_scalings && rate_scalings[i] > 0)
          sum[j] *= scale_minlh[rate_scalings[i]-1];
      }

      clvc += states;
      sum += states;
    }
  }

  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_repeats_generic_sse(unsigned int states,
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

  double * sum = sumtable;

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;
  unsigned int span_padded = rate_cats * states_padded;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m128d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  double * tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_SSE);

  /* transposed padded inv_eigenvecs */
  double * tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_SSE);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    const double *t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    unsigned int pid = PLL_GET_ID(parent_site_id, n);
    unsigned int cid = PLL_GET_ID(child_site_id, n);
    const double * t_clvp = &clvp[pid * span_padded];
    const double * t_clvc = &clvc[cid * span_padded];
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
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
      for (j = 0; j < states; j += 2)
      {
        /* point to the two rows of the eigenvecs matrix */
        const double * em0 = c_eigenvecs;
        const double * em1 = em0 + states_padded;
        c_eigenvecs += 2*states_padded;

        /* point to the two rows of the inv_eigenvecs matrix */
        const double * im0 = ct_inv_eigenvecs;
        const double * im1 = im0 + states_padded;
        ct_inv_eigenvecs += 2*states_padded;

        __m128d v_lefterm0 = _mm_setzero_pd ();
        __m128d v_righterm0 = _mm_setzero_pd ();
        __m128d v_lefterm1 = _mm_setzero_pd ();
        __m128d v_righterm1 = _mm_setzero_pd ();

        __m128d v_eigen;
        __m128d v_clvp;
        __m128d v_clvc;

        for (k = 0; k < states_padded; k += 2)
        {
          v_clvp = _mm_load_pd (t_clvp + k);
          v_clvc = _mm_load_pd (t_clvc + k);

          /* row 0 */
          v_eigen = _mm_load_pd (im0 + k);
          v_lefterm0 = _mm_add_pd (v_lefterm0,
                                   _mm_mul_pd (v_eigen, v_clvp));

          v_eigen = _mm_load_pd (em0 + k);
          v_righterm0 = _mm_add_pd (v_righterm0,
                                    _mm_mul_pd (v_eigen, v_clvc));

          /* row 1 */
          v_eigen = _mm_load_pd (im1 + k);
          v_lefterm1 = _mm_add_pd (v_lefterm1,
                                   _mm_mul_pd (v_eigen, v_clvp));

          v_eigen = _mm_load_pd (em1 + k);
          v_righterm1 = _mm_add_pd (v_righterm1,
                                    _mm_mul_pd (v_eigen, v_clvc));
        }

        /* update sum */
        __m128d v_lefterm_sum = _mm_hadd_pd(v_lefterm0,v_lefterm1);
        __m128d v_righterm_sum = _mm_hadd_pd(v_righterm0,v_righterm1);
        __m128d v_prod = _mm_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm_store_pd (sum + j, v_prod);
      }
      t_clvc += states_padded;
      t_clvp += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ii_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const double * child_clv,
                                               const unsigned int * parent_scaler,
                                               const unsigned int * child_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs,
                                               double * sumtable,
                                               unsigned int attrib)
{
  unsigned int i, j, k, n;

  const double * clvp = parent_clv;
  const double * clvc = child_clv;
  double * sum = sumtable;

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  __m128d v_scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      v_scale_minlh[i] = _mm_set1_pd(scale_factor);
    }
  }

  /* padded eigenvecs */
  double * tt_eigenvecs = (double *) pll_aligned_alloc (
        (states_padded * states_padded * rate_cats) * sizeof(double),
        PLL_ALIGNMENT_SSE);

  /* transposed padded inv_eigenvecs */
  double * tt_inv_eigenvecs = (double *) pll_aligned_alloc (
      (states_padded * states_padded * rate_cats) * sizeof(double),
      PLL_ALIGNMENT_SSE);

  if (!tt_eigenvecs || !tt_inv_eigenvecs)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200, "Cannot allocate memory for tt_inv_eigenvecs");
    return PLL_FAILURE;
  }

  memset(tt_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));
  memset(tt_inv_eigenvecs, 0, (states_padded * states_padded * rate_cats) * sizeof(double));

  /* add padding to eigenvecs matrices and multiply with frequencies */
  for (i = 0; i < rate_cats; ++i)
  {
    const double *t_freqs = freqs[i];
    for (j = 0; j < states; ++j)
      for (k = 0; k < states; ++k)
      {
        tt_inv_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = inv_eigenvecs[i][k * states_padded + j] * t_freqs[k];
        tt_eigenvecs[i * states_padded * states_padded + j * states_padded
            + k] = eigenvecs[i][j * states_padded + k];
      }
  }

  /* build sumtable */
  for (n = 0; n < sites; n++)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
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

      for (j = 0; j < states; j += 2)
      {
        /* point to the two rows of the eigenvecs matrix */
        const double * em0 = c_eigenvecs;
        const double * em1 = em0 + states_padded;
        c_eigenvecs += 2*states_padded;

        /* point to the two rows of the inv_eigenvecs matrix */
        const double * im0 = ct_inv_eigenvecs;
        const double * im1 = im0 + states_padded;
        ct_inv_eigenvecs += 2*states_padded;

        __m128d v_lefterm0 = _mm_setzero_pd ();
        __m128d v_righterm0 = _mm_setzero_pd ();
        __m128d v_lefterm1 = _mm_setzero_pd ();
        __m128d v_righterm1 = _mm_setzero_pd ();

        __m128d v_eigen;
        __m128d v_clvp;
        __m128d v_clvc;

        for (k = 0; k < states_padded; k += 2)
        {
          v_clvp = _mm_load_pd (clvp + k);
          v_clvc = _mm_load_pd (clvc + k);

          /* row 0 */
          v_eigen = _mm_load_pd (im0 + k);
          v_lefterm0 = _mm_add_pd (v_lefterm0,
                                   _mm_mul_pd (v_eigen, v_clvp));

          v_eigen = _mm_load_pd (em0 + k);
          v_righterm0 = _mm_add_pd (v_righterm0,
                                    _mm_mul_pd (v_eigen, v_clvc));

          /* row 1 */
          v_eigen = _mm_load_pd (im1 + k);
          v_lefterm1 = _mm_add_pd (v_lefterm1,
                                   _mm_mul_pd (v_eigen, v_clvp));

          v_eigen = _mm_load_pd (em1 + k);
          v_righterm1 = _mm_add_pd (v_righterm1,
                                    _mm_mul_pd (v_eigen, v_clvc));
        }

        /* update sum */
        __m128d v_lefterm_sum = _mm_hadd_pd(v_lefterm0,v_lefterm1);
        __m128d v_righterm_sum = _mm_hadd_pd(v_righterm0,v_righterm1);
        __m128d v_prod = _mm_mul_pd (v_lefterm_sum, v_righterm_sum);

        /* apply per-rate scalers */
        if (rate_scalings && rate_scalings[i] > 0)
        {
          v_prod = _mm_mul_pd(v_prod, v_scale_minlh[rate_scalings[i]-1]);
        }

        _mm_store_pd (sum + j, v_prod);
      }

      clvc += states_padded;
      clvp += states_padded;
      sum += states_padded;
    }
  }

  pll_aligned_free (tt_inv_eigenvecs);
  pll_aligned_free (tt_eigenvecs);
  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_sumtable_ti_sse(unsigned int states,
                                               unsigned int sites,
                                               unsigned int rate_cats,
                                               const double * parent_clv,
                                               const unsigned char * left_tipchars,
                                               const unsigned int * parent_scaler,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               double * const * freqs_indices,
                                               const pll_state_t * tipmap,
                                               double * sumtable,
                                               unsigned int attrib)
{
  unsigned int i,j,k,n;
  pll_state_t tipstate;
  double lterm = 0;
  double rterm = 0;

  double * sum = sumtable;
  const double * clvc = parent_clv;
  const double * ev;
  const double * invev;
  const double * freqs;

  if (states == 4)
  {
    return core_update_sumtable_ti_4x4_sse(sites,
                                           rate_cats,
                                           parent_clv,
                                           left_tipchars,
                                           parent_scaler,
                                           eigenvecs,
                                           inv_eigenvecs,
                                           freqs_indices,
                                           sumtable,
                                           attrib);
  }

  unsigned int states_padded = (states+1) & 0xFFFFFFFE;

  unsigned int min_scaler;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  /* build sumtable: non-vectorized version, general case */
  for (n = 0; n < sites; n++)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
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

    for (i = 0; i < rate_cats; ++i)
    {
      ev    = eigenvecs[i];
      invev = inv_eigenvecs[i];
      freqs = freqs_indices[i];

      for (j = 0; j < states; ++j)
      {
        tipstate = tipmap[(unsigned int)left_tipchars[n]];
        lterm = 0;
        rterm = 0;

        for (k = 0; k < states; ++k)
        {
          lterm += (tipstate & 1) * freqs[k] * invev[k*states_padded+j];
          rterm += ev[j*states_padded+k] * clvc[k];
          tipstate >>= 1;
        }
        sum[j] = lterm*rterm;

        if (rate_scalings && rate_scalings[i] > 0)
          sum[j] *= scale_minlh[rate_scalings[i]-1];
      }

      clvc += states_padded;
      sum += states_padded;
    }
  }

  if (rate_scalings)
    free(rate_scalings);

  return PLL_SUCCESS;
}
