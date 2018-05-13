/*
    Copyright (C) 2015 Tomas Flouri

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

PLL_EXPORT double pll_core_root_loglikelihood_avx512f(unsigned int states,
                                                      unsigned int sites,
                                                      unsigned int rate_cats,
                                                      const double *clv,
                                                      const unsigned int *scaler,
                                                      double *const *frequencies,
                                                      const double *rate_weights,
                                                      const unsigned int *pattern_weights,
                                                      const double *invar_proportion,
                                                      const int *invar_indices,
                                                      const unsigned int *freqs_indices,
                                                      double *persite_lnl) {
  unsigned int i, j, k;
  double logl = 0;
  double prop_invar = 0;

  const double *freqs = NULL;

  double term, term_r;
  double inv_site_lk;

  unsigned int states_padded = (states + 7) & (0xFFFFFFFF - 7);

  __m512d xmm0, xmm1, xmm3;

  double reg[ELEM_PER_AVX512_REGISTER] __attribute__( ( aligned ( PLL_ALIGNMENT_AVX512F ) ) ) ;

  for (i = 0; i < sites; ++i) {
    term = 0;
    for (j = 0; j < rate_cats; ++j) {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm512_setzero_pd();

      for (k = 0; k < states_padded; k += ELEM_PER_AVX512_REGISTER) {
        /* load frequencies for current rate matrix */
        xmm0 = _mm512_load_pd(freqs);

        /* load clv */
        xmm1 = _mm512_load_pd(clv);

        /* multiply with frequencies */
        xmm3 = _mm512_fmadd_pd(xmm0, xmm1, xmm3);

        freqs += ELEM_PER_AVX512_REGISTER;
        clv += ELEM_PER_AVX512_REGISTER;
      }

      /* add up the elements of xmm2 */

      _mm512_store_pd(reg, xmm3);
      term_r = reg[0] + reg[1] + reg[2] + reg[3] + reg[4] + reg[5] + reg[6] + reg[7];

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[j]] : 0;
      if (prop_invar > 0) {
        freqs = frequencies[freqs_indices[j]];
        inv_site_lk = (invar_indices[i] == -1) ?
                      0 : freqs[invar_indices[i]];
        term += rate_weights[j] * (term_r * (1 - prop_invar) +
                                   inv_site_lk * prop_invar);
      } else {
        term += term_r * rate_weights[j];
      }
    }

    /* compute site log-likelihood and scale if necessary */
    term = log(term);
    if (scaler && scaler[i])
      term += scaler[i] * log(PLL_SCALE_THRESHOLD);

    term *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = term;

    logl += term;
  }
  return logl;
}

PLL_EXPORT double pll_core_root_loglikelihood_repeats_avx512f(unsigned int states,
                                                              unsigned int sites,
                                                              unsigned int rate_cats,
                                                              const double * clv,
                                                              const unsigned int * site_id,
                                                              const unsigned int * scaler,
                                                              double * const * frequencies,
                                                              const double * rate_weights,
                                                              const unsigned int * pattern_weights,
                                                              const double * invar_proportion,
                                                              const int * invar_indices,
                                                              const unsigned int * freqs_indices,
                                                              double * persite_lnl)
{
	//TODO: not implemented
	assert(0);
	return 0;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_20x20_avx512f(unsigned int sites,
                                                    unsigned int rate_cats,
                                                    const double *parent_clv,
                                                    const unsigned int *parent_scaler,
                                                    const unsigned char *tipchars,
                                                    const pll_state_t *tipmap,
                                                    unsigned int tipmap_size,
                                                    const double *pmatrix,
                                                    double *const *frequencies,
                                                    const double *rate_weights,
                                                    const unsigned int *pattern_weights,
                                                    const double *invar_proportion,
                                                    const int *invar_indices,
                                                    const unsigned int *freqs_indices,
                                                    double *persite_lnl,
                                                    unsigned int attrib) {
  //TODO not implemented
  assert(!(attrib & PLL_ATTRIB_PATTERN_TIP));
  return 0;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ii_avx512f(unsigned int states,
                                              unsigned int sites,
                                              unsigned int rate_cats,
                                              const double *parent_clv,
                                              const unsigned int *parent_scaler,
                                              const double *child_clv,
                                              const unsigned int *child_scaler,
                                              const double *pmatrix,
                                              double *const *frequencies,
                                              const double *rate_weights,
                                              const unsigned int *pattern_weights,
                                              const double *invar_proportion,
                                              const int *invar_indices,
                                              const unsigned int *freqs_indices,
                                              double *persite_lnl,
                                              unsigned int attrib) {
  unsigned int n, i, j, k;
  double logl = 0;
  double prop_invar = 0;

  const double *clvp = parent_clv;
  const double *clvc = child_clv;
  const double *pmat;
  const double *freqs = NULL;

  double terma, terma_r;
  double site_lk, inv_site_lk;

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

  __m512d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7;

  size_t displacement = (states_padded - states) * (states_padded);

  double reg[ELEM_PER_AVX512_REGISTER] __attribute__( ( aligned ( PLL_ALIGNMENT_AVX512F ) ) ) ;

  /* scaling stuff */
  unsigned int site_scalings;
  unsigned int *rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling) {
    rate_scalings = (unsigned int *) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings) {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }

    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i) {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  for (n = 0; n < sites; ++n) {
    pmat = pmatrix;
    terma = 0;

    if (per_rate_scaling) {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i) {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n * rate_cats + i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n * rate_cats + i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i) {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    } else {
      /* count number of scaling factors to account for */
      site_scalings = (parent_scaler) ? parent_scaler[n] : 0;
      site_scalings += (child_scaler) ? child_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i) {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;

      /* iterate over octuples of rows */
      for (j = 0; j < states_padded; j += ELEM_PER_AVX512_REGISTER) {
        xmm0 = _mm512_setzero_pd();
        xmm1 = _mm512_setzero_pd();
        xmm2 = _mm512_setzero_pd();
        xmm3 = _mm512_setzero_pd();
        xmm4 = _mm512_setzero_pd();
        xmm5 = _mm512_setzero_pd();
        xmm6 = _mm512_setzero_pd();
        xmm7 = _mm512_setzero_pd();

        /* point to the four rows */
        const double *row0 = pmat;
        const double *row1 = row0 + states_padded;
        const double *row2 = row1 + states_padded;
        const double *row3 = row2 + states_padded;
        const double *row4 = row3 + states_padded;
        const double *row5 = row4 + states_padded;
        const double *row6 = row5 + states_padded;
        const double *row7 = row6 + states_padded;

        /* iterate quadruples of columns */
        for (k = 0; k < states_padded; k += ELEM_PER_AVX512_REGISTER) {
          __m512d v_clvc = _mm512_load_pd(clvc + k);

          /* row 0 */
          __m512d v_row0 = _mm512_load_pd(row0);
          xmm0 = _mm512_fmadd_pd(v_row0, v_clvc, xmm0);
          row0 += ELEM_PER_AVX512_REGISTER;

          /* row 1 */
          __m512d v_row1 = _mm512_load_pd(row1);
          xmm1 = _mm512_fmadd_pd(v_row1, v_clvc, xmm1);
          row1 += ELEM_PER_AVX512_REGISTER;

          /* row 2 */
          __m512d v_row2 = _mm512_load_pd(row2);
          xmm2 = _mm512_fmadd_pd(v_row2, v_clvc, xmm2);
          row2 += ELEM_PER_AVX512_REGISTER;

          /* row 3 */
          __m512d v_row3 = _mm512_load_pd(row3);
          xmm3 = _mm512_fmadd_pd(v_row3, v_clvc, xmm3);
          row3 += ELEM_PER_AVX512_REGISTER;

          /* row 4 */
          __m512d v_row4 = _mm512_load_pd(row4);
          xmm4 = _mm512_fmadd_pd(v_row4, v_clvc, xmm4);
          row4 += ELEM_PER_AVX512_REGISTER;

          /* row 5 */
          __m512d v_row5 = _mm512_load_pd(row5);
          xmm5 = _mm512_fmadd_pd(v_row5, v_clvc, xmm5);
          row5 += ELEM_PER_AVX512_REGISTER;

          /* row 6 */
          __m512d v_row6 = _mm512_load_pd(row6);
          xmm6 = _mm512_fmadd_pd(v_row6, v_clvc, xmm6);
          row6 += ELEM_PER_AVX512_REGISTER;

          /* row 7 */
          __m512d v_row7 = _mm512_load_pd(row7);
          xmm7 = _mm512_fmadd_pd(v_row7, v_clvc, xmm7);
          row7 += ELEM_PER_AVX512_REGISTER;
        }

        /* point pmatrix to the next four rows */
        pmat = row7;

        /* create a vector containing the sums of xmm0 - xmm7 */
        __m512d ymm0 = _mm512_add_pd(_mm512_unpackhi_pd(xmm0, xmm1),
                                     _mm512_unpacklo_pd(xmm0, xmm1));
        __m512d ymm1 = _mm512_add_pd(_mm512_unpackhi_pd(xmm2, xmm3),
                                     _mm512_unpacklo_pd(xmm2, xmm3));
        __m512d ymm2 = _mm512_add_pd(_mm512_unpackhi_pd(xmm4, xmm5),
                                     _mm512_unpacklo_pd(xmm4, xmm5));
        __m512d ymm3 = _mm512_add_pd(_mm512_unpackhi_pd(xmm6, xmm7),
                                     _mm512_unpacklo_pd(xmm6, xmm7));

        __m512d zmm0 = _mm512_add_pd(_mm512_permutex2var_pd(ymm0, permute_mask, ymm2),
                                     _mm512_mask_blend_pd(0xF0, ymm0, ymm2));

        __m512d zmm1 = _mm512_add_pd(_mm512_permutex2var_pd(ymm1, permute_mask, ymm3),
                                     _mm512_mask_blend_pd(0xF0, ymm1, ymm3));

        __m512d v_sum = _mm512_add_pd(_mm512_permutex2var_pd(zmm0,
                                                             permute_mask_final_stage,
                                                             zmm1),
                                      _mm512_mask_blend_pd(0xCC, zmm0, zmm1));

        /* multiply with frequencies */
        xmm1 = _mm512_load_pd(freqs);
        xmm2 = _mm512_mul_pd(v_sum, xmm1);

        /* multiply with clvp */
        xmm0 = _mm512_load_pd(clvp);
        xmm1 = _mm512_mul_pd(xmm2, xmm0);

        /* add up the elements of xmm1 */
        _mm512_store_pd(reg, xmm1);
        terma_r += reg[0] + reg[1] + reg[2] + reg[3] + reg[4] + reg[5] + reg[6] + reg[7];

        freqs += ELEM_PER_AVX512_REGISTER;
        clvp += ELEM_PER_AVX512_REGISTER;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0) {
        terma_r *= scale_minlh[rate_scalings[i] - 1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0) {
        freqs = frequencies[freqs_indices[i]];
        inv_site_lk = (invar_indices[n] == -1) ?
                      0 : freqs[invar_indices[n]];
        terma += rate_weights[i] * (terma_r * (1 - prop_invar) +
                                    inv_site_lk * prop_invar);
      } else {
        terma += terma_r * rate_weights[i];
      }

      clvc += states_padded;
      pmat -= displacement;
    }

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(terma);
    if (site_scalings)
      site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_repeats_generic_avx512f(unsigned int states,
                                                           unsigned int sites,
                                                           const unsigned int child_sites,
                                                           unsigned int rate_cats,
                                                           const double * parent_clv,
                                                           const unsigned int * parent_scaler,
                                                           const double * child_clv,
                                                           const unsigned int * child_scaler,
                                                           const double * pmatrix,
                                                           double ** frequencies,
                                                           const double * rate_weights,
                                                           const unsigned int * pattern_weights,
                                                           const double * invar_proportion,
                                                           const int * invar_indices,
                                                           const unsigned int * freqs_indices,
                                                           double * persite_lnl,
                                                           const unsigned int * parent_site_id,
                                                           const unsigned int * child_site_id,
                                                           double * bclv,
                                                           unsigned int attrib)
{
	//TODO: not implemented
	assert(!(PLL_ATTRIB_SITE_REPEATS & attrib));
	return 0;
}
