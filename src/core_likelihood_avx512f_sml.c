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

__m512d __mm512_log_pd(__m512d x);

//TODO: Only temporary
void print_256d(const __m256d v) {
  const double *val = (const double *) &v;
  printf("% .3e % .3e % .3e % .3e\n",
         val[3], val[2], val[1], val[0]);
}

void print_256i(const __m256i v) {
  const int32_t *val = (const int32_t *) &v;
  printf("%d %d %d %d %d %d %d %d\n",
         val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void print_256i_u(const __m256i v) {
  const uint32_t *val = (const uint32_t *) &v;
  printf("%u %u %u %u %u %u %u %u\n",
         val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void print_512d(const __m512d v) {
  const double *val = (const double *) &v;
  printf("% .3e % .3e % .3e % .3e % .3e % .3e % .3e % .3e\n",
         val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void print_512i(const __m512i v) {
  const int64_t *val = (const int64_t *) &v;
  printf("%ld %ld %ld %ld %ld %ld %ld %ld\n",
         val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void print_512i_u(const __m512i v) {
  const uint64_t *val = (const uint64_t *) &v;
  printf("%lu %lu %lu %lu %lu %lu %lu %lu\n",
         val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
}

void print_512d_half(const __m512d v, size_t half) {
  const double *val = (const double *) &v;
  printf("% .3e % .3e % .3e % .3e\n",
         val[3 + 4 * half], val[2 + 4 * half], val[1 + 4 * half], val[0 + 4 * half]);
}

inline double reduce_add_pd(const __m512d zmm) {
  __m256d low = _mm512_castpd512_pd256(zmm);
  __m256d high = _mm512_extractf64x4_pd(zmm, 1);

  __m256d a = _mm256_add_pd(low, high);
  __m256d t1 = _mm256_hadd_pd(a, a);
  __m128d t2 = _mm256_extractf128_pd(t1, 1);
  __m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1), t2);
  return _mm_cvtsd_f64(t3);
}

#define SAVE_LOG(val) (val)<=0 ? 0 : log((val))

inline __m512d __mm512_log_pd(__m512d x) {
  const double *val = (const double *) &x;
  return _mm512_set_pd(SAVE_LOG(val[7]), SAVE_LOG(val[6]), SAVE_LOG(val[5]), SAVE_LOG(val[4]),
                       SAVE_LOG(val[3]), SAVE_LOG(val[2]), SAVE_LOG(val[1]), SAVE_LOG(val[0]));
}

inline __m512i __mm512_load_epi32_to_epi64(const int * x) {
  return _mm512_set_epi64(x[7], x[6], x[5], x[4],
                          x[3], x[2], x[1], x[0]);
}


inline __m512d __mm512_load_pattern_weights(const unsigned int *x) {
  return _mm512_set_pd(x[7], x[6], x[5], x[4],
                       x[3], x[2], x[1], x[0]);
}

PLL_EXPORT double pll_core_root_loglikelihood_avx512f_sml(unsigned int states,
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
  assert(0); // TODO not supported yet :(
  unsigned int i, j, k;
  double logl = 0;
  double prop_invar = 0;

  const double *freqs = NULL;

  double term, term_r;
  double inv_site_lk;

  unsigned int states_padded = (states + 7) & (0xFFFFFFFF - 7);

  __m512d xmm0, xmm1, xmm3;

  for (i = 0; i < sites; ++i) {
    term = 0;
    for (j = 0; j < rate_cats; ++j) {
      freqs = frequencies[freqs_indices[j]];
      xmm3 = _mm512_setzero_pd();

      for (k = 0; k < states_padded; k += ELEM_PER_AVX515_REGISTER) {
        /* load frequencies for current rate matrix */
        xmm0 = _mm512_load_pd(freqs);

        /* load clv */
        xmm1 = _mm512_load_pd(clv);

        /* multiply with frequencies */
        xmm3 = _mm512_fmadd_pd(xmm0, xmm1, xmm3);

        freqs += ELEM_PER_AVX515_REGISTER;
        clv += ELEM_PER_AVX515_REGISTER;
      }

      /* add up the elements of xmm2 */
      term_r = reduce_add_pd(xmm3);

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

PLL_EXPORT double pll_core_root_loglikelihood_repeats_avx512f_sml(unsigned int states,
                                                                  unsigned int sites,
                                                                  unsigned int rate_cats,
                                                                  const double *clv,
                                                                  const unsigned int *site_id,
                                                                  const unsigned int *scaler,
                                                                  double *const *frequencies,
                                                                  const double *rate_weights,
                                                                  const unsigned int *pattern_weights,
                                                                  const double *invar_proportion,
                                                                  const int *invar_indices,
                                                                  const unsigned int *freqs_indices,
                                                                  double *persite_lnl) {
  //TODO: not implemented
  assert(0);
  return 0;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_20x20_avx512f_sml(unsigned int sites,
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
double pll_core_edge_loglikelihood_ii_avx512f_sml(unsigned int states,
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
  __m512d v_logl = _mm512_setzero_pd();

  const double * clvp = parent_clv;
  const double * clvc = child_clv;

  unsigned int states_padded = (states + 7) & (0xFFFFFFFF - 7);

  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    double scale_factor = 1.0;
    for (unsigned int i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }

  unsigned int site_scalings;
  for (unsigned int n = 0; n < sites; n+= ELEM_PER_AVX515_REGISTER)
  {
    if (per_rate_scaling)
    {
      assert(0); // TODO not implemented yet :(
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (unsigned int i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (unsigned int i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
      site_scalings += (child_scaler) ? child_scaler[n] : 0;
    }

    const double *pmat = pmatrix;

    __m512d terma = _mm512_setzero_pd();
    for (unsigned int i = 0; i < rate_cats; ++i)
    {
      const double *freqs = frequencies[freqs_indices[i]];

      __m512d terma_r = _mm512_setzero_pd();
      for (unsigned int j = 0; j < states; ++j)
      {
        __m512d termb = _mm512_setzero_pd();
        for (unsigned int k = 0; k < states; ++k)
        {
          __m512d v_pmat = _mm512_set1_pd(pmat[k]);

          __m512d v_clvc = _mm512_load_pd(clvc + ELEM_PER_AVX515_REGISTER * k);
          termb = _mm512_fmadd_pd(v_pmat, v_clvc, termb);
        }

        __m512d v_clvp = _mm512_load_pd(clvp + ELEM_PER_AVX515_REGISTER * j);
        __m512d v_freqs = _mm512_set1_pd(freqs[j]);
        terma_r = _mm512_fmadd_pd(_mm512_mul_pd(v_clvp, v_freqs), termb, terma_r);

        pmat += states_padded;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        assert(0); //TODO not implemented yet :(
        //terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      double prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {

        __m512i v_invar_indices =  __mm512_load_epi32_to_epi64(invar_indices + n);
        __mmask8 valid_invar_indices = _mm512_cmpneq_epi64_mask(v_invar_indices, _mm512_set1_epi64(-1));

        __m512d v_inv_site_lk = _mm512_mask_i64gather_pd(_mm512_setzero_pd(),
                                                         valid_invar_indices,
                                                         v_invar_indices,
                                                         freqs,
                                                         sizeof(double));
        terma = _mm512_fmadd_pd(_mm512_set1_pd(rate_weights[i]),
                                _mm512_add_pd(_mm512_mul_pd(terma_r, _mm512_set1_pd(1-prop_invar)),
                                              _mm512_mul_pd(v_inv_site_lk, _mm512_set1_pd(prop_invar))),
                                terma);
      }
      else
      {
        terma = _mm512_fmadd_pd(terma_r, _mm512_set1_pd(rate_weights[i]), terma);
      }

      clvp += states * ELEM_PER_AVX515_REGISTER;
      clvc += states * ELEM_PER_AVX515_REGISTER;
    }

    /* compute site log-likelihood and scale if necessary */
    __m512d site_lk = __mm512_log_pd(terma);
    if (site_scalings) {
      assert(0); // TODO: not supported yet :(
      //site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);
    }

    __m512d v_pattern_weights = __mm512_load_pattern_weights(pattern_weights + n);

    site_lk = _mm512_mul_pd(site_lk, v_pattern_weights);

    /* store per-site log-likelihood */
    if (persite_lnl) {
      assert(0); // TODO: Have to assume site padding here
      _mm512_store_pd(persite_lnl + n*ELEM_PER_AVX515_REGISTER, site_lk);
    }

    v_logl = _mm512_add_pd(v_logl, site_lk);
  }

  if (rate_scalings)
    free(rate_scalings);

  return reduce_add_pd(v_logl);
}


PLL_EXPORT
double pll_core_edge_loglikelihood_repeats_generic_avx512f_sml(unsigned int states,
                                                               unsigned int sites,
                                                               const unsigned int child_sites,
                                                               unsigned int rate_cats,
                                                               const double *parent_clv,
                                                               const unsigned int *parent_scaler,
                                                               const double *child_clv,
                                                               const unsigned int *child_scaler,
                                                               const double *pmatrix,
                                                               double **frequencies,
                                                               const double *rate_weights,
                                                               const unsigned int *pattern_weights,
                                                               const double *invar_proportion,
                                                               const int *invar_indices,
                                                               const unsigned int *freqs_indices,
                                                               double *persite_lnl,
                                                               const unsigned int *parent_site_id,
                                                               const unsigned int *child_site_id,
                                                               double *bclv,
                                                               unsigned int attrib) {
  //TODO: not implemented
  assert(!(PLL_ATTRIB_SITE_REPEATS & attrib));
  return 0;
}
