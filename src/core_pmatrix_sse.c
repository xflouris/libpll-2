/*
    Copyright (C) 2015 Tomas Flouri, Alexey Kozlov

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

#include "pll.h"

#define ONESTEP4(x)                                                            \
    /* compute pmat row x/4 */                                                 \
    xmm12 = _mm_load_pd(inv_evecs+x);                                          \
    xmm13 = _mm_load_pd(inv_evecs+x+2);                                        \
    xmm12 = _mm_mul_pd(xmm12,xmm1);          /* temp row x/4 (0-1) */          \
    xmm13 = _mm_mul_pd(xmm13,xmm2);          /* temp row x/4 (2-3) */          \
                                                                               \
    /* multiply with row 0 of transposed eigenvector */                        \
    xmm14 = _mm_mul_pd(xmm12,xmm4);                                            \
    xmm15 = _mm_mul_pd(xmm13,xmm8);                                            \
    xmm14 = _mm_add_pd(xmm14,xmm15);                                           \
                                                                               \
    /* multiply with row 1 of transposed eigenvector */                        \
    xmm15 = _mm_mul_pd(xmm12,xmm5);                                            \
    xmm16 = _mm_mul_pd(xmm13,xmm9);                                            \
    xmm15 = _mm_add_pd(xmm15,xmm16);                                           \
                                                                               \
    xmm16 = _mm_hadd_pd(xmm14,xmm15);                                          \
    _mm_store_pd(pmat+x,xmm16);                                                \
                                                                               \
    /* multiply with row 2 of transposed eigenvector */                        \
    xmm14 = _mm_mul_pd(xmm12,xmm6);                                            \
    xmm15 = _mm_mul_pd(xmm13,xmm10);                                           \
    xmm14 = _mm_add_pd(xmm14,xmm15);                                           \
                                                                               \
    /* multiply with row 3 of transposed eigenvector */                        \
    xmm15 = _mm_mul_pd(xmm12,xmm7);                                            \
    xmm16 = _mm_mul_pd(xmm13,xmm11);                                           \
    xmm15 = _mm_add_pd(xmm15,xmm16);                                           \
                                                                               \
    xmm16 = _mm_hadd_pd(xmm14,xmm15);                                          \
    _mm_store_pd(pmat+x+2,xmm16);                                              \


#define ONESTEP20(x,baseptr)                                 \
            ymm0 = _mm_load_pd(baseptr+0);                   \
            ymm1 = _mm_load_pd(baseptr+2);                   \
            ymm2 = _mm_load_pd(baseptr+4);                   \
            ymm3 = _mm_load_pd(baseptr+6);                   \
            ymm4 = _mm_load_pd(baseptr+8);                   \
            ymm5 = _mm_load_pd(baseptr+10);                  \
            ymm6 = _mm_load_pd(baseptr+12);                  \
            ymm7 = _mm_load_pd(baseptr+14);                  \
            ymm8 = _mm_load_pd(baseptr+16);                  \
            ymm9 = _mm_load_pd(baseptr+18);                  \
                                                             \
            ymm0 = _mm_mul_pd(xmm0,ymm0);                    \
            ymm1 = _mm_mul_pd(xmm1,ymm1);                    \
            ymm2 = _mm_mul_pd(xmm2,ymm2);                    \
            ymm3 = _mm_mul_pd(xmm3,ymm3);                    \
            ymm4 = _mm_mul_pd(xmm4,ymm4);                    \
            ymm5 = _mm_mul_pd(xmm5,ymm5);                    \
            ymm6 = _mm_mul_pd(xmm6,ymm6);                    \
            ymm7 = _mm_mul_pd(xmm7,ymm7);                    \
            ymm8 = _mm_mul_pd(xmm8,ymm8);                    \
            ymm9 = _mm_mul_pd(xmm9,ymm9);                    \
                                                             \
            x = _mm_add_pd(ymm0,ymm1);                       \
            x = _mm_add_pd(x,ymm2);                          \
            x = _mm_add_pd(x,ymm3);                          \
            x = _mm_add_pd(x,ymm4);                          \
            x = _mm_add_pd(x,ymm5);                          \
            x = _mm_add_pd(x,ymm6);                          \
            x = _mm_add_pd(x,ymm7);                          \
            x = _mm_add_pd(x,ymm8);                          \
            x = _mm_add_pd(x,ymm9);                          \

PLL_EXPORT int pll_core_update_pmatrix_4x4_sse(double ** pmatrix,
                                               unsigned int rate_cats,
                                               const double * rates,
                                               const double * branch_lengths,
                                               const unsigned int * matrix_indices,
                                               const unsigned int * params_indices,
                                               const double * prop_invar,
                                               double * const * eigenvals,
                                               double * const * eigenvecs,
                                               double * const * inv_eigenvecs,
                                               unsigned int count)
{
  unsigned int i,j,n;
  double * expd;

  double pinvar;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;

  expd = (double *)pll_aligned_alloc(4*sizeof(double), PLL_ALIGNMENT_SSE);

  if (!expd)
  {
    if (expd) pll_aligned_free(expd);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  __m128d xmm0, xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9;
  __m128d xmm10, xmm11, xmm12, xmm13, xmm14, xmm15, xmm16;

  xmm0 = _mm_setzero_pd();

  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);
    
    xmm3 = _mm_set1_pd(branch_lengths[i]);
    pmat = pmatrix[matrix_indices[i]];

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pinvar = prop_invar[params_indices[n]];
      evecs = eigenvecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        _mm_store_pd(pmat+0, xmm0);
        _mm_store_pd(pmat+2, xmm0);
        _mm_store_pd(pmat+4, xmm0);
        _mm_store_pd(pmat+6, xmm0);
        _mm_store_pd(pmat+8, xmm0);
        _mm_store_pd(pmat+10,xmm0);
        _mm_store_pd(pmat+12,xmm0);
        _mm_store_pd(pmat+14,xmm0);

        pmat[0] = pmat[5] = pmat[10] = pmat[15] = 1;
      }
      else
      {
        /* exponentiate eigenvalues */

        /* 1) load eigenvalues and 2) load rate into all slots of register */
        xmm1 = _mm_load_pd(evals+0);
        xmm2 = _mm_load_pd(evals+2);
        xmm4 = _mm_set1_pd(rates[n]);

        /* multiply eigenvalues with rate */
        xmm5 = _mm_mul_pd(xmm1,xmm4);
        xmm6 = _mm_mul_pd(xmm2,xmm4);

        /* multiply product with  branch length */
        xmm7 = _mm_mul_pd(xmm5,xmm3);
        xmm8 = _mm_mul_pd(xmm6,xmm3);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm1 = _mm_set1_pd(1.0 - pinvar);
          xmm7 = _mm_div_pd(xmm7,xmm1);
          xmm8 = _mm_div_pd(xmm8,xmm1);
        }
          
        /* TODO: implement a vectorized double-precision exponentiation */
        //xmm1 = _mm_exp_pd(xmm7);     /* expd */
        //xmm2 = _mm_exp_pd(xmm8);     /* expd */

        /* for now exponentiate non-vectorized */
        _mm_store_pd(expd+0,xmm7);
        _mm_store_pd(expd+2,xmm8);

        /* transpose eigenvector */
        xmm1 = _mm_load_pd(evecs+0);
        xmm2 = _mm_load_pd(evecs+4);
        xmm4 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 0 (0,1) */
        xmm5 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 1 (0,1) */

        xmm1 = _mm_load_pd(evecs+2);
        xmm2 = _mm_load_pd(evecs+6);
        xmm6 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 2 (0,1) */
        xmm7 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 3 (0,1) */

        xmm1 = _mm_load_pd(evecs+8);
        xmm2 = _mm_load_pd(evecs+12);
        xmm8 = _mm_unpacklo_pd(xmm1,xmm2);     /* row 0 (2,3) */
        xmm9 = _mm_unpackhi_pd(xmm1,xmm2);     /* row 1 (2,3) */

        xmm1 = _mm_load_pd(evecs+10);
        xmm2 = _mm_load_pd(evecs+14);
        xmm10 = _mm_unpacklo_pd(xmm1,xmm2);    /* row 2 (2,3) */
        xmm11 = _mm_unpackhi_pd(xmm1,xmm2);    /* row 3 (2,3) */

        /* NOTE: in order to deal with numerical issues in cases when Qt -> 0, we
         * use a trick suggested by Ben Redelings and explained here:
         * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
         * In short, we use expm1() to compute (exp(Qt) - I), and then correct
         * for this by adding an identity matrix I in the very end */

        /* load exponentiated eigenvalues */
        xmm1 = _mm_set_pd(expm1(expd[1]), expm1(expd[0]));
        xmm2 = _mm_set_pd(expm1(expd[3]), expm1(expd[2]));

        /* compute pmatrix */
        ONESTEP4(0);
        ONESTEP4(4);
        ONESTEP4(8);
        ONESTEP4(12);

        /* add identity matrix */
        for (j = 0; j < 4; ++j)
        {
          pmat[j] += 1.0;
          pmat += 4;
        }
        pmat -= 16;
      }
      #ifdef DEBUG
      unsigned int j,k;
      for (j = 0; j < 4; ++j)
        for (k = 0; k < 4; ++k)
          assert(pmat[j*4+k] >= 0);
      #endif
      pmat = pmat+16;
    }
  }

  pll_aligned_free(expd);
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_core_update_pmatrix_20x20_sse(double ** pmatrix,
                                                 unsigned int rate_cats,
                                                 const double * rates,
                                                 const double * branch_lengths,
                                                 const unsigned int * matrix_indices,
                                                 const unsigned int * params_indices,
                                                 const double * prop_invar,
                                                 double * const * eigenvals,
                                                 double * const * eigenvecs,
                                                 double * const * inv_eigenvecs,
                                                 unsigned int count)
{
  unsigned int i,n,j,k;
  double pinvar;

  int * transposed;
  double * evecs;
  double * inv_evecs;
  double * evals;
  double * pmat;
  double * expd;
  double * temp;
  double ** tran_evecs;

  expd = (double *)pll_aligned_alloc(20*sizeof(double), PLL_ALIGNMENT_SSE);
  temp = (double *)pll_aligned_alloc(400*sizeof(double), PLL_ALIGNMENT_SSE);

  /* transposed eigen vectors */
  transposed = (int *)calloc((size_t)rate_cats, sizeof(int));
  tran_evecs= (double **)calloc((size_t)rate_cats, sizeof(double *));

  if (!expd || !temp || !transposed || !tran_evecs)
  {
    if (expd) pll_aligned_free(expd);
    if (temp) pll_aligned_free(temp);
    if (transposed) free(transposed);
    if (tran_evecs) free(tran_evecs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* transpose eigenvectors */
  /* TODO: The same trick can be applied for exponentiations */
  for (n = 0; n < rate_cats; ++n)
  {
    int index = params_indices[n];

    if (!transposed[index])
    {
      /* allocate space for transposed eigenvectors and check that
         allocation succeeds */
      double * tran = (double *)pll_aligned_alloc(400*sizeof(double),
                                                  PLL_ALIGNMENT_SSE);
      if (!tran)
      {
        pll_aligned_free(expd);
        pll_aligned_free(temp);
        free(transposed);
        for (i = 0; i < n; ++i)
          if (tran_evecs[i]) pll_aligned_free(tran_evecs[i]);
        free(tran_evecs);

        pll_errno = PLL_ERROR_MEM_ALLOC;
        snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
        return PLL_FAILURE;
      }

      /* transpose eigen vectors */
      evecs = eigenvecs[index];
      for (i = 0; i < 20; ++i)
      {
        for (j = 0; j < 20; ++j)
          tran[i*20+j] = evecs[j*20+i];
      }

      /* update pointers and indicate that the eigen vector for the current
         rate matrix with index was updated */
      tran_evecs[index] = tran;
      transposed[index] = 1;
    }
  }
  free(transposed);

  __m128d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6,xmm7,xmm8,xmm9;
  __m128d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6,ymm7,ymm8,ymm9;
  __m128d zmm0,zmm1,zmm2;
  __m128d brlen, rate;

  double * tran = NULL;
  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    brlen = _mm_set1_pd(branch_lengths[i]);
    pmat = pmatrix[matrix_indices[i]];

    /* compute effective pmatrix location */
    for (n = 0; n < rate_cats; ++n)
    {
      pinvar = prop_invar[params_indices[n]];
      tran = tran_evecs[params_indices[n]];
      inv_evecs = inv_eigenvecs[params_indices[n]];
      evals = eigenvals[params_indices[n]];

      /* if branch length is zero then set the p-matrix to identity matrix */
      if (!branch_lengths[i])
      {
        xmm0 = _mm_setzero_pd();
        for (j = 0; j < 20; ++j)
        {
          _mm_store_pd(pmat+0,xmm0);
          _mm_store_pd(pmat+2,xmm0);
          _mm_store_pd(pmat+4,xmm0);
          _mm_store_pd(pmat+6,xmm0);
          _mm_store_pd(pmat+8,xmm0);
          _mm_store_pd(pmat+10,xmm0);
          _mm_store_pd(pmat+12,xmm0);
          _mm_store_pd(pmat+14,xmm0);
          _mm_store_pd(pmat+16,xmm0);
          _mm_store_pd(pmat+18,xmm0);
          pmat[j] = 1;
          pmat += 20;
        }
        continue;
      }

      rate = _mm_set1_pd(rates[n]);

      if (pinvar > PLL_MISC_EPSILON)
        xmm6 = _mm_set1_pd(1.0 - pinvar);

      for (k = 0; k < 20; k += 2)
      {
        xmm1 = _mm_load_pd(evals+k);

        /* scalar multiplication with rates */
        xmm4 = _mm_mul_pd(xmm1,rate);

        /* scalar multiplication with branch lengths */
        xmm5 = _mm_mul_pd(xmm4,brlen);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm5 = _mm_div_pd(xmm5,xmm6);
        }

        _mm_store_pd(expd+k,xmm5);
      }

      /* NOTE: in order to deal with numerical issues in cases when Qt -> 0, we
       * use a trick suggested by Ben Redelings and explained here:
       * https://github.com/xflouris/libpll/issues/129#issuecomment-304004005
       * In short, we use expm1() to compute (exp(Qt) - I), and then correct
       * for this by adding an identity matrix I in the very end */


      /* exponentiate eigenvalues */
      for (k = 0; k < 20; ++k)
        expd[k] = expm1(expd[k]);

      /* load expd */
      xmm0 = _mm_load_pd(expd+0);
      xmm1 = _mm_load_pd(expd+2);
      xmm2 = _mm_load_pd(expd+4);
      xmm3 = _mm_load_pd(expd+6);
      xmm4 = _mm_load_pd(expd+8);
      xmm5 = _mm_load_pd(expd+10);
      xmm6 = _mm_load_pd(expd+12);
      xmm7 = _mm_load_pd(expd+14);
      xmm8 = _mm_load_pd(expd+16);
      xmm9 = _mm_load_pd(expd+18);

      /* compute temp matrix */
      for (k = 0; k < 400; k += 20)
      {
        ymm0 = _mm_load_pd(inv_evecs+k+0);
        ymm1 = _mm_load_pd(inv_evecs+k+2);
        ymm2 = _mm_load_pd(inv_evecs+k+4);
        ymm3 = _mm_load_pd(inv_evecs+k+6);
        ymm4 = _mm_load_pd(inv_evecs+k+8);
        ymm5 = _mm_load_pd(inv_evecs+k+10);
        ymm6 = _mm_load_pd(inv_evecs+k+12);
        ymm7 = _mm_load_pd(inv_evecs+k+14);
        ymm8 = _mm_load_pd(inv_evecs+k+16);
        ymm9 = _mm_load_pd(inv_evecs+k+18);

        ymm0 = _mm_mul_pd(xmm0,ymm0);
        ymm1 = _mm_mul_pd(xmm1,ymm1);
        ymm2 = _mm_mul_pd(xmm2,ymm2);
        ymm3 = _mm_mul_pd(xmm3,ymm3);
        ymm4 = _mm_mul_pd(xmm4,ymm4);
        ymm5 = _mm_mul_pd(xmm5,ymm5);
        ymm6 = _mm_mul_pd(xmm6,ymm6);
        ymm7 = _mm_mul_pd(xmm7,ymm7);
        ymm8 = _mm_mul_pd(xmm8,ymm8);
        ymm9 = _mm_mul_pd(xmm9,ymm9);

        _mm_store_pd(temp+k+0,ymm0);
        _mm_store_pd(temp+k+2,ymm1);
        _mm_store_pd(temp+k+4,ymm2);
        _mm_store_pd(temp+k+6,ymm3);
        _mm_store_pd(temp+k+8,ymm4);
        _mm_store_pd(temp+k+10,ymm5);
        _mm_store_pd(temp+k+12,ymm6);
        _mm_store_pd(temp+k+14,ymm7);
        _mm_store_pd(temp+k+16,ymm8);
        _mm_store_pd(temp+k+18,ymm9);
      }

      for (j = 0; j < 400; j += 20)
      {
        xmm0 = _mm_load_pd(temp+j+0);
        xmm1 = _mm_load_pd(temp+j+2);
        xmm2 = _mm_load_pd(temp+j+4);
        xmm3 = _mm_load_pd(temp+j+6);
        xmm4 = _mm_load_pd(temp+j+8);
        xmm5 = _mm_load_pd(temp+j+10);
        xmm6 = _mm_load_pd(temp+j+12);
        xmm7 = _mm_load_pd(temp+j+14);
        xmm8 = _mm_load_pd(temp+j+16);
        xmm9 = _mm_load_pd(temp+j+18);

        /* process two rows at a time */
        for (k = 0; k < 400; k += 40)
        {
          /* row 0 */
          ONESTEP20(zmm0,tran+k+0);

          /* row 1 */
          ONESTEP20(zmm1,tran+k+20);

          zmm2 = _mm_hadd_pd(zmm0,zmm1);

          _mm_store_pd(pmat,zmm2);

          pmat += 2;
        }
      }

      /* add identity matrix */
      pmat -= 400;
      for (j = 0; j < 20; ++j)
      {
        pmat[j] += 1.0;
        pmat += 20;
      }
    }
  }

  pll_aligned_free(expd);
  pll_aligned_free(temp);

  for (i = 0; i < rate_cats; ++i)
    if (tran_evecs[i]) pll_aligned_free(tran_evecs[i]);

  free(tran_evecs);
  return PLL_SUCCESS;
}
