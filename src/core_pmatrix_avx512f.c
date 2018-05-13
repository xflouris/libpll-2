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

#include "pll.h"

#define STATES 20
#define STATES_PADDED 24

#define ONESTEP(x,baseptr)                                      \
            ymm0 = _mm512_load_pd(baseptr+0);                   \
            ymm1 = _mm512_load_pd(baseptr+8);                   \
            ymm2 = _mm512_load_pd(baseptr+16);                  \
                                                                \
            x = _mm512_mul_pd(xmm4,ymm0);                       \
            x = _mm512_fmadd_pd(xmm5,ymm1,x);                   \
            x = _mm512_fmadd_pd(xmm6,ymm2,x);                   \

PLL_EXPORT
int pll_core_update_pmatrix_20x20_avx512f(double ** pmatrix,
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

  expd = (double *)pll_aligned_alloc(STATES_PADDED*sizeof(double), PLL_ALIGNMENT_AVX512F);
  temp = (double *)pll_aligned_alloc(STATES_PADDED*STATES_PADDED*sizeof(double), PLL_ALIGNMENT_AVX512F);

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
      double * tran = (double *)pll_aligned_alloc(STATES_PADDED*STATES_PADDED*sizeof(double),
                                                  PLL_ALIGNMENT_AVX512F);
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
      for (i = 0; i < STATES; ++i)
      {
        for (j = 0; j < STATES; ++j)
          tran[i*STATES_PADDED+j] = evecs[j*STATES_PADDED+i];
      }
      
      /* update pointers and indicate that the eigen vector for the current
         rate matrix with index was updated */
      tran_evecs[index] = tran;
      transposed[index] = 1;
    }
  }
  free(transposed);

  __m512i permute_mask = _mm512_setr_epi64(0|4, 0|5, 0|6, 0|7, 8|0, 8|1, 8|2, 8|3);
  __m512i permute_mask_final_stage = _mm512_setr_epi64(0|2, 0|3, 8|0, 8|1, 0|6, 0|7, 8|4, 8|5);

  __m512d xmm0,xmm1,xmm2,xmm3,xmm4,xmm5,xmm6;
  __m512d ymm0,ymm1,ymm2,ymm3,ymm4,ymm5,ymm6;
  __m512d zmm0,zmm1,zmm2,zmm3,zmm4,zmm5,zmm6,zmm7;

  double * tran = NULL;
  for (i = 0; i < count; ++i)
  {
    assert(branch_lengths[i] >= 0);

    xmm3 = _mm512_set1_pd(branch_lengths[i]);
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
        xmm0 = _mm512_setzero_pd();
        for (j = 0; j < STATES; ++j)
        {
          _mm512_store_pd(pmat+0,xmm0);
          _mm512_store_pd(pmat+8,xmm0);
          _mm512_store_pd(pmat+16,xmm0);
          pmat[j] = 1;
          pmat += STATES_PADDED;
        }
        continue;
      }

      /* exponentiate eigenvalues */
      xmm2 = _mm512_set1_pd(rates[n]);

      if (pinvar > PLL_MISC_EPSILON)
        xmm6 = _mm512_set1_pd(1.0 - pinvar);

      for (k = 0; k < STATES_PADDED/ELEM_PER_AVX512_REGISTER; ++k)
      {
        xmm1 = _mm512_load_pd(evals+k*ELEM_PER_AVX512_REGISTER);

        /* scalar multiplication with rates */
        xmm4 = _mm512_mul_pd(xmm1,xmm2);

        /* scalar multiplication with branch lengths */
        xmm5 = _mm512_mul_pd(xmm4,xmm3);

        if (pinvar > PLL_MISC_EPSILON)
        {
          xmm5 = _mm512_div_pd(xmm5,xmm6);
        }

        _mm512_store_pd(expd+k*ELEM_PER_AVX512_REGISTER, xmm5);
      }

      for (k = 0; k < STATES; ++k)
        expd[k] = expm1(expd[k]);

      /* load expd */
      xmm4 = _mm512_load_pd(expd+0);
      xmm5 = _mm512_load_pd(expd+8);
      xmm6 = _mm512_load_pd(expd+16);

      /* compute temp matrix */
      for (k = 0; k < STATES*STATES_PADDED; k += STATES_PADDED)
      {
        ymm0 = _mm512_load_pd(inv_evecs+k+0);
        ymm1 = _mm512_load_pd(inv_evecs+k+8);
        ymm2 = _mm512_load_pd(inv_evecs+k+16);

        ymm4 = _mm512_mul_pd(xmm4,ymm0);
        ymm5 = _mm512_mul_pd(xmm5,ymm1);
        ymm6 = _mm512_mul_pd(xmm6,ymm2);

        _mm512_store_pd(temp+k+0,ymm4);
        _mm512_store_pd(temp+k+8,ymm5);
        _mm512_store_pd(temp+k+16,ymm6);
      }

      for (j = 0; j < STATES_PADDED*STATES; j += STATES_PADDED)
      {
        xmm4 = _mm512_load_pd(temp+j+0);
        xmm5 = _mm512_load_pd(temp+j+8);
        xmm6 = _mm512_load_pd(temp+j+16);

        /* process four rows at a time */
        for (k = 0; k < STATES_PADDED*STATES_PADDED; k += STATES_PADDED*8)
        {
          /* row 0 */
          ONESTEP(zmm0,tran+k+0);

          /* row 1 */
          ONESTEP(zmm1,tran+k+STATES_PADDED);

          /* row 2 */
          ONESTEP(zmm2,tran+k+STATES_PADDED*2);

          /* row 3 */
          ONESTEP(zmm3,tran+k+STATES_PADDED*3);

          /* row 4 */
          ONESTEP(zmm4,tran+k+STATES_PADDED*4);

          /* row 5 */
          ONESTEP(zmm5,tran+k+STATES_PADDED*5);

          /* row 6 */
          ONESTEP(zmm6,tran+k+STATES_PADDED*6);

          /* row 7 */
          ONESTEP(zmm7,tran+k+STATES_PADDED*7);

          /* create a vector with the sums of zmm0, zmm1, zmm2, zmm3, zmm4, zmm5, zmm6, zmm7 */
          ymm0 = _mm512_add_pd(_mm512_unpackhi_pd(zmm0,zmm1),
                               _mm512_unpacklo_pd(zmm0,zmm1));
          ymm1 = _mm512_add_pd(_mm512_unpackhi_pd(zmm2,zmm3),
                               _mm512_unpacklo_pd(zmm2,zmm3));
          ymm2 = _mm512_add_pd(_mm512_unpackhi_pd(zmm4,zmm5),
                               _mm512_unpacklo_pd(zmm4,zmm5));
          ymm3 = _mm512_add_pd(_mm512_unpackhi_pd(zmm6,zmm7),
                               _mm512_unpacklo_pd(zmm6,zmm7));

          zmm0 = _mm512_add_pd(_mm512_permutex2var_pd(ymm0, permute_mask, ymm2),
                               _mm512_mask_blend_pd(0xF0, ymm0, ymm2));

          zmm1 = _mm512_add_pd(_mm512_permutex2var_pd(ymm1, permute_mask, ymm3),
                               _mm512_mask_blend_pd(0xF0, ymm1, ymm3));


          zmm2 = _mm512_add_pd(_mm512_permutex2var_pd(zmm0,
                                                      permute_mask_final_stage,
                                                      zmm1),
                               _mm512_mask_blend_pd(0xCC, zmm0, zmm1));

          _mm512_store_pd(pmat, zmm2);

          pmat += ELEM_PER_AVX512_REGISTER;
        }
      }

      /* add identity matrix */
      pmat -= STATES_PADDED*STATES;
      for (j = 0; j < 20; ++j)
      {
        pmat[j] += 1.0;
        pmat += STATES_PADDED;
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
