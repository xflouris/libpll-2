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

PLL_EXPORT double pll_core_root_loglikelihood(unsigned int states,
                                              unsigned int sites,
                                              unsigned int rate_cats,
                                              const double * clv,
                                              const unsigned int * scaler,
                                              double * const * frequencies,
                                              const double * rate_weights,
                                              const double * pattern_weights,
                                              const double * invar_proportion,
                                              const int * invar_indices,
                                              const unsigned int * freqs_indices,
                                              double * persite_lnl,
                                              unsigned int attrib)
{
  unsigned int i,j,k;
  double logl = 0;
  const double * freqs = NULL;

  double prop_invar = 0;

  double term, term_r;
  double site_lk, inv_site_lk;

  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_sse(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 invar_proportion,
                                                 invar_indices,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_sse(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             invar_proportion,
                                             invar_indices,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 invar_proportion,
                                                 invar_indices,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_avx(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             invar_proportion,
                                             invar_indices,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
  {
    if (states == 4)
    {
      return pll_core_root_loglikelihood_4x4_avx(sites,
                                                 rate_cats,
                                                 clv,
                                                 scaler,
                                                 frequencies,
                                                 rate_weights,
                                                 pattern_weights,
                                                 invar_proportion,
                                                 invar_indices,
                                                 freqs_indices,
                                                 persite_lnl);
    }
    else
    {
      return pll_core_root_loglikelihood_avx2(states,
                                             sites,
                                             rate_cats,
                                             clv,
                                             scaler,
                                             frequencies,
                                             rate_weights,
                                             pattern_weights,
                                             invar_proportion,
                                             invar_indices,
                                             freqs_indices,
                                             persite_lnl);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif


  /* iterate through sites */
  for (i = 0; i < sites; ++i)
  {
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clv[k] * freqs[k];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[j]] : 0;
      if (prop_invar > 0)
      {
        inv_site_lk = (invar_indices[i] == -1) ?
                           0 : freqs[invar_indices[i]];
        term += rate_weights[j] * (term_r * (1 - prop_invar) +
                                   inv_site_lk*prop_invar);
      }
      else
      {
        term += term_r * rate_weights[j];
      }

      clv += states_padded;
    }

    site_lk = term;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(site_lk);
    if (scaler && scaler[i])
      site_lk += scaler[i] * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = site_lk;

    logl += site_lk;
  }
  return logl;
}

PLL_EXPORT double pll_core_root_loglikelihood_repeats_generic(unsigned int states,
                                                              unsigned int sites,
                                                              unsigned int rate_cats,
                                                              const double * clv,
                                                              const unsigned int * site_id,
                                                              const unsigned int * scaler,
                                                              double * const * frequencies,
                                                              const double * rate_weights,
                                                              const double * pattern_weights,
                                                              const double * invar_proportion,
                                                              const int * invar_indices,
                                                              const unsigned int * freqs_indices,
                                                              double * persite_lnl)
{
  unsigned int i,j,k;
  double logl = 0;
  const double * freqs = NULL;

  double prop_invar = 0;

  double term, term_r;
  double site_lk, inv_site_lk;

  unsigned int states_padded = states;
  unsigned int span = states_padded * rate_cats;

  /* iterate through sites */
  for (i = 0; i < sites; ++i)
  {
    unsigned int id = PLL_GET_ID(site_id, i);
    const double *clvp = &clv[id * span];
    term = 0;
    for (j = 0; j < rate_cats; ++j)
    {
      freqs = frequencies[freqs_indices[j]];
      term_r = 0;
      for (k = 0; k < states; ++k)
      {
        term_r += clvp[k] * freqs[k];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[j]] : 0;
      if (prop_invar > 0)
      {
        inv_site_lk = (invar_indices[i] == -1) ?
                           0 : freqs[invar_indices[i]];
        term += rate_weights[j] * (term_r * (1 - prop_invar) +
                                   inv_site_lk*prop_invar);
      }
      else
      {
        term += term_r * rate_weights[j];
      }

      clvp += states_padded;
    }

    site_lk = term;

    /* compute site log-likelihood and scale if necessary */
    site_lk = log(site_lk);
    if (scaler && scaler[id])
      site_lk += scaler[id] * log(PLL_SCALE_THRESHOLD);

    site_lk *= pattern_weights[i];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[i] = site_lk;

    logl += site_lk;
  }
  return logl;
}

PLL_EXPORT double pll_core_root_loglikelihood_repeats(unsigned int states,
                                                      unsigned int sites,
                                                      unsigned int rate_cats,
                                                      const double * clv,
                                                      const unsigned int * site_id,
                                                      const unsigned int * scaler,
                                                      double * const * frequencies,
                                                      const double * rate_weights,
                                                      const double * pattern_weights,
                                                      const double * invar_proportion,
                                                      const int * invar_indices,
                                                      const unsigned int * freqs_indices,
                                                      double * persite_lnl,
                                                      unsigned int attrib)
{
  double (*core_root_loglikelihood) (unsigned int states,
                                     unsigned int sites,
                                     unsigned int rate_cats,
                                     const double * clv,
                                     const unsigned int * site_id,
                                     const unsigned int * scaler,
                                     double * const * frequencies,
                                     const double * rate_weights,
                                     const double * pattern_weights,
                                     const double * invar_proportion,
                                     const int * invar_indices,
                                     const unsigned int * freqs_indices,
                                     double * persite_lnl) = 0x0;

  core_root_loglikelihood = pll_core_root_loglikelihood_repeats_generic;
#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
  {
    core_root_loglikelihood = pll_core_root_loglikelihood_repeats_avx;
  }
#endif
#ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
  {
    core_root_loglikelihood = pll_core_root_loglikelihood_repeats_sse;
  }
#endif
#ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
  {
    core_root_loglikelihood = pll_core_root_loglikelihood_repeats_avx2;
    // TODO call 4x4 avx (not avx2) functions when implemented
  }
#endif
    return core_root_loglikelihood(states,
                                  sites,
                                  rate_cats,
                                  clv,
                                  site_id,
                                  scaler,
                                  frequencies,
                                  rate_weights,
                                  pattern_weights,
                                  invar_proportion,
                                  invar_indices,
                                  freqs_indices,
                                  persite_lnl);
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti_4x4(unsigned int sites,
                                          unsigned int rate_cats,
                                          const double * parent_clv,
                                          const unsigned int * parent_scaler,
                                          const unsigned char * tipchars,
                                          const double * pmatrix,
                                          double * const * frequencies,
                                          const double * rate_weights,
                                          const double * pattern_weights,
                                          const double * invar_proportion,
                                          const int * invar_indices,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl,
                                          unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;

  const double * clvp = parent_clv;
  double prop_invar = 0;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r, termb, terminv;
  double site_lk, inv_site_lk;

  unsigned int cstate;

  unsigned int states = 4;
  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
  {
    return pll_core_edge_loglikelihood_ti_4x4_sse(sites,
                                                  rate_cats,
                                                  parent_clv,
                                                  parent_scaler,
                                                  tipchars,
                                                  pmatrix,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  invar_proportion,
                                                  invar_indices,
                                                  freqs_indices,
                                                  persite_lnl,
                                                  attrib);
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
  {
    return pll_core_edge_loglikelihood_ti_4x4_avx(sites,
                                                  rate_cats,
                                                  parent_clv,
                                                  parent_scaler,
                                                  tipchars,
                                                  pmatrix,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  invar_proportion,
                                                  invar_indices,
                                                  freqs_indices,
                                                  persite_lnl,
                                                  attrib);
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
  {
    return pll_core_edge_loglikelihood_ti_4x4_avx(sites,
                                                  rate_cats,
                                                  parent_clv,
                                                  parent_scaler,
                                                  tipchars,
                                                  pmatrix,
                                                  frequencies,
                                                  rate_weights,
                                                  pattern_weights,
                                                  invar_proportion,
                                                  invar_indices,
                                                  freqs_indices,
                                                  persite_lnl,
                                                  attrib);
  }
  #endif

  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || invar_proportion)
  {
    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat = pmatrix;
    terma = 0;
    terminv = 0;

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      prop_invar = invar_proportion[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        cstate = (unsigned int) (*tipchars);
        for (k = 0; k < states; ++k)
        {
          if (cstate & 1)
            termb += pmat[k];
          cstate >>= 1;
        }
        terma_r += clvp[j] * freqs[j] * termb;
        pmat += states_padded;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      if (prop_invar > 0)
      {
        terma += rate_weights[i] * terma_r * (1. - prop_invar);
        if (invar_indices[n] != -1)
        {
          freqs = frequencies[freqs_indices[i]];
          inv_site_lk = freqs[invar_indices[n]];
          terminv += rate_weights[i] * inv_site_lk * prop_invar;
        }
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvp += states_padded;
    }

    /* compute site log-likelihood and scale if necessary */
    if (site_scalings)
    {
      if (terminv > 0.)
      {
        /* IMPORTANT: undoing the scaling for non-invariant likelihood term only! */
        unsigned int capped_scalings = PLL_MIN(site_scalings, PLL_SCALE_RATE_MAXDIFF);
        double scale_factor = scale_minlh[capped_scalings-1];
        site_lk = log(terma * scale_factor + terminv);
      }
      else
      {
        site_lk = log(terma);
        site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);
      }
    }
    else
    {
      site_lk = log(terma + terminv);
    }

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;

    tipchars++;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_ti(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      const double * parent_clv,
                                      const unsigned int * parent_scaler,
                                      const unsigned char * tipchars,
                                      const pll_state_t * tipmap,
                                      unsigned int tipmap_size,
                                      const double * pmatrix,
                                      double * const * frequencies,
                                      const double * rate_weights,
                                      const double * pattern_weights,
                                      const double * invar_proportion,
                                      const int * invar_indices,
                                      const unsigned int * freqs_indices,
                                      double * persite_lnl,
                                      unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;

  const double * clvp = parent_clv;
  double prop_invar = 0;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r, termb, terminv;
  double site_lk, inv_site_lk;

  pll_state_t cstate;

  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ti_4x4_sse(sites,
                                                    rate_cats,
                                                    parent_clv,
                                                    parent_scaler,
                                                    tipchars,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ti_sse(states,
                                                sites,
                                                rate_cats,
                                                parent_clv,
                                                parent_scaler,
                                                tipchars,
                                                tipmap,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ti_4x4_avx(sites,
                                                    rate_cats,
                                                    parent_clv,
                                                    parent_scaler,
                                                    tipchars,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else if (states == 20)
    {
      return pll_core_edge_loglikelihood_ti_20x20_avx(sites,
                                                      rate_cats,
                                                      parent_clv,
                                                      parent_scaler,
                                                      tipchars,
                                                      tipmap,
                                                      tipmap_size,
                                                      pmatrix,
                                                      frequencies,
                                                      rate_weights,
                                                      pattern_weights,
                                                      invar_proportion,
                                                      invar_indices,
                                                      freqs_indices,
                                                      persite_lnl,
                                                      attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ti_avx(states,
                                                sites,
                                                rate_cats,
                                                parent_clv,
                                                parent_scaler,
                                                tipchars,
                                                tipmap,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ti_4x4_avx(sites,
                                                    rate_cats,
                                                    parent_clv,
                                                    parent_scaler,
                                                    tipchars,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else if (states == 20)
    {
      return pll_core_edge_loglikelihood_ti_20x20_avx2(sites,
                                                      rate_cats,
                                                      parent_clv,
                                                      parent_scaler,
                                                      tipchars,
                                                      tipmap,
                                                      tipmap_size,
                                                      pmatrix,
                                                      frequencies,
                                                      rate_weights,
                                                      pattern_weights,
                                                      invar_proportion,
                                                      invar_indices,
                                                      freqs_indices,
                                                      persite_lnl,
                                                      attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ti_avx(states,
                                                sites,
                                                rate_cats,
                                                parent_clv,
                                                parent_scaler,
                                                tipchars,
                                                tipmap,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || invar_proportion)
  {
    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[n] : 0;
    }

    pmat = pmatrix;
    terma = 0;
    terminv = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      prop_invar = invar_proportion[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        cstate = tipmap[(unsigned int)(*tipchars)];
        for (k = 0; k < states; ++k)
        {
          if (cstate & 1)
            termb += pmat[k];
          cstate >>= 1;
        }
        terma_r += clvp[j] * freqs[j] * termb;
        pmat += states_padded;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      if (prop_invar > 0)
      {
        terma += rate_weights[i] * terma_r * (1. - prop_invar);
        if (invar_indices[n] != -1)
        {
          freqs = frequencies[freqs_indices[i]];
          inv_site_lk = freqs[invar_indices[n]];
          terminv += rate_weights[i] * inv_site_lk * prop_invar;
        }
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvp += states_padded;
    }

    /* compute site log-likelihood and scale if necessary */
    if (site_scalings)
    {
      if (terminv > 0.)
      {
        /* IMPORTANT: undoing the scaling for non-variant likelihood term only! */
        unsigned int capped_scalings = PLL_MIN(site_scalings, PLL_SCALE_RATE_MAXDIFF);
        double scale_factor = scale_minlh[capped_scalings-1];
        site_lk = log(terma * scale_factor + terminv);
      }
      else
      {
        site_lk = log(terma);
        site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);
      }
    }
    else
    {
      site_lk = log(terma + terminv);
    }

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl)
      persite_lnl[n] = site_lk;

    logl += site_lk;

    tipchars++;
  }

  if (rate_scalings)
    free(rate_scalings);

  return logl;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_repeats(unsigned int states,
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
                                           const double * pattern_weights,
                                           const double * invar_proportion,
                                           const int * invar_indices,
                                           const unsigned int * freqs_indices,
                                           double * persite_lnl,
                                           const unsigned int * parent_site_id,
                                           const unsigned int * child_site_id,
                                           double * bclv,
                                           unsigned int attrib)
{
  double (*core_edge_loglikelihood)(unsigned int states,
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
                                    const double * pattern_weights,
                                    const double * invar_proportion,
                                    const int * invar_indices,
                                    const unsigned int * freqs_indices,
                                    double * persite_lnl,
                                    const unsigned int * parent_site_id,
                                    const unsigned int * child_site_id,
                                    double * bclv,
                                    unsigned int attrib) = 0x0;

  unsigned int use_bclv = bclv && (child_sites < (sites * 2) / 3);
  core_edge_loglikelihood = pll_core_edge_loglikelihood_repeats_generic;
  
#ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX &&  PLL_STAT(avx_present))
  {
    core_edge_loglikelihood = pll_core_edge_loglikelihood_repeats_generic_avx;
    if (states == 4)
    {
      if (use_bclv)
        core_edge_loglikelihood =  pll_core_edge_loglikelihood_repeatsbclv_4x4_avx;
      else
        core_edge_loglikelihood = pll_core_edge_loglikelihood_repeats_4x4_avx;
    }
  } 
#endif
#ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE &&  PLL_STAT(sse3_present))
  {
    core_edge_loglikelihood = pll_core_edge_loglikelihood_repeats_generic_sse;
  }
#endif
#ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 &&  PLL_STAT(avx2_present))
  {
    core_edge_loglikelihood = pll_core_edge_loglikelihood_repeats_generic_avx2;
  }
#endif
  return core_edge_loglikelihood(states,
                                 sites,
                                 child_sites,
                                 rate_cats,
                                 parent_clv,
                                 parent_scaler,
                                 child_clv,
                                 child_scaler,
                                 pmatrix,
                                 frequencies,
                                 rate_weights,
                                 pattern_weights,
                                 invar_proportion,
                                 invar_indices,
                                 freqs_indices,
                                 persite_lnl, 
                                 parent_site_id,
                                 child_site_id,
                                 bclv,
                                 attrib);
  return 0.0;
}

PLL_EXPORT
double pll_core_edge_loglikelihood_repeats_generic(unsigned int states,
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
                                                   const double * pattern_weights,
                                                   const double * invar_proportion,
                                                   const int * invar_indices,
                                                   const unsigned int * freqs_indices,
                                                   double * persite_lnl,
                                                   const unsigned int * parent_site_id,
                                                   const unsigned int * child_site_id,
                                                   double * bclv,
                                                   unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;

  double prop_invar = 0;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r, termb,terminv;
  double site_lk, inv_site_lk;

  unsigned int span = states * rate_cats;
  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || invar_proportion)
  {
    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    unsigned int pid = PLL_GET_ID(parent_site_id, n);
    unsigned int cid = PLL_GET_ID(child_site_id, n);
    const double *clvp = &parent_clv[pid * span];
    const double *clvc = &child_clv[cid * span];
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[pid*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[cid*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = PLL_MIN(rate_scalings[i] - site_scalings,
                                   PLL_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings =  (parent_scaler) ? parent_scaler[pid] : 0;
      site_scalings += (child_scaler) ? child_scaler[cid] : 0;
    }

    pmat = pmatrix;
    terma = 0;
    terminv = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
        {
          termb += pmat[k] * clvc[k];
        }

        terma_r += clvp[j] * freqs[j] * termb;
        pmat += states;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        terma += rate_weights[i] * terma_r * (1. - prop_invar);
        if (invar_indices[n] != -1)
        {
          freqs = frequencies[freqs_indices[i]];
          inv_site_lk = freqs[invar_indices[n]];
          terminv += rate_weights[i] * inv_site_lk * prop_invar;
        }
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvp += states;
      clvc += states;
    }

    /* compute site log-likelihood and scale if necessary */
    if (site_scalings)
    {
      if (terminv > 0.)
      {
        /* IMPORTANT: undoing the scaling for non-variant likelihood term only! */
        unsigned int capped_scalings = PLL_MIN(site_scalings, PLL_SCALE_RATE_MAXDIFF);
        double scale_factor = scale_minlh[capped_scalings-1];
        site_lk = log(terma * scale_factor + terminv);
      }
      else
      {
        site_lk = log(terma);
        site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);
      }
    }
    else
    {
      site_lk = log(terma + terminv);
    }

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
double pll_core_edge_loglikelihood_ii(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      const double * parent_clv,
                                      const unsigned int * parent_scaler,
                                      const double * child_clv,
                                      const unsigned int * child_scaler,
                                      const double * pmatrix,
                                      double * const * frequencies,
                                      const double * rate_weights,
                                      const double * pattern_weights,
                                      const double * invar_proportion,
                                      const int * invar_indices,
                                      const unsigned int * freqs_indices,
                                      double * persite_lnl,
                                      unsigned int attrib)
{
  unsigned int n,i,j,k;
  double logl = 0;

  const double * clvp = parent_clv;
  const double * clvc = child_clv;
  double prop_invar = 0;
  const double * pmat;
  const double * freqs = NULL;

  double terma, terma_r, termb,terminv;
  double site_lk, inv_site_lk;

  /* TODO: We need states_padded in the AVX/SSE implementations
  */
  unsigned int states_padded = states;

  #ifdef HAVE_SSE3
  if (attrib & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ii_4x4_sse(sites,
                                                    rate_cats,
                                                    clvp,
                                                    parent_scaler,
                                                    clvc,
                                                    child_scaler,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ii_sse(states,
                                                sites,
                                                rate_cats,
                                                clvp,
                                                parent_scaler,
                                                clvc,
                                                child_scaler,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+1) & 0xFFFFFFFE;
  }
  #endif
  #ifdef HAVE_AVX
  if (attrib & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ii_4x4_avx(sites,
                                                    rate_cats,
                                                    clvp,
                                                    parent_scaler,
                                                    clvc,
                                                    child_scaler,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ii_avx(states,
                                                sites,
                                                rate_cats,
                                                clvp,
                                                parent_scaler,
                                                clvc,
                                                child_scaler,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif
  #ifdef HAVE_AVX2
  if (attrib & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present))
  {
    if (states == 4)
    {
      return pll_core_edge_loglikelihood_ii_4x4_avx(sites,
                                                    rate_cats,
                                                    clvp,
                                                    parent_scaler,
                                                    clvc,
                                                    child_scaler,
                                                    pmatrix,
                                                    frequencies,
                                                    rate_weights,
                                                    pattern_weights,
                                                    invar_proportion,
                                                    invar_indices,
                                                    freqs_indices,
                                                    persite_lnl,
                                                    attrib);
    }
    else
    {
      return pll_core_edge_loglikelihood_ii_avx2(states,
                                                sites,
                                                rate_cats,
                                                clvp,
                                                parent_scaler,
                                                clvc,
                                                child_scaler,
                                                pmatrix,
                                                frequencies,
                                                rate_weights,
                                                pattern_weights,
                                                invar_proportion,
                                                invar_indices,
                                                freqs_indices,
                                                persite_lnl,
                                                attrib);
    }
    /* this line is never called, but should we disable the else case above,
       then states_padded must be set to this value */
    states_padded = (states+3) & 0xFFFFFFFC;
  }
  #endif

  unsigned int site_scalings;
  unsigned int * rate_scalings = NULL;
  int per_rate_scaling = (attrib & PLL_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[PLL_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || invar_proportion)
  {
    double scale_factor = 1.0;
    for (i = 0; i < PLL_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= PLL_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int*) calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate space for rate scalers.");
      return -INFINITY;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = (parent_scaler) ? parent_scaler[n*rate_cats+i] : 0;
        rate_scalings[i] += (child_scaler) ? child_scaler[n*rate_cats+i] : 0;
        if (rate_scalings[i] < site_scalings)
          site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
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

    pmat = pmatrix;
    terma = 0;
    terminv = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      freqs = frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
        {
          termb += pmat[k] * clvc[k];
        }

        terma_r += clvp[j] * freqs[j] * termb;
        pmat += states_padded;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i]-1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        terma += rate_weights[i] * terma_r * (1. - prop_invar);
        if (invar_indices[n] != -1)
        {
          freqs = frequencies[freqs_indices[i]];
          inv_site_lk = freqs[invar_indices[n]];
          terminv += rate_weights[i] * inv_site_lk * prop_invar;
        }
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvp += states_padded;
      clvc += states_padded;
    }

    /* compute site log-likelihood and scale if necessary */
    if (site_scalings)
    {
      if (terminv > 0.)
      {
        /* IMPORTANT: undoing the scaling for non-variant likelihood term only! */
        unsigned int capped_scalings = PLL_MIN(site_scalings, PLL_SCALE_RATE_MAXDIFF);
        double scale_factor = scale_minlh[capped_scalings-1];
        site_lk = log(terma * scale_factor + terminv);
      }
      else
      {
        site_lk = log(terma);
        site_lk += site_scalings * log(PLL_SCALE_THRESHOLD);
      }
    }
    else
    {
      site_lk = log(terma + terminv);
    }

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

