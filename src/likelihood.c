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

#include "pll.h"

static double compute_asc_bias_correction(double logl_base,
                                          unsigned int sum_w,
                                          unsigned int sum_w_inv,
                                          int asc_bias_type)
{
  double logl_correction = 0.0;
  switch (asc_bias_type)
  {
    case PLL_ATTRIB_AB_LEWIS:
      logl_correction = -(sum_w*log(1 - logl_base));
      break;
    case PLL_ATTRIB_AB_STAMATAKIS:
      /* no need to add anything here */
      logl_correction = logl_base;
      break;
    case PLL_ATTRIB_AB_FELSENSTEIN:
      logl_correction = sum_w_inv * log(logl_base);
      break;
    default:
      pll_errno = PLL_ERROR_AB_INVALIDMETHOD;
      snprintf(pll_errmsg, 200, "Illegal ascertainment bias algorithm");
      return -INFINITY;
  }
  return logl_correction;
}

static double root_loglikelihood_asc_bias(pll_partition_t * partition,
                                          unsigned int sites,
                                          const double * clv,
                                          unsigned int * scaler,
                                          const unsigned int * freqs_indices)
{
   unsigned int i,j,k;
   double logl = 0;
   double term_r, term;
   double site_lk;

   const double * freqs = NULL;
   unsigned int states = partition->states;
   unsigned int states_padded = partition->states_padded;
   unsigned int scale_factors;
   double * pattern_weights = partition->pattern_weights;
   double * rate_weights = partition->rate_weights;

   double logl_correction = 0;
   unsigned int sum_w_inv = 0;
   int asc_bias_type = partition->attributes & PLL_ATTRIB_AB_MASK;

   /* point clvp to state sites */
   clv += sites * partition->rate_cats * partition->states_padded;

   /* 1. compute per-site logl for each state */
   for (i = 0; i < states; ++i)
   {
     term = 0;
     for (j = 0; j < partition->rate_cats; ++j)
     {
       freqs = partition->frequencies[freqs_indices[j]];
       term_r = 0;
       for (k = 0; k < states; ++k)
       {
         term_r += clv[k] * freqs[k];
       }
       term += term_r * rate_weights[j];
       clv += states_padded;
     }

     /* count number of scaling factors to acount for */
     scale_factors = scaler ? scaler[sites + i] : 0;

     sum_w_inv += pattern_weights[sites + i];
     if (asc_bias_type == PLL_ATTRIB_AB_STAMATAKIS)
     {
       /* 2a. site_lk is the lnl weighted by the number of occurences */
       site_lk = log(term) * partition->pattern_weights[sites + i];
       if (scale_factors)
         site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);
     }
     else
     {
       /* 2b. site_lk is the actual likelihood */
         site_lk = term * pow(PLL_SCALE_THRESHOLD, scale_factors);
     }

     /* logl_correction is the sum of weighted log likelihoods if
        asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
     logl_correction += site_lk;
   }

   /* 3. apply correction to the lnl score */
   logl += compute_asc_bias_correction(logl_correction,
                                       partition->pattern_weight_sum,
                                       sum_w_inv,
                                       asc_bias_type);

   return logl;
}

PLL_EXPORT double pll_compute_root_loglikelihood(pll_partition_t * partition,
                                                 unsigned int clv_index,
                                                 int scaler_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl)
{
  double logl = 0;
  unsigned int * scaler;
  unsigned int identifiers;
  /* get scaler array if specified */
  if (scaler_index == PLL_SCALE_BUFFER_NONE)
    scaler = NULL;
  else
    scaler = partition->scale_buffer[scaler_index];

  /* compute log-likelihood via the core function */
  if (pll_repeats_enabled(partition) &&
      partition->repeats->pernode_ids[clv_index]) 
  {
    logl = pll_core_root_loglikelihood_repeats(partition->states,
                                     partition->sites,
                                     partition->rate_cats,
                                     partition->clv[clv_index],
                                     partition->repeats->pernode_site_id[clv_index],
                                     scaler,
                                     partition->frequencies,
                                     partition->rate_weights,
                                     partition->pattern_weights,
                                     partition->prop_invar,
                                     partition->invariant,
                                     freqs_indices,
                                     persite_lnl,
                                     partition->attributes);
  }
  else
  {
    /* compute log-likelihood via the core function */
    logl = pll_core_root_loglikelihood(partition->states,
                                     partition->sites,
                                     partition->rate_cats,
                                     partition->clv[clv_index],
                                     scaler,
                                     partition->frequencies,
                                     partition->rate_weights,
                                     partition->pattern_weights,
                                     partition->prop_invar,
                                     partition->invariant,
                                     freqs_indices,
                                     persite_lnl,
                                     partition->attributes);
  }

  /* ascertainment bias correction */
  if (partition->attributes & PLL_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    identifiers = pll_get_sites_number(partition, clv_index) - partition->rate_cats;
    logl += root_loglikelihood_asc_bias(partition,
                                        identifiers,
                                        partition->clv[clv_index],
                                        scaler,
                                        freqs_indices);
  }

  return logl;
}

static double edge_loglikelihood_asc_bias_ti(pll_partition_t * partition,
                                             unsigned int parent_clv_index,
                                             unsigned int * parent_scaler,
                                             unsigned int matrix_index,
                                             const unsigned int * freqs_indices)
{
  unsigned int n,i,j;
  double logl = 0;
  double terma, terma_r;
  double site_lk;

  const double * clvp = partition->clv[parent_clv_index];
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int scale_factors;
  double * pattern_weights = partition->pattern_weights;
  double * rate_weights = partition->rate_weights;

  double logl_correction = 0;
  unsigned int sum_w_inv = 0;
  int asc_bias_type = partition->attributes & PLL_ATTRIB_AB_MASK;

  /* point clvp to state sites */
  clvp += partition->sites * partition->rate_cats * partition->states_padded;
  pattern_weights += partition->sites;
  if (parent_scaler)
    parent_scaler += partition->sites;

  /* 1. compute per-site logl for each state */
  for (n = 0; n < partition->states; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      freqs = partition->frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        terma_r += clvp[j] * freqs[j] * pmatrix[n];
        pmatrix += states_padded;
      }

      terma += terma_r * rate_weights[i];
      clvp += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

    sum_w_inv += pattern_weights[n];
    if (asc_bias_type == PLL_ATTRIB_AB_STAMATAKIS)
    {
      /* 2a. site_lk is the lnl weighted by the number of occurences */
      site_lk = log(terma) * pattern_weights[n];
      if (scale_factors)
        site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);
    }
    else
    {
      /* 2b. site_lk is the actual likelihood */
      site_lk = terma * pow(PLL_SCALE_THRESHOLD, scale_factors);
    }

    /* logl_correction is the sum of weighted log likelihoods if
       asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
    logl_correction += site_lk;
  }

  /* 3. apply correction to the lnl score */
  logl += compute_asc_bias_correction(logl_correction,
                                      partition->pattern_weight_sum,
                                      sum_w_inv,
                                      asc_bias_type);
  return logl;
}

static double edge_loglikelihood_tipinner(pll_partition_t * partition,
                                          unsigned int parent_clv_index,
                                          int parent_scaler_index,
                                          unsigned int child_clv_index,
                                          unsigned int matrix_index,
                                          const unsigned int * freqs_indices,
                                          double * persite_lnl)
{
  double logl = 0;
  unsigned int states = partition->states;

  unsigned int * parent_scaler;

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  if (states == 4)
  {
    logl = pll_core_edge_loglikelihood_ti_4x4(partition->sites,
                                              partition->rate_cats,
                                              partition->clv[parent_clv_index],
                                              parent_scaler,
                                              partition->tipchars[child_clv_index],
                                              partition->pmatrix[matrix_index],
                                              partition->frequencies,
                                              partition->rate_weights,
                                              partition->pattern_weights,
                                              partition->prop_invar,
                                              partition->invariant,
                                              freqs_indices,
                                              persite_lnl,
                                              partition->attributes);
  }
  else
  {
    logl = pll_core_edge_loglikelihood_ti(partition->states,
                                          partition->sites,
                                          partition->rate_cats,
                                          partition->clv[parent_clv_index],
                                          parent_scaler,
                                          partition->tipchars[child_clv_index],
                                          partition->tipmap,
                                          partition->maxstates,
                                          partition->pmatrix[matrix_index],
                                          partition->frequencies,
                                          partition->rate_weights,
                                          partition->pattern_weights,
                                          partition->prop_invar,
                                          partition->invariant,
                                          freqs_indices,
                                          persite_lnl,
                                          partition->attributes);
  }

  /* ascertainment bias correction */
  if (partition->attributes & PLL_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ti(partition,
                                           parent_clv_index,
                                           parent_scaler,
                                           matrix_index,
                                           freqs_indices);
  }

  return logl;
}

static double edge_loglikelihood_asc_bias_ii(pll_partition_t * partition,
                                             const double * clvp,
                                             unsigned int * parent_scaler,
                                             unsigned int parent_sites,
                                             const double * clvc,
                                             unsigned int * child_scaler,
                                             unsigned int child_sites,
                                             unsigned int matrix_index,
                                             const unsigned int * freqs_indices)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk;

  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int scale_factors;
  double * pattern_weights = partition->pattern_weights;
  double * rate_weights = partition->rate_weights;

  /* point clv, clvc, scalers and pattern weights to state sites */
  unsigned int parent_offset = parent_sites * 
                        partition->rate_cats *
                        partition->states_padded;
  unsigned int child_offset = child_sites * 
                        partition->rate_cats *
                        partition->states_padded;

  pattern_weights += partition->sites;
  clvp += parent_offset;
  clvc += child_offset;
  if (parent_scaler)
    parent_scaler += parent_sites;
  if (child_scaler)
    child_scaler += child_sites;

  double logl_correction = 0;
  unsigned int sum_w_inv = 0;
  int asc_bias_type = partition->attributes & PLL_ATTRIB_AB_MASK;

  /* 1. compute per-site logl for each state */
  for (n = 0; n < partition->states; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      freqs = partition->frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
        {
          termb += pmatrix[k] * clvc[k];
        }
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states_padded;
      }

      terma += terma_r * rate_weights[i];
      clvp += states_padded;
      clvc += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    sum_w_inv += pattern_weights[n];
    if (asc_bias_type == PLL_ATTRIB_AB_STAMATAKIS)
    {
      /* 2a. site_lk is the lnl weighted by the number of occurences */
      site_lk = log(terma) * pattern_weights[n];
      if (scale_factors)
        site_lk += scale_factors * log(PLL_SCALE_THRESHOLD);
    }
    else
    {
      /* 2b. site_lk is the actual likelihood */
      site_lk = terma * pow(PLL_SCALE_THRESHOLD, scale_factors);
    }

    /* logl_correction is the sum of weighted log likelihoods if
       asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
    logl_correction += site_lk;
  }

  /* 3. apply correction to the lnl score */
  logl += compute_asc_bias_correction(logl_correction,
                                      partition->pattern_weight_sum,
                                      sum_w_inv,
                                      asc_bias_type);

  return logl;
}

static double edge_loglikelihood(pll_partition_t * partition,
                                 unsigned int parent_clv_index,
                                 int parent_scaler_index,
                                 unsigned int child_clv_index,
                                 int child_scaler_index,
                                 unsigned int matrix_index,
                                 const unsigned int * freqs_indices,
                                 double * persite_lnl)
{
  double logl = 0;

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];

  unsigned int * parent_scaler;
  unsigned int * child_scaler;

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  /* compute log-likelihood via the core function */
  logl = pll_core_edge_loglikelihood_ii(partition->states,
                                        partition->sites,
                                        partition->rate_cats,
                                        clvp,
                                        parent_scaler,
                                        clvc,
                                        child_scaler,
                                        partition->pmatrix[matrix_index],
                                        partition->frequencies,
                                        partition->rate_weights,
                                        partition->pattern_weights,
                                        partition->prop_invar,
                                        partition->invariant,
                                        freqs_indices,
                                        persite_lnl,
                                        partition->attributes);

  /* ascertainment bias correction */
  if (partition->attributes & PLL_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ii(partition,
                                           clvp,
                                           parent_scaler,
                                           partition->sites,
                                           clvc,
                                           child_scaler,
                                           partition->sites,
                                           matrix_index,
                                           freqs_indices);
  }
  return logl;
}

static double edge_loglikelihood_repeats(pll_partition_t * partition,
                                 unsigned int parent_clv_index,
                                 int parent_scaler_index,
                                 unsigned int child_clv_index,
                                 int child_scaler_index,
                                 unsigned int matrix_index,
                                 const unsigned int * freqs_indices,
                                 double * persite_lnl)
{
  double logl = 0;

  const double * clvp = partition->clv[parent_clv_index];
  const double * clvc = partition->clv[child_clv_index];

  const unsigned int * parent_site_id =
    pll_get_site_id(partition, parent_clv_index);
  const unsigned int * child_site_id = 
    pll_get_site_id(partition, child_clv_index);
  unsigned int parent_sites = 
    pll_get_sites_number(partition, parent_clv_index);
  unsigned int child_sites = 
    pll_get_sites_number(partition, child_clv_index);
  unsigned int inv = parent_sites > child_sites;
  
  unsigned int * parent_scaler;
  unsigned int * child_scaler;

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  /* compute log-likelihood via the core function */
  logl = pll_core_edge_loglikelihood_repeats(partition->states,
                                        partition->sites,
                                        !inv ? parent_sites : child_sites, 
                                        partition->rate_cats,
                                        inv ? clvp : clvc,
                                        inv ? parent_scaler : child_scaler,
                                        !inv ? clvp : clvc,
                                        !inv ? parent_scaler : child_scaler,
                                        partition->pmatrix[matrix_index],
                                        partition->frequencies,
                                        partition->rate_weights,
                                        partition->pattern_weights,
                                        partition->prop_invar,
                                        partition->invariant,
                                        freqs_indices,
                                        persite_lnl, 
                                        inv ? parent_site_id : child_site_id,
                                        !inv ? parent_site_id : child_site_id,
                                        partition->repeats->bclv_buffer,
                                        partition->attributes);

  /* ascertainment bias correction */
  if (partition->attributes & PLL_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ii(partition,
                                           clvp,
                                           parent_scaler,
                                           parent_sites - partition->states,
                                           clvc,
                                           child_scaler,
                                           child_sites - partition->states,
                                           matrix_index,
                                           freqs_indices);
  }
  return logl;
}



PLL_EXPORT double pll_compute_edge_loglikelihood(pll_partition_t * partition,
                                                 unsigned int parent_clv_index,
                                                 int parent_scaler_index,
                                                 unsigned int child_clv_index,
                                                 int child_scaler_index,
                                                 unsigned int matrix_index,
                                                 const unsigned int * freqs_indices,
                                                 double * persite_lnl)
{
  double logl;
  
  if (pll_repeats_enabled(partition) 
      && (partition->repeats->pernode_ids[parent_clv_index] 
          || partition->repeats->pernode_ids[child_clv_index]))
  {
    return edge_loglikelihood_repeats(partition,
                            parent_clv_index,
                            parent_scaler_index,
                            child_clv_index,
                            child_scaler_index,
                            matrix_index,
                            freqs_indices,
                            persite_lnl);

  }

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    if ((parent_clv_index < partition->tips) ||
        (child_clv_index < partition->tips))
      return edge_loglikelihood_tipinner(partition,
                                        (parent_clv_index < partition->tips) ?
                                            child_clv_index : parent_clv_index,
                                        (parent_clv_index < partition->tips) ?
                                            child_scaler_index : parent_scaler_index,
                                        (parent_clv_index < partition->tips) ?
                                            parent_clv_index : child_clv_index,
                                        matrix_index,
                                        freqs_indices,
                                        persite_lnl);

  logl = edge_loglikelihood(partition,
                            parent_clv_index,
                            parent_scaler_index,
                            child_clv_index,
                            child_scaler_index,
                            matrix_index,
                            freqs_indices,
                            persite_lnl);

  return logl;
}


PLL_EXPORT int pll_compute_node_ancestral_extbuf(pll_partition_t * partition,
                                                 unsigned int node_clv_index,
                                                 int node_scaler_index,
                                                 unsigned int other_clv_index,
                                                 int other_scaler_index,
                                                 unsigned int pmatrix_index,
                                                 const unsigned int * freqs_indices,
                                                 double * ancestral,
                                                 double * temp_clv,
                                                 unsigned int * temp_scaler,
                                                 double * ident_pmat)
{
  if (!partition || !ancestral)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200,
             "Parameter value is NULL!");
    return PLL_FAILURE;
  }

  if (!temp_clv || !temp_scaler || !ident_pmat)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200,
             "NULL buffer pointer");
    return PLL_FAILURE;
  }

  if (pll_repeats_enabled(partition))
  {
    pll_errno = PLL_ERROR_EINVAL;
    snprintf(pll_errmsg, 200,
             "Site repeats are not compatible with ancestral state reconstruction!");
    return PLL_FAILURE;
  }

  unsigned int n, i, j;
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites = partition->sites;
  unsigned int rate_cats = partition->rate_cats;

  const double * pmat = partition->pmatrix[pmatrix_index];

  const double * node_clv = partition->clv[node_clv_index];

  unsigned int * node_scaler = (node_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                            NULL : partition->scale_buffer[node_scaler_index];

  if (other_clv_index < partition->tips &&
      (partition->attributes & PLL_ATTRIB_PATTERN_TIP))
  {
    pll_core_update_partial_ti(states,
                               sites,
                               rate_cats,
                               temp_clv,
                               temp_scaler,
                               partition->tipchars[other_clv_index],
                               node_clv,
                               pmat,
                               ident_pmat,
                               node_scaler,
                               partition->tipmap,
                               partition->maxstates,
                               partition->attributes);
  }
  else
  {
    const double * other_clv = partition->clv[other_clv_index];
    unsigned int * other_scaler = (other_scaler_index == PLL_SCALE_BUFFER_NONE) ?
                              NULL : partition->scale_buffer[other_scaler_index];

    pll_core_update_partial_ii(states,
                               sites,
                               rate_cats,
                               temp_clv,
                               temp_scaler,
                               node_clv,
                               other_clv,
                               ident_pmat,
                               pmat,
                               node_scaler,
                               other_scaler,
                               partition->attributes);
  }

  double * clvp = temp_clv;
  double * ancp = ancestral;

  memset(ancestral, 0, sites * states * sizeof(double));

  for (n = 0; n < sites; ++n)
  {
    for (i = 0; i < rate_cats; ++i)
    {
      double * freqs = partition->frequencies[freqs_indices[i]];
      double rate_weight = partition->rate_weights[i];

      for (j = 0; j < states; ++j)
      {
        ancp[j] += clvp[j] * freqs[j] * rate_weight;
      }

      clvp += states_padded;
    }

    // normalize probs
    double sum = 0;
    for (j = 0; j < states; ++j)
      sum += ancp[j];
    for (j = 0; j < states; ++j)
      ancp[j] /= sum;

    ancp += states;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_compute_node_ancestral(pll_partition_t * partition,
                                          unsigned int node_clv_index,
                                          int node_scaler_index,
                                          unsigned int other_clv_index,
                                          int other_scaler_index,
                                          unsigned int matrix_index,
                                          const unsigned int * freqs_indices,
                                          double * ancestral)
{
  int retval = PLL_FAILURE;
  unsigned int i, j, k;

  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites = partition->sites;
  unsigned int rate_cats = partition->rate_cats;

  unsigned int clv_size = pll_get_clv_size(partition, node_clv_index) * sizeof(double);
  double * temp_clv = (double *) pll_aligned_alloc(clv_size, partition->alignment);
  unsigned int scaler_size = ((partition->attributes & PLL_ATTRIB_RATE_SCALERS) ?
                             sites*rate_cats : sites) * sizeof(unsigned int);
  unsigned int * temp_scaler = (unsigned int *) pll_aligned_alloc(scaler_size, partition->alignment);
  unsigned int displacement = (states_padded - states) * (states_padded) * sizeof(double);
  unsigned int pmat_size = states_padded * states * rate_cats * sizeof(double) + displacement;
  double * ident_pmat = (double *) pll_aligned_alloc(pmat_size, partition->alignment);

  if (!temp_clv || !temp_scaler || !ident_pmat)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory");
    goto cleanup;
  }

  // init identity matrix
  memset(ident_pmat, 0, pmat_size);
  double * pmat = ident_pmat;
  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states; ++j)
    {
      for (k = 0; k < states_padded; ++k)
        pmat[j*states_padded + k] = (j == k) ? 1 : 0;
    }

    pmat +=  states*states_padded;
  }

  retval = pll_compute_node_ancestral_extbuf(partition,
                                           node_clv_index,
                                           node_scaler_index,
                                           other_clv_index,
                                           other_scaler_index,
                                           matrix_index,
                                           freqs_indices,
                                           ancestral,
                                           temp_clv,
                                           temp_scaler,
                                           ident_pmat);

cleanup:
  pll_aligned_free(temp_clv);
  pll_aligned_free(temp_scaler);
  pll_aligned_free(ident_pmat);

  return retval;
}
