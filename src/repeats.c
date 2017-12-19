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

const unsigned int EMPTY_ELEMENT = (unsigned int) -1;


// map in charmap each char to a unique char identifier, according to map
static void repeats_fill_charmap(const pll_state_t *map, char *charmap)
{
  unsigned int i,j;
  char maxChar = 0;
  for (i = 0; i < PLL_ASCII_SIZE; ++i) 
  {
    for (j = 0; j < i; ++j) 
    {
      if (map[i] == map[j]) 
      {
        charmap[i] = charmap[j];
        break;
      }
    }
    if (!charmap[i]) 
      charmap[i] = ++maxChar;
  }
}

PLL_EXPORT int pll_repeats_enabled(const pll_partition_t *partition)
{
  return PLL_ATTRIB_SITE_REPEATS & partition->attributes;
}

PLL_EXPORT void pll_resize_repeats_lookup(pll_partition_t *partition, size_t size)
{
  if (!size)
    return;
  partition->repeats->lookup_buffer_size = size;
  free(partition->repeats->lookup_buffer);
  partition->repeats->lookup_buffer = 
    malloc(size * sizeof(unsigned int));
  memset(partition->repeats->lookup_buffer, EMPTY_ELEMENT, partition->repeats->lookup_buffer_size * sizeof(unsigned int));
}

PLL_EXPORT unsigned int pll_get_sites_number(const pll_partition_t * partition,
                                             unsigned int clv_index)
{
  unsigned int sites = partition->attributes & PLL_ATTRIB_SITE_REPEATS ?
      partition->repeats->pernode_ids[clv_index] : 0;
  sites = sites ? sites : partition->sites;
  sites += partition->asc_bias_alloc ? partition->states : 0;
  return sites;
}

PLL_EXPORT unsigned int pll_get_clv_size(const pll_partition_t * partition,
                                             unsigned int clv_index)
{
  return pll_get_sites_number(partition, clv_index) * 
    partition->states_padded * partition->rate_cats;
}

PLL_EXPORT unsigned int * pll_get_site_id(const pll_partition_t *partition,
                                                  unsigned int clv_index)
{
  unsigned int *site_id = 0;
  if (pll_repeats_enabled(partition) 
      && partition->repeats->pernode_ids[clv_index])
    site_id = partition->repeats->pernode_site_id[clv_index];
  return site_id;
}

PLL_EXPORT unsigned int * pll_get_id_site(const pll_partition_t *partition,
                                                  unsigned int clv_index)
{
  unsigned int *id_site = 0;
  if (pll_repeats_enabled(partition) 
      && partition->repeats->pernode_ids[clv_index])
    id_site = partition->repeats->pernode_id_site[clv_index];
  return id_site;
}

PLL_EXPORT unsigned int pll_default_enable_repeats(pll_partition_t *partition,
    unsigned int left_clv,
    unsigned int right_clv)
{
  pll_repeats_t * repeats = partition->repeats;
  unsigned int min_size = repeats->pernode_ids[left_clv] 
                          * repeats->pernode_ids[right_clv];
  return !(!min_size || (repeats->lookup_buffer_size <= min_size)
      || (repeats->pernode_ids[left_clv] > (partition->sites / 2))
      || (repeats->pernode_ids[right_clv] > (partition->sites / 2)));
}

PLL_EXPORT unsigned int pll_no_enable_repeats(pll_partition_t *partition,
    unsigned int left_clv,
    unsigned int right_clv)
{
  return 0;
}


PLL_EXPORT int pll_repeats_initialize(pll_partition_t *partition)
{
  int sites_alloc = partition->asc_additional_sites + partition->sites;
  unsigned int i;
  partition->repeats = malloc(sizeof(pll_repeats_t));
  if (!partition->repeats) 
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for repeats structure.");
    return PLL_FAILURE;
  }
  memset(partition->repeats, 0, sizeof(pll_repeats_t));
  pll_repeats_t *repeats = partition->repeats;
  repeats->enable_repeats = pll_default_enable_repeats;
  repeats->reallocate_repeats = pll_default_reallocate_repeats;
  repeats->pernode_site_id = calloc(partition->nodes, sizeof(unsigned int*));
  repeats->pernode_id_site = calloc(partition->nodes, sizeof(unsigned int*));
  if (!repeats->pernode_site_id || !repeats->pernode_id_site) 
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for repeats identifiers.");
    return PLL_FAILURE;
  }
  for (i = 0; i < partition->nodes; ++i) 
  {
    repeats->pernode_site_id[i] = calloc(sites_alloc, 
                                         sizeof(unsigned int));
    repeats->pernode_id_site[i] = calloc(sites_alloc, 
                                         sizeof(unsigned int));
    if (!repeats->pernode_site_id[i]) 
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg,
               200,
               "Unable to allocate enough memory for repeats identifiers.");
      return PLL_FAILURE;
    }
  }
  repeats->pernode_ids = calloc(partition->nodes, sizeof(unsigned int));
  repeats->perscale_ids = calloc(partition->scale_buffers, sizeof(unsigned int));
  repeats->pernode_allocated_clvs = 
    calloc(partition->nodes, sizeof(unsigned int));
  repeats->lookup_buffer = 0;
  repeats->lookup_buffer_size = 0;
  repeats->toclean_buffer = malloc(sites_alloc * sizeof(unsigned int));
  repeats->id_site_buffer = malloc(sites_alloc * sizeof(unsigned int));
  repeats->bclv_buffer = pll_aligned_alloc(sites_alloc 
      * partition->rate_cats * partition->states_padded
      * sizeof(double), partition->alignment);
  repeats->charmap = calloc(PLL_ASCII_SIZE, sizeof(char));
  if (!(repeats->pernode_ids
       && repeats->pernode_allocated_clvs && repeats->bclv_buffer
       && repeats->toclean_buffer && repeats->id_site_buffer 
       && repeats->charmap))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
          200,
          "Unable to allocate enough memory for one of the repeats buffer.");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_update_repeats_tips(pll_partition_t * partition,
                                  unsigned int tip_index,
                                  const pll_state_t * map,
                                  const char * sequence)
{
  if (!partition->repeats->lookup_buffer)
    pll_resize_repeats_lookup(partition, PLL_REPEATS_LOOKUP_SIZE);
   
  unsigned int s;
  pll_repeats_t * repeats = partition->repeats;
  unsigned int ** id_site = repeats->pernode_id_site;
  unsigned int additional_sites = 
    partition->asc_bias_alloc ? partition->states : 0;

  repeats_fill_charmap(map, repeats->charmap);
  repeats->pernode_ids[tip_index] = 0;
  unsigned int curr_id = 0;
  /* fill pernode_site_id */
  for (s = 0; s < partition->sites; ++s) 
  {
    unsigned int index_lookup = repeats->charmap[(int)sequence[s]];
    if (EMPTY_ELEMENT == repeats->lookup_buffer[index_lookup]) 
    {
      repeats->toclean_buffer[curr_id] = index_lookup;
      repeats->id_site_buffer[curr_id] = s;
      repeats->lookup_buffer[index_lookup] = curr_id++;
    }
    repeats->pernode_site_id[tip_index][s] = repeats->lookup_buffer[index_lookup];
  }
  unsigned int ids = curr_id;
  repeats->pernode_ids[tip_index] = ids;
  free(id_site[tip_index]);
  id_site[tip_index] = malloc(sizeof(unsigned int) 
      * (ids + additional_sites));
  for (s = 0; s < ids; ++s) 
  {
    id_site[tip_index][s] = repeats->id_site_buffer[s];
    repeats->lookup_buffer[repeats->toclean_buffer[s]] = EMPTY_ELEMENT;
  }
  for (s = 0; s < additional_sites; ++s) 
  {
    id_site[tip_index][ids + s] = partition->sites + s;
    repeats->pernode_site_id[tip_index][partition->sites + s] = ids + s;
  }
  unsigned int sizealloc = (ids + additional_sites) * partition->states_padded * 
                          partition->rate_cats * sizeof(double);
  free(partition->clv[tip_index]); 
  
  partition->clv[tip_index] = pll_aligned_alloc(sizealloc,
                                        partition->alignment);
  if (!partition->clv[tip_index]) 
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for repeats structure.");
    return PLL_FAILURE;
  }
  /* zero-out CLV vectors to avoid valgrind warnings when using odd number of
       states with vectorized code */
   repeats->pernode_allocated_clvs[tip_index] = ids;
  memset(partition->clv[tip_index],
          0,
          sizealloc);
  return PLL_SUCCESS;
}

PLL_EXPORT void pll_default_reallocate_repeats(pll_partition_t * partition,
                              unsigned int parent,
                              int scaler_index,
                              unsigned int sites_to_alloc)
{
  pll_repeats_t * repeats = partition->repeats;
  if (sites_to_alloc == repeats->pernode_allocated_clvs[parent]) 
    return;
  repeats->pernode_allocated_clvs[parent] = sites_to_alloc; 
  unsigned int ** id_site = repeats->pernode_id_site;
  // reallocate clvs
  pll_aligned_free(partition->clv[parent]);  
  partition->clv[parent] = pll_aligned_alloc(
      sites_to_alloc * partition->states_padded 
      * partition->rate_cats * sizeof(double), 
      partition->alignment);
  
  if (!partition->clv[parent]) 
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for repeats structure.");
    return;
  }
  // reallocate scales
  if (PLL_SCALE_BUFFER_NONE != scaler_index) 
  {
    unsigned int scaler_size = sites_to_alloc;
    if (partition->attributes & PLL_ATTRIB_RATE_SCALERS) 
      scaler_size *= partition->rate_cats;
    free(partition->scale_buffer[scaler_index]);
    partition->scale_buffer[scaler_index] = calloc(scaler_size, 
        sizeof(unsigned int));
  }
  // reallocate id to site lookup  
  free(id_site[parent]);
  id_site[parent] = malloc(sites_to_alloc * sizeof(unsigned int));
  // avoid valgrind errors
  memset(partition->clv[parent], 0, sites_to_alloc);
}

/* Fill the repeat structure in partition for the parent node of op */
PLL_EXPORT void pll_update_repeats(pll_partition_t * partition,
                    const pll_operation_t * op) 
{
  if (!partition->repeats->lookup_buffer)
    pll_resize_repeats_lookup(partition, PLL_REPEATS_LOOKUP_SIZE);
  
  pll_repeats_t * repeats = partition->repeats;
  unsigned int left = op->child1_clv_index;
  unsigned int right = op->child2_clv_index;
  unsigned int parent = op->parent_clv_index;
  unsigned int ** site_ids = repeats->pernode_site_id;
  unsigned int * site_id_parent = site_ids[parent];
  const unsigned int * site_id_left = site_ids[left];
  const unsigned int * site_id_right = site_ids[right];
  const unsigned int ids_left = repeats->pernode_ids[left];
  unsigned int ** id_site = repeats->pernode_id_site;
  unsigned int * toclean_buffer = repeats->toclean_buffer;
  unsigned int * id_site_buffer = repeats->id_site_buffer;
  unsigned int curr_id = 0;
  unsigned int additional_sites = partition->asc_bias_alloc ?
    partition->states : 0;
  unsigned int sites_to_alloc;
  unsigned int s;
  unsigned int ids = 0;
  // in case site repeats is activated but not used for this node
  if (!partition->repeats->enable_repeats(partition, left, right))
  {
    sites_to_alloc = partition->sites + additional_sites;
    repeats->pernode_ids[parent] = 0;
    if (op->parent_scaler_index != PLL_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = 0;
  } 
  else
  {
    // fill the parent repeats identifiers
    for (s = 0; s < partition->sites; ++s) 
    {
      unsigned int index_lookup = site_id_left[s] +
        site_id_right[s] * ids_left;
      unsigned int id = repeats->lookup_buffer[index_lookup];
      if (EMPTY_ELEMENT == id) 
      {
        toclean_buffer[curr_id] = index_lookup;
        id_site_buffer[curr_id] = s;
        id = curr_id;
        repeats->lookup_buffer[index_lookup] = curr_id++;
      }
      site_id_parent[s] = id;
    }
    ids = curr_id;
    for (s = 0; s < additional_sites; ++s) 
    {
      site_id_parent[s + partition->sites] = ids + s;
    }
    repeats->pernode_ids[parent] = ids;
    if (op->parent_scaler_index != PLL_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = ids;
    sites_to_alloc = ids + additional_sites;
  }

  repeats->reallocate_repeats( partition, 
                          op->parent_clv_index, 
                          op->parent_scaler_index, 
                          sites_to_alloc);

  // there is no repeats. Set pernode_ids to 0
  // to force the core functions not to use repeats
  if (sites_to_alloc >= partition->sites + additional_sites) {
    repeats->pernode_ids[parent] = 0;
    if (op->parent_scaler_index != PLL_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = 0;
  }

  // set id to site lookups
  for (s = 0; s < ids; ++s) 
  {
    id_site[parent][s] = id_site_buffer[s];
    repeats->lookup_buffer[toclean_buffer[s]] = EMPTY_ELEMENT;
  }
  for (s = 0; s < additional_sites; ++s) 
  {
    id_site[parent][s + ids] = partition->sites + s;
  }
}

PLL_EXPORT void pll_disable_bclv(pll_partition_t *partition)
{
  if (!pll_repeats_enabled(partition))
    return;
  pll_aligned_free(partition->repeats->bclv_buffer);
  partition->repeats->bclv_buffer = 0;
} 

PLL_EXPORT void pll_fill_parent_scaler_repeats(unsigned int sites,
                                       unsigned int * parent_scaler,
                                       const unsigned int * psites,
                                       const unsigned int * left_scaler,
                                       const unsigned int * lids,
                                       const unsigned int * right_scaler,
                                       const unsigned int * rids)
{
  // no repeats
  if (!lids && !rids) 
  {
    pll_fill_parent_scaler(sites, parent_scaler, left_scaler, right_scaler);
    return;
  }
  
  // no scalers
  if (!left_scaler && !right_scaler) 
  {
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
    return;
  }
  
  unsigned int i;
  if (!psites) {
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
    if (left_scaler) 
    {
      if (lids) 
      {
        for (i = 0; i < sites; ++i) 
          parent_scaler[i] += left_scaler[lids[i]];
      }
      else
      {
        for (i = 0; i < sites; ++i)
          parent_scaler[i] += left_scaler[i];
      }
    }
    if (right_scaler) 
    {
      if (rids) 
      {
        for (i = 0; i < sites; ++i) 
          parent_scaler[i] += right_scaler[rids[i]];
      }
      else
      {
        for (i = 0; i < sites; ++i)
          parent_scaler[i] += right_scaler[i];
      }
    }
  }
  else 
  {
    if (left_scaler && right_scaler) 
    {
      for (i = 0; i < sites; ++i) 
        parent_scaler[i] = left_scaler[lids[psites[i]]] + right_scaler[rids[psites[i]]];
    }
    else if (left_scaler) 
    {
      for (i = 0; i < sites; ++i) 
        parent_scaler[i] = left_scaler[lids[psites[i]]];
    }
    else 
    {
      for (i = 0; i < sites; ++i) 
        parent_scaler[i] = right_scaler[rids[psites[i]]];
    } 
  }
}

PLL_EXPORT void pll_fill_parent_scaler_repeats_per_rate(unsigned int sites,
                                       unsigned int rates,
                                       unsigned int * parent_scaler,
                                       const unsigned int * psites,
                                       const unsigned int * left_scaler,
                                       const unsigned int * lids,
                                       const unsigned int * right_scaler,
                                       const unsigned int * rids)
{
  unsigned int total_size = sites * rates;
  unsigned int cpy_size = rates * sizeof(unsigned int);
  unsigned int total_cpy_size = total_size * sizeof(unsigned int);
  // no repeats
  if (!lids && !rids) 
  {
    pll_fill_parent_scaler(total_size, parent_scaler, left_scaler, right_scaler);
    return;
  }
  
  // no scalers
  if (!left_scaler && !right_scaler) 
  {
    memset(parent_scaler, 0, total_cpy_size);
    return;
  }
  
  unsigned int i, j;
  if (!psites) {
    memset(parent_scaler, 0, total_cpy_size);
    if (left_scaler) 
    {
      if (lids) 
      {
        for (i = 0; i < sites; ++i) 
          memcpy(&parent_scaler[i * rates], &left_scaler[lids[i] * rates], cpy_size);
      }
      else
      {
        memcpy(parent_scaler, left_scaler, total_cpy_size);
      }
    }
    if (right_scaler) 
    {
      if (rids) 
      {
        for (i = 0; i < sites; ++i)
          for (j = 0; j < rates; ++j)
            parent_scaler[i * rates + j] += right_scaler[rids[i] * rates + j];
      }
      else
      {
        for (i = 0; i < total_size; ++i)
          parent_scaler[i] += right_scaler[i];
      }
    }
  }
  else 
  {
    if (left_scaler && right_scaler) 
    {
      for (i = 0; i < sites; ++i) 
          for (j = 0; j < rates; ++j)
            parent_scaler[i * rates + j] = left_scaler[lids[psites[i]] * rates + j] 
              + right_scaler[rids[psites[i]] * rates + j];
    }
    else if (left_scaler) 
    {
      for (i = 0; i < sites; ++i) 
        memcpy(&parent_scaler[i * rates], &left_scaler[lids[psites[i]] * rates], cpy_size);
    }
    else 
    {
      for (i = 0; i < sites; ++i) 
        memcpy(&parent_scaler[i * rates], &right_scaler[rids[psites[i]] * rates], cpy_size);
    } 
  }
}

