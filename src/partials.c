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

static void case_tiptip(pll_partition_t * partition,
                        const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];

  double * parent_clv = pll_get_clv_writing(partition, op->parent_clv_index);
  assert(parent_clv);

  unsigned int * parent_scaler;
  unsigned int sites = partition->sites;

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  /* precompute lookup table */
  pll_core_create_lookup(partition->states,
                         partition->rate_cats,
                         partition->ttlookup,
                         left_matrix,
                         right_matrix,
                         partition->tipmap,
                         partition->maxstates,
                         partition->attributes);


  /* and update CLV at inner node */
  pll_core_update_partial_tt(partition->states,
                             sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             partition->tipchars[op->child1_clv_index],
                             partition->tipchars[op->child2_clv_index],
                             partition->tipmap,
                             partition->maxstates,
                             partition->ttlookup,
                             partition->attributes);
}

static void case_tipinner(pll_partition_t * partition,
                          const pll_operation_t * op)
{
  double * parent_clv = pll_get_clv_writing(partition, op->parent_clv_index);
  assert(parent_clv);
  unsigned int tip_clv_index;
  unsigned int inner_clv_index;
  unsigned int tip_matrix_index;
  unsigned int inner_matrix_index;
  unsigned int * right_scaler;
  unsigned int * parent_scaler;
  unsigned int sites = partition->sites;

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  /* find which of the two child nodes is the tip */
  if (op->child1_clv_index < partition->tips)
  {
    tip_clv_index = op->child1_clv_index;
    tip_matrix_index = op->child1_matrix_index;
    inner_clv_index = op->child2_clv_index;
    inner_matrix_index = op->child2_matrix_index;
    if (op->child2_scaler_index == PLL_SCALE_BUFFER_NONE)
      right_scaler = NULL;
    else
      right_scaler = partition->scale_buffer[op->child2_scaler_index];
  }
  else
  {
    tip_clv_index = op->child2_clv_index;
    tip_matrix_index = op->child2_matrix_index;
    inner_clv_index = op->child1_clv_index;
    inner_matrix_index = op->child1_matrix_index;
    if (op->child1_scaler_index == PLL_SCALE_BUFFER_NONE)
      right_scaler = NULL;
    else
      right_scaler = partition->scale_buffer[op->child1_scaler_index];
  }

  const double * inner_clv = pll_get_clv_reading(partition, inner_clv_index);
  assert(inner_clv);

  pll_core_update_partial_ti(partition->states,
                             sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             partition->tipchars[tip_clv_index],
                             inner_clv,
                             partition->pmatrix[tip_matrix_index],
                             partition->pmatrix[inner_matrix_index],
                             right_scaler,
                             partition->tipmap,
                             partition->maxstates,
                             partition->attributes);
}

static void case_innerinner(pll_partition_t * partition,
                            const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv      = pll_get_clv_writing(partition,
                                                 op->parent_clv_index);
  const double * left_clv  = pll_get_clv_reading(partition,
                                                 op->child1_clv_index);
  const double * right_clv = pll_get_clv_reading(partition,
                                                 op->child2_clv_index);
  assert(parent_clv);
  assert(left_clv);
  assert(right_clv);
  unsigned int * parent_scaler;
  unsigned int * left_scaler;
  unsigned int * right_scaler;
  unsigned int sites = partition->sites;

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  if (op->child1_scaler_index != PLL_SCALE_BUFFER_NONE)
    left_scaler = partition->scale_buffer[op->child1_scaler_index];
  else
    left_scaler = NULL;

  /* if child2 has a scaler add its values to the parent scaler */
  if (op->child2_scaler_index != PLL_SCALE_BUFFER_NONE)
    right_scaler = partition->scale_buffer[op->child2_scaler_index];
  else
    right_scaler = NULL;

  pll_core_update_partial_ii(partition->states,
                             sites,
                             partition->rate_cats,
                             parent_clv,
                             parent_scaler,
                             left_clv,
                             right_clv,
                             left_matrix,
                             right_matrix,
                             left_scaler,
                             right_scaler,
                             partition->attributes);
}

static void case_repeats(pll_partition_t * partition,
                            const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
  double * parent_clv      = pll_get_clv_writing(partition,
                                                 op->parent_clv_index);
  const double * left_clv  = pll_get_clv_reading(partition,
                                                 op->child1_clv_index);
  const double * right_clv = pll_get_clv_reading(partition,
                                                 op->child2_clv_index);
  assert(parent_clv);
  assert(left_clv);
  assert(right_clv);
  unsigned int * parent_scaler;
  unsigned int * left_scaler;
  unsigned int * right_scaler;
  unsigned int parent_sites         = pll_get_sites_number(partition,
                                                           op->parent_clv_index);
  const unsigned int * parent_id_site = pll_get_id_site(partition,
                                                        op->parent_clv_index);
  const unsigned int * left_site_id   = pll_get_site_id(partition,
                                                        op->child1_clv_index);
  const unsigned int * right_site_id  = pll_get_site_id(partition,
                                                        op->child2_clv_index);
  unsigned int left_sites          = pll_get_sites_number(partition,
                                                          op->child1_clv_index);
  unsigned int right_sites         = pll_get_sites_number(partition,
                                                          op->child2_clv_index);
  double * bclv_buffer             = partition->repeats
                                    ? partition->repeats->bclv_buffer : 0;
  unsigned int inv                 = left_sites < right_sites;


  /* get parent scaler */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];

  if (op->child1_scaler_index != PLL_SCALE_BUFFER_NONE)
    left_scaler = partition->scale_buffer[op->child1_scaler_index];
  else
    left_scaler = NULL;

  /* if child2 has a scaler add its values to the parent scaler */
  if (op->child2_scaler_index != PLL_SCALE_BUFFER_NONE)
    right_scaler = partition->scale_buffer[op->child2_scaler_index];
  else
    right_scaler = NULL;

  /* call the function with the shortest clv on the left */
  pll_core_update_partial_repeats(partition->states,
                                parent_sites,
                                inv  ? left_sites   : right_sites,
                                !inv ? left_sites   : right_sites,
                                partition->rate_cats,
                                parent_clv,
                                parent_scaler,
                                inv  ? left_clv     : right_clv,
                                !inv ? left_clv     : right_clv,
                                inv  ? left_matrix  : right_matrix,
                                !inv ? left_matrix  : right_matrix,
                                inv  ? left_scaler  : right_scaler,
                                !inv ? left_scaler  : right_scaler,
                                parent_id_site,
                                inv  ? left_site_id : right_site_id,
                                !inv ? left_site_id : right_site_id,
                                bclv_buffer,
                                partition->attributes);
}

PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count)
{
  pll_update_partials_rep(partition, operations, count, 1);
}


PLL_EXPORT void pll_update_partials_rep(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count,
                                    unsigned int update_repeats)
{
  unsigned int i;
  const pll_operation_t * op;

  for (i = 0; i < count; ++i)
  {
    op = &(operations[i]);

    // in the memory saver case, make sure to pin the clv we are about to
    // calculate before deciding which slot it should be in (avoids RC)
    if (pll_clv_manager_enabled(partition))
      pll_pin_clv( partition, op->parent_clv_index );

    if (pll_repeats_enabled(partition) && update_repeats) 
      pll_update_repeats(partition, op);

    if (pll_repeats_enabled(partition)
        && (partition->repeats->pernode_ids[op->child1_clv_index]
            ||  partition->repeats->pernode_ids[op->child2_clv_index]))
    {
      case_repeats(partition, op);
    }
    else if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      if ((op->child1_clv_index < partition->tips) &&
          (op->child2_clv_index < partition->tips))
      {
        /* tip-tip case */
        case_tiptip(partition,op);
      }
      else if ((operations[i].child1_clv_index < partition->tips) ||
               (operations[i].child2_clv_index < partition->tips))
      {
        /* tip-inner */
        case_tipinner(partition,op);
      }
      else
      {
        /* inner-inner */
        case_innerinner(partition,op);
      }
    }
    else
    {
      /* inner-inner */
      case_innerinner(partition,op);
    }

    // memory saver: reading of the children CLV definitely done now, so unpin
    if (pll_clv_manager_enabled(partition))
    {
      pll_unpin_clv( partition, op->child1_clv_index );
      pll_unpin_clv( partition, op->child2_clv_index );
    }
  }
}

