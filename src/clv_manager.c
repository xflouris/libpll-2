/*
    Copyright (C) 2020 Pierre Barbera, HITS gGmbH

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

    Contact: Alexandros Stamatakis <alexandros.stamatakis@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

/**
 * Get the address of the specified CLV for reading. If partition is
 * memory-managed and the specified CLV is not currently in memory, return NULL.
 *
 * @param  partition the partition the CLV belongs to
 * @param  clv_index the index of the CLV
 * @return           address of the CLV, or NULL if CLV doesn't reside in memory
 */
PLL_EXPORT const double * pll_get_clv_reading(
                                        const pll_partition_t * const partition,
                                        const unsigned int clv_index)
{
  assert(partition);

  if ( !(partition->attributes & PLL_ATTRIB_LIMIT_MEMORY) )
  {
    return partition->clv[clv_index];
  }
  else
  {
    pll_clv_manager_t * clv_man = partition->clv_man;
    assert(clv_man);

    unsigned int slot = clv_man->slot_of_clvid[clv_index];

    return (slot != PLL_CLV_NODE_UNPINNED) ? partition->clv[slot] : NULL;
  }
}

/**
 * Get the address of the specified CLV for writing.
 *
 * If partition is memory-managed and the specified CLV is not currently in
 * memory, first tries to pick a CLV slot to overwrite and return its address.
 * If no overwritable slot exists, apply replacement strategy to get a slot.
 *
 * @param  partition the partition the CLV belongs to
 * @param  clv_index the index of the CLV
 * @return           address of the CLV
 */
PLL_EXPORT double * pll_get_clv_writing(pll_partition_t * const partition,
                                        const unsigned int clv_index)
{
  assert(partition);

  if ( !(partition->attributes & PLL_ATTRIB_LIMIT_MEMORY) )
  {
    return partition->clv[clv_index];
  }
  else
  {
    pll_clv_manager_t * clv_man = partition->clv_man;
    assert(clv_man);

    unsigned int slot = clv_man->slot_of_clvid[clv_index];
    // check if clv is pinned
    if (slot != PLL_CLV_NODE_UNPINNED)
    {
      return partition->clv[slot];
    }
    else
    {
      // CLV not slotted, check if any slots are available to be overwritten
      if (!clv_man->unpinnable->empty)
      {
        return partition->clv[pll_uint_stack_pop(clv_man->unpinnable)];
      }
      else
      {
        // no slots available, need to run the replacement strategy
        return clv_man->replacer->replace(partition, clv_man);
      }
    }
  }
}

static void* alloc_and_set(const size_t n, const size_t size, const int val)
{
  void* data = malloc(n * size);

  if (!data)
  {
    return PLL_FAILURE;
  }

  memset(data, val, n * size);

  return data;
}

void dealloc_clv_manager(pll_clv_manager_t * clv_man)
{
  if (clv_man)
  {
    free(clv_man->clvid_of_slot);
    free(clv_man->slot_of_clvid);
    free(clv_man->unpinnable);
    free(clv_man->replacer);
    free(clv_man);
  }
}

/**
 * Initializes the memory manager according to the selected maximum number of
 * CLVs that should be in memory concurrently.
 *
 * @return                 PLL_SUCCESS or PLL_FAILURE
 */
PLL_EXPORT int pll_clv_manager_init(pll_partition_t * const partition,
                                    const size_t concurrent_clvs,
                                    pll_clv_manager_strategy_t * strategy)
{
  assert(partition);
  assert(strategy);

  const size_t addressable_clvs = partition->clv_buffers;
  // const size_t concurrent_clvs = partition->clv_buffers;

  assert(concurrent_clvs <= addressable_clvs);

  // set the pll attribute
  // partition->attributes |= PLL_ATTRIB_LIMIT_MEMORY;
  // this happens at partition creation time now

  if (pll_repeats_enabled(partition))
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg,
             200,
             "Memory management not yet possible together with site repeats");
    return PLL_FAILURE;
  }


  /**
   * =============== MEMORY ALLOCATION ===============
   */

  // alloc the manager struct
  pll_clv_manager_t * clv_man = partition->clv_man =
                        (pll_clv_manager_t *)malloc(sizeof(pll_clv_manager_t));
  if (!clv_man)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for clv_manager.");
    return PLL_FAILURE;
  }

  // set member fields to defaults
  clv_man->size = concurrent_clvs;

  clv_man->clvid_of_slot = (unsigned int *)alloc_and_set(concurrent_clvs,
                                                         sizeof(unsigned int),
                                                         PLL_CLV_SLOT_UNUSED);
  if (!clv_man->clvid_of_slot)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for clvid_of_slot.");
    return PLL_FAILURE;
  }


  clv_man->slot_of_clvid = (unsigned int *)alloc_and_set(addressable_clvs,
                                                         sizeof(unsigned int),
                                                         PLL_CLV_NODE_UNPINNED);
  if (!clv_man->slot_of_clvid)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for slot_of_clvid.");
    return PLL_FAILURE;
  }

  clv_man->unpinnable = pll_uint_stack_create(concurrent_clvs);

  if (!clv_man->unpinnable)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for unpinnable.");
    return PLL_FAILURE;
  }

  // initialize the unpinnable array: all slots are ready to be used
  for( size_t i = 0; i < concurrent_clvs; ++i )
  {
    pll_uint_stack_push(clv_man->unpinnable, i);
  }

  clv_man->replacer = strategy;

  // alloc the CLVs of the partition!
  if(!alloc_clvs(partition, concurrent_clvs))
  {
    dealloc_clv_manager(clv_man);
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}
