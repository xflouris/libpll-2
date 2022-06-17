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

// Helper functions
static bool is_addressable( const pll_partition_t * const partition,
                            const unsigned int clv_index)
{
  assert(partition->clv_man);
  return (clv_index >= partition->clv_man->addressable_begin)
      && (clv_index < partition->clv_man->addressable_end);
}

static bool is_tip( const pll_partition_t * const partition,
                    const unsigned int clv_index )
{
  return clv_index < partition->tips;
}

///////////////////////
// GETTERS / SETTERS //
///////////////////////

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

  if ( is_tip(partition, clv_index) || !pll_clv_manager_enabled(partition) )
  {
    return partition->clv[clv_index];
  }
  else
  {
    if (!is_addressable(partition, clv_index))
    {
      pll_errno = PLL_ERROR_CLV_MANAGER_FAIL;
      snprintf(pll_errmsg,
               200,
               "clv_index not in addressable range (tipchar?).");
      return NULL;
    }

    pll_clv_manager_t * clv_man = partition->clv_man;
    assert(clv_man);
    const size_t offset = clv_man->addressable_begin;

    unsigned int slot = clv_man->slot_of_clvid[clv_index];

    return (slot != PLL_CLV_CLV_UNSLOTTED) 
            ? partition->clv[offset + slot] : NULL;
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

  if ( is_tip(partition, clv_index) || !pll_clv_manager_enabled(partition))
  {
    return partition->clv[clv_index];
  }
  else
  {
    if (!is_addressable(partition, clv_index))
    {
      pll_errno = PLL_ERROR_CLV_MANAGER_FAIL;
      snprintf(pll_errmsg,
               200,
               "clv_index not in addressable range (tipchar?).");
      return NULL;
    }

    pll_clv_manager_t * clv_man = partition->clv_man;
    assert(clv_man);

    unsigned int slot = clv_man->slot_of_clvid[clv_index];
    const size_t offset = clv_man->addressable_begin;

    // check if clv is slotted
    if (slot != PLL_CLV_CLV_UNSLOTTED)
    {
      return partition->clv[offset + slot];
    }
    else
    {
      // CLV not slotted, check if any slots are completely unused so far
      if (pll_uint_stack_empty(clv_man->unused_slots))
      {
        // no slots available, need to run the replacement strategy
        slot = clv_man->strat_replace(clv_man);
      }
      else
      {
        // get a free slot
        slot = pll_uint_stack_pop(clv_man->unused_slots);
      }

      // finally, associate the requested CLV with that slot
      pll_clv_manager_update_slot(clv_man, slot, clv_index);

      return partition->clv[offset + slot];
    }
  }
}

/**
 * Check if given clv is slotted
 *
 * @param  partition the partition the clv_index belongs to
 * @param  clv_index the index
 * @return           whether the clv_index is slotted
 */
PLL_EXPORT bool pll_clv_is_slotted( const pll_partition_t * const partition,
                                    const unsigned int clv_index)
{
  if ( !(partition->attributes & PLL_ATTRIB_LIMIT_MEMORY) )
  {
    return true;
  }
  else
  {
    assert(partition->clv_man);
    return  partition->clv_man->slot_of_clvid[clv_index] != PLL_CLV_CLV_UNSLOTTED;
  }
}

////////////////////////////
// HANDLING PINNING       //
////////////////////////////

static int pinning_impl(pll_clv_manager_t* const clv_man,
                        const unsigned int clv_index,
                        const bool pin)
{
  if (!clv_man)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg,
             200,
             "Pinning only makes sense when using memory saver.");
    return PLL_FAILURE;
  }

  // pin / unpin
  const bool was_pinned = clv_man->is_pinned[ clv_index ];
  clv_man->is_pinned[ clv_index ] = pin;

  // keep track of number of concurrently pinned clvs
  clv_man->num_pinned += -(int)was_pinned + pin;

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_pin_clv( pll_partition_t * const partition,
                            const unsigned int clv_index)
{
  return pinning_impl( partition->clv_man, clv_index, true );
}

PLL_EXPORT int pll_unpin_clv( pll_partition_t * const partition,
                              const unsigned int clv_index)
{
  return pinning_impl( partition->clv_man, clv_index, false );
}


////////////////////////////
// CREATION / DESTRUCTION //
////////////////////////////

void dealloc_clv_manager(pll_clv_manager_t * clv_man)
{
  if (clv_man)
  {
    free(clv_man->clvid_of_slot);
    free(clv_man->slot_of_clvid);
    free(clv_man->is_pinned);
    pll_uint_stack_destroy(clv_man->unused_slots);
    if(clv_man->strat_data_dealloc)
      clv_man->strat_data_dealloc(clv_man->repl_strat_data);
    if (clv_man->repl_strat_data)
      free(clv_man->repl_strat_data);
    free(clv_man);
  }
}

static void* alloc_and_set(const size_t n, const size_t size, const int val)
{
  void* data = malloc(n * size);

  if (!data)
    return PLL_FAILURE;

  memset(data, val, n * size);

  return data;
}

/**
 * Initializes the memory manager according to the selected maximum number of
 * CLVs that should be in memory concurrently.
 * 
 * @param  partition       partition for which the manager should be initialized
 * @param  concurrent_clvs number of slots/clvs to be allocated
 * @param  cb_replace      replacement strategy replace callback (NULL = MRC)
 * @param  cb_update       replacement strategy slot update callback (NULL = MRC)
 * @param  cb_dealloc      replacement strategy deallocation callback (NULL = MRC)
 * @return                 PLL_{SUCCESS|FAILURE}
 */
PLL_EXPORT int pll_clv_manager_init(pll_partition_t * const partition,
                                    const size_t concurrent_clvs,
                                    pll_clv_manager_replace_cb cb_replace,
                                    pll_clv_manager_update_cb cb_update,
                                    pll_clv_manager_dealloc_cb cb_dealloc)
{
  assert(partition);

  const size_t addressable_clvs = partition->nodes;
  const bool pattern_tip = (partition->attributes & PLL_ATTRIB_PATTERN_TIP);

  const size_t partition_clv_bufs = partition->clv_buffers
                                  + (pattern_tip ? 0 : partition->tips);

  assert(concurrent_clvs <= partition_clv_bufs);

  // if (pll_repeats_enabled(partition))
  // {
  //   pll_errno = PLL_ERROR_PARAM_INVALID;
  //   snprintf(pll_errmsg,
  //            200,
  //            "Memory management not yet possible together with site repeats");
  //   return PLL_FAILURE;
  // }

  ///////////////////////
  // MEMORY ALLOCATION //
  ///////////////////////

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
  clv_man->slottable_size     = concurrent_clvs;
  clv_man->addressable_end    = addressable_clvs;
  clv_man->addressable_begin  = pattern_tip ? partition->tips : 0;
  // if replacement func was null, use default
  clv_man->strat_replace      = cb_replace ? cb_replace : &MRC_replace_cb;
  clv_man->strat_update_slot  = cb_update ? cb_update : &MRC_update_slot_cb;
  clv_man->strat_data_dealloc = cb_dealloc ? cb_dealloc : &MRC_dealloc_cb;

  // alloc everything else

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
                                                         PLL_CLV_CLV_UNSLOTTED);
  if (!clv_man->slot_of_clvid)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for slot_of_clvid.");
    return PLL_FAILURE;
  }

  clv_man->is_pinned = (bool *)alloc_and_set(addressable_clvs,
                                             sizeof(bool),
                                             false);
  if (!clv_man->is_pinned)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for is_pinned.");
    return PLL_FAILURE;
  }

  clv_man->num_pinned = 0;

  clv_man->unused_slots = pll_uint_stack_create(concurrent_clvs);

  if (!clv_man->unused_slots)
  {
    dealloc_clv_manager(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for unused_slots.");
    return PLL_FAILURE;
  }

  // initialize the unused_slots array: all slots are ready to be used
  for( size_t i = 0; i < concurrent_clvs; ++i )
  {
    int ret = pll_uint_stack_push(clv_man->unused_slots, i);
    assert(ret == PLL_SUCCESS);
  }

  // alloc the CLVs of the partition!
  if(!alloc_clvs(partition, concurrent_clvs))
  {
    dealloc_clv_manager(clv_man);
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT bool pll_clv_manager_enabled(pll_partition_t const * const partition)
{
  return (partition->attributes & PLL_ATTRIB_LIMIT_MEMORY);
}

/**
 * Update the tracked information about a given slot with a given clv_index.
 *
 * Calls update function of the replacement strategy if one is set.
 *
 * @param clv_man   pointer to the clv manager struct to be updated
 * @param slot      slot to be updated
 * @param clv_index clv_index that will occupy the given slot
 */
void pll_clv_manager_update_slot(pll_clv_manager_t * clv_man,
                                 const unsigned int slot,
                                 const unsigned int clv_index)
{
  // unslot the old index
  const unsigned int old_clvid = clv_man->clvid_of_slot[slot];
  if (old_clvid != PLL_CLV_SLOT_UNUSED)
    clv_man->slot_of_clvid[old_clvid] = PLL_CLV_CLV_UNSLOTTED;
  // update to the new
  clv_man->clvid_of_slot[slot]      = clv_index;
  clv_man->slot_of_clvid[clv_index] = slot;
  // if custom replacement strat update function is set, call it
  if (clv_man->strat_update_slot)
    clv_man->strat_update_slot(clv_man, slot, clv_index);
}

//////////////////////////
// REPLACEMENT STRATEGY //
//////////////////////////

typedef struct mrc_data
{
  unsigned int* cost_of_clvid;
    // <addressable_size> entries, provides the approx. recomputation cost value
    // of each clvid
  unsigned int* cost_of_slot;
    // <slottable_size> entries, provides the approx. recomputation cost value
    // of each slot. Has to be updated during replacement
} mrc_data_t;

/**
 * MRC specific callback to determine which slot should be overwritten/replaced
 *
 * @param  clv_man pointer to the clv manager
 * @return         slot to be overwritten, according to this strategy
 */
unsigned int MRC_replace_cb(pll_clv_manager_t* clv_man)
{
  assert(clv_man);

  // this function should only be called if there is no free slot
  assert(pll_uint_stack_empty(clv_man->unused_slots));

  // get a pointer to the data and cast it correctly
  // in this case our "data" is the mrc_data struct, through which we can access
  // the cost of each slot and clvid
  mrc_data_t * mrc = (mrc_data_t *) clv_man->repl_strat_data;
  assert(mrc);

  const size_t slots = clv_man->slottable_size;

  // get the cheapest non-pinned slot
  size_t cheapest_slot        = -1u;
  unsigned int cheapest_cost  = -1u;
  for (size_t slot_id = 0; slot_id < slots; ++slot_id)
  {
    assert(clv_man->clvid_of_slot[slot_id] != PLL_CLV_SLOT_UNUSED);

    if (!clv_man->is_pinned[clv_man->clvid_of_slot[slot_id]]
     && (mrc->cost_of_slot[slot_id] < cheapest_cost))
    {
      cheapest_slot = slot_id;
      cheapest_cost = mrc->cost_of_slot[slot_id];
    }
  }

  assert(cheapest_slot != -1u);
  assert(cheapest_cost != -1u);

  return cheapest_slot;
}

/**
 * MRC strategy specific slot update function
 *
 * @param clv_man
 * @param slot      slot to be updated
 * @param clv_index index to update to
 */
void MRC_update_slot_cb(pll_clv_manager_t * clv_man,
                        const unsigned int slot,
                        const unsigned int clv_index)
{
  mrc_data_t * mrc = (mrc_data_t *) clv_man->repl_strat_data;
  assert(mrc);

  mrc->cost_of_slot[slot] = mrc->cost_of_clvid[clv_index];
}

/**
 * Deallocates the custom data stored for the MRC strategy
 *
 * @param data the data pointer
 */
void MRC_dealloc_cb(void* data)
{
  mrc_data_t* mrcd = (mrc_data_t*) data;
  free(mrcd->cost_of_clvid);
  free(mrcd->cost_of_slot);
}

static int MRC_strategy_init(pll_clv_manager_t * clv_man)
{
  const size_t addr_size  = clv_man->addressable_end;
  const size_t slots      = clv_man->slottable_size;

  // firstly, allocate the mrc_data struct
  mrc_data_t * mrc = clv_man->repl_strat_data = 
                                      (mrc_data_t *)malloc(sizeof(mrc_data_t));
  if (!mrc)
  {
    MRC_dealloc_cb(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for mrc_data.");
    return PLL_FAILURE;
  }

  // then, its parts
  mrc->cost_of_clvid = (unsigned int*)calloc(addr_size, sizeof(unsigned int));
  if (!mrc->cost_of_clvid)
  {
    MRC_dealloc_cb(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for cost_of_clvid.");
    return PLL_FAILURE;
  }

  mrc->cost_of_slot = (unsigned int*)alloc_and_set(slots,
                                                   sizeof(unsigned int),
                                                   -1);
  if (!mrc->cost_of_slot)
  {
    MRC_dealloc_cb(clv_man);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for cost_of_slot.");
    return PLL_FAILURE;
  }
  return PLL_SUCCESS;
}

/**
 * Initializes the replacement strategy for the default MRC replacement.
 *
 * Basically just allocs / fills the array holding the clv indices sorted by
 * their recomputation cost. Does so according to the structure of the tree:
 * each node's recomputation cost is essentially their subtree size.
 * 
 * @param  clv_man        the pll_clv_manager struct
 * @param  tree           (utree) tree structure
 * @param  subtree_sizes  array mapping from node_index to subtree size
 *                        (see pll_utree_get_subtree_sizes)
 *
 * @return                PLL_{SUCCESS|FAILURE}
 */
PLL_EXPORT int pll_clv_manager_MRC_strategy_init(pll_clv_manager_t * clv_man,
                                        pll_utree_t * const tree,
                                        unsigned int const * const subtree_sizes)
{
  assert(subtree_sizes);

  // initialize tree-independent data structures
  if (!MRC_strategy_init(clv_man))
    return PLL_FAILURE;

  mrc_data_t * mrc = clv_man->repl_strat_data;

  // next, we do the initial cost computation
  // TODOCLV this should be also done in a separate MRC function, that recomputes
  // the sizes and sets them in the MRC size arrays

  const size_t nodes_count = tree->tip_count + tree->inner_count;

  // go through all the nodes of the tree and fill the cost array
  for (size_t i = 0; i < nodes_count; ++i)
  {
    pll_unode_t* node = tree->nodes[ i ];
    mrc->cost_of_clvid[ node->clv_index ] = subtree_sizes[ node->node_index ];
  }

  // finally, set up the cost per slot array
  // as this is used here as part of the init, probably they will all be unused 
  const unsigned int * const clvid_of_slot = clv_man->clvid_of_slot;
  for (size_t slot_id = 0; slot_id < clv_man->slottable_size; ++slot_id)
  {
    unsigned int clv_id = clvid_of_slot[slot_id];
    if (clv_id != PLL_CLV_SLOT_UNUSED)
    {
      mrc->cost_of_slot[ slot_id ] = mrc->cost_of_clvid[ clv_id ];
    }
  }

  return PLL_SUCCESS;
}

/**
 * Initializes the replacement strategy for the default MRC replacement, rtree version.
 *
 * Basically just allocs / fills the array holding the clv indices sorted by
 * their recomputation cost. Does so according to the structure of the tree:
 * each node's recomputation cost is essentially their subtree size.
 * 
 * @param  clv_man        the pll_clv_manager struct
 * @param  tree           (rtree) tree structure
 * @param  subtree_sizes  array mapping from node_index to subtree size
 *                        (see pll_utree_get_subtree_sizes)
 *
 * @return                PLL_{SUCCESS|FAILURE}
 */
PLL_EXPORT int pll_clv_manager_MRC_strategy_rtree_init(
                                        pll_clv_manager_t * clv_man,
                                        pll_rtree_t * const tree,
                                        unsigned int const * const subtree_sizes)
{
  assert(subtree_sizes);

  // initialize tree-independent data structures
  if (!MRC_strategy_init(clv_man))
    return PLL_FAILURE;

  mrc_data_t * mrc = clv_man->repl_strat_data;

  // next, we do the initial cost computation
  // TODOCLV this should be also done in a separate MRC function, that recomputes
  // the sizes and sets them in the MRC size arrays

  const size_t nodes_count = tree->tip_count + tree->inner_count;

  // go through all the nodes of the tree and fill the cost array
  for (size_t i = 0; i < nodes_count; ++i)
  {
    pll_rnode_t* node = tree->nodes[ i ];
    mrc->cost_of_clvid[ node->clv_index ] = subtree_sizes[ node->node_index ];
  }

  // finally, set up the cost per slot array
  // as this is used here as part of the init, probably they will all be unused 
  const unsigned int * const clvid_of_slot = clv_man->clvid_of_slot;
  for (size_t slot_id = 0; slot_id < clv_man->slottable_size; ++slot_id)
  {
    unsigned int clv_id = clvid_of_slot[slot_id];
    if (clv_id != PLL_CLV_SLOT_UNUSED)
    {
      mrc->cost_of_slot[ slot_id ] = mrc->cost_of_clvid[ clv_id ];
    }
  }

  return PLL_SUCCESS;
}

