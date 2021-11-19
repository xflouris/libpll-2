/*
    Copyright (C) 2016-2020 Tomas Flouri, Alexey Kozlov

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

/* simulate exactly the non-reentrant glibc srandom() function */
#define RAND_STATE_SIZE 128

typedef struct
{
  int clv_valid;
} node_info_t;

typedef struct
{
  pll_parsimony_t ** pars_list;
  unsigned int pars_count;
  pll_unode_t ** travbuffer;
  pll_pars_buildop_t * parsops;
  unsigned int ops_count;
  unsigned int traversal_size;
} pars_info_t;

static char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)malloc(len+1);
  if (!p)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Memory allocation failed");
    return NULL;
  }
  return strcpy(p,s);
}

/* Fisher-Yates shuffle */
static unsigned int * create_shuffled(unsigned int n, unsigned int seed)
{
  unsigned int i,j;
  char * statebuf;
  struct pll_random_data * buf;
  
  unsigned int * x = (unsigned int *)malloc(n*sizeof(unsigned int));
  if (!x)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  for (i=0; i<n; ++i)
    x[i] = i;

  /* if seed == 0 then do not shuffle! */
  if (!seed)
    return x;

  /* init re-entrant randomizer */
  buf = (struct pll_random_data *)calloc(1, sizeof(struct pll_random_data));
  statebuf = (char *)calloc(RAND_STATE_SIZE,sizeof(char));

  pll_initstate_r(seed,statebuf,RAND_STATE_SIZE,buf);
  pll_srandom_r(seed,buf);

  /* perform Fisher-Yates shuffle */
  if (n > 1)
  {
    i = n - 1;
    while (1)
    {
      int rint;
      pll_random_r(buf,&rint);
      double r = ((double)rint / RAND_MAX);
      j = (unsigned int)(r * (i+1));

      PLL_SWAP(x[i],x[j]);

      if (i == 0) break;
      --i;
    }
  }

  /* dealloc and return shuffled array */
  free(statebuf);
  free(buf);
  return x;
}

static void dealloc_data_onenode(pll_unode_t * node)
{
  if (node->data)
  {
    free(node->data);
    node->data = NULL;
  }
}

static void dealloc_data(pll_unode_t * node)
{
  if (node->next)
  {
    dealloc_data_onenode(node);
    dealloc_data_onenode(node->next);
    dealloc_data_onenode(node->next->next);
  }
}

static pars_info_t * create_pars_info(pll_parsimony_t ** pars_list,
                                      unsigned int pars_count)
{
  if (!pars_list)
    return NULL;

  unsigned int tip_count = pars_list[0]->tips;


  pars_info_t * pars_info = (pars_info_t *) calloc(1, sizeof(pars_info_t));

  if (!pars_info)
    return NULL;

  pars_info->pars_list = pars_list;
  pars_info->pars_count = pars_count;

  pars_info->travbuffer = (pll_unode_t **)malloc((2*tip_count-2) * sizeof(pll_unode_t *));

  pars_info->parsops = (pll_pars_buildop_t *)malloc((tip_count-2)*
                                         sizeof(pll_pars_buildop_t));

  return pars_info;
}

static void destroy_pars_info(pars_info_t * pars_info)
{
  if (pars_info)
  {
    free(pars_info->travbuffer);
    free(pars_info->parsops);
    free(pars_info);
  }
}

static void invalidate_node(pll_unode_t * node)
{
  node_info_t * info;

  assert(node->data);

  info = (node_info_t *)(node->data);
  info->clv_valid = 0;
  info = (node_info_t *)(node->next->data);
  info->clv_valid = 0;
  info = (node_info_t *)(node->next->next->data);
  info->clv_valid = 0;
}


/* a callback function for performing a partial traversal */
static int cb_partial_traversal(pll_unode_t * node)
{
  node_info_t * node_info;

  /* if we don't want tips in the traversal we must return 0 here. For now,
     allow tips */
  if (!node->next) return 1;

  /* get the data element from the node -> all inner nodes must have
   * pre-allocated data */
  node_info = (node_info_t *)(node->data);

  assert(node_info);

  /* if the CLV is valid, we instruct the traversal routine not to
     traverse the subtree rooted in this node/direction by returning 0 */
  if (node_info->clv_valid) return 0;

  /* otherwise, mark CLV as valid since it will be updated */
  node_info->clv_valid = 1;

  return 1;
}

static int cb_validate(pll_unode_t * node)
{
  if (node->data)
  {
    node_info_t * node_info = (node_info_t *)(node->data);
    node_info->clv_valid = 1;
  }

  return 1;
}

static int cb_invalidate(pll_unode_t * node)
{
  if (!node->next)
    node = node->back;

  invalidate_node(node);

  return 1;
}


static int cb_full(pll_unode_t * node)
{
  return 1;
}

static int cb_full_subtree(pll_unode_t * node)
{
  /* ignore subtrees which "dead-end" subtrees with empty back pointers */
  return !node->next || (node->next->back && node->next->next->back);
}


static pll_unode_t * utree_inner_create(unsigned int i, unsigned int tip_count)
{
  pll_unode_t * node = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node)
    return NULL;
  
  node->next = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node->next)
  {
    free(node);
    return NULL;
  }
  node->next->next = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  if (!node->next->next)
  {
    free(node->next);
    free(node);
    return NULL;
  }

  /* allocate data element */
  node->data             = (node_info_t *)calloc(1,sizeof(node_info_t));
  node->next->data       = (node_info_t *)calloc(1,sizeof(node_info_t));
  node->next->next->data = (node_info_t *)calloc(1,sizeof(node_info_t));

  if (!node->data || !node->next->data || !node->next->next->data)
  {
    free(node->next->next->data);
    free(node->next->data);
    free(node->data);
    free(node->next->next);
    free(node->next);
    free(node);
    return NULL;
  }

  node->next->next->next = node;

  unsigned int clv_id = tip_count + i;
  node->clv_index = clv_id;
  node->next->clv_index = clv_id;
  node->next->next->clv_index = clv_id;

  unsigned int node_id = tip_count + i*3;
  node->node_index = node_id;
  node->next->node_index = node_id + 1;
  node->next->next->node_index = node_id + 2;

  return node;
}

static pll_unode_t * utree_tip_create(unsigned int i)
{
  pll_unode_t * node = (pll_unode_t *)calloc(1,sizeof(pll_unode_t));
  node->next = NULL;
  node->clv_index = i;
  node->node_index = i;

  return node;
}

static void utree_link(pll_unode_t * a, pll_unode_t * b)
{
  /*

    *               *               *                * 
     \             /                 \              /
      *---*   *---*        -->        *---*-----*--*
     /    a   b    \                 /    a     b   \
    *               *               *                *

  */

  a->back = b;
  b->back = a;
  b->pmatrix_index = a->pmatrix_index;
}

static void utree_edgesplit(pll_unode_t * a, pll_unode_t * b, pll_unode_t * c)
{
  /*
                *                                      *
                |                                      |
                *                                      *
               / \                                    / \
            b *   * c                              b *   * c
                                                    /     \
    *                      *      -->      *       /       \      * 
     \                    /                 \     /         \    /
      *---*----------*---*                   *---*           *--*
     /    a          d    \                 /    a           d   \
    *                      *               *                      *

  */

  /* link d<->c */
  utree_link(c, a->back);

  /* link a<->b */
  utree_link(a,b);
}

static pll_unode_t * utree_prune(pll_unode_t * p)
{
  pll_unode_t * a = p->next->back;
  pll_unode_t * b = p->next->next->back;
  utree_link(a, b);
  p->next->back = p->next->next->back = NULL;
  return a;
}

static int utree_is_tip(pll_unode_t * node)
{
  return (node->next == NULL);
}

static int utree_collect_edges(pll_unode_t * root,
                               pll_unode_t ** edge_list,
                               unsigned int * edge_count)
{
  unsigned int i;

  if (!pll_utree_traverse(root,
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_full,
                          edge_list,
                          edge_count))
    return PLL_FAILURE;

  for (i = 0; i < *edge_count; ++i)
  {
    if (!edge_list[i]->next)
      edge_list[i] = edge_list[i]->back;
  }

  // root edge is traversed twice -> correct for this
  (*edge_count)--;

  return PLL_SUCCESS;
}

static int utree_update_pars_vectors(pll_unode_t * root,
                                     pars_info_t * pars_info,
                                     int trav_type)
{
  unsigned int i;

  if (trav_type == PLL_TREE_TRAVERSE_NONE)
  {
    /* update a single CLV at the root */
    pars_info->travbuffer[0] = root;
    pars_info->traversal_size = 1;
  }
  else
  {
    /* update all CLVs (full traversal) or invalid CLVs only (partial) */
    if (!pll_utree_traverse(root,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            trav_type == PLL_TREE_TRAVERSE_PARTIAL ?
                                cb_partial_traversal : cb_full_subtree,
                            pars_info->travbuffer,
                            &pars_info->traversal_size))
      return PLL_FAILURE;
  }

  /* create parsimony operations */
  pll_utree_create_pars_buildops(pars_info->travbuffer,
                                 pars_info->traversal_size,
                                 pars_info->parsops,
                                 &pars_info->ops_count);

  for (i = 0; i < pars_info->pars_count; ++i)
  {
    /* update parsimony vectors */
    pll_fastparsimony_update_vectors(pars_info->pars_list[i],
                                     pars_info->parsops,
                                     pars_info->ops_count);
  }

  return PLL_SUCCESS;
}

static unsigned int utree_pars_edge_score(pll_unode_t * edge,
                                          pars_info_t * pars_info)
{
  unsigned int i;
  unsigned int cost = 0;

  for (i = 0; i < pars_info->pars_count; ++i)
  {
    /* get parsimony score */
    cost += pll_fastparsimony_edge_score(pars_info->pars_list[i],
                                         edge->node_index,
                                         edge->back->node_index);
  }

  return cost;
}


static unsigned int utree_insert_best(pars_info_t * pars_info,
                                      pll_unode_t ** edge_list,
                                      unsigned int edge_count,
                                      pll_unode_t * inner_node,
                                      const unsigned int * constraint,
                                      pll_unode_t * prune_edge)
{
  unsigned int i;
  unsigned int min_cost;
  unsigned int best_index;
  unsigned int cost;
  size_t total_ops = 0;

  /* subtree must be pruned  / not inserted yet */
  assert(!inner_node->next->back && !inner_node->next->next->back);

  /* set min cost to maximum possible value */
  min_cost = ~0u;

  /* find first empty slot in edge_list */
  pll_unode_t ** empty_slot = edge_list + edge_count;

  /* fill *all* CLV vectors in all directions, ie 3 CLVs per inner nodes ->
   * this way, we can do avoid unnecessary CLV recomputation when
   * evaluating insertion branches in the loop below */
  for (i = 0; i < edge_count; ++i)
  {
    pll_unode_t * root = edge_list[i]->next ? edge_list[i] : edge_list[i]->back;

    /* traverse from every OUTER branch */
    if (root->back->next)
      continue;

    /* make a partial traversal */
    utree_update_pars_vectors(root, pars_info, PLL_TREE_TRAVERSE_PARTIAL);

    total_ops += pars_info->ops_count;
  }

  /* if we insert a subtree, recompute all CLVs in this subtree
   * in the direction of re-insertion point */
  if (!utree_is_tip(inner_node->back))
    utree_update_pars_vectors(inner_node->back, pars_info, PLL_TREE_TRAVERSE_FULL);

  cost = pll_fastparsimony_edge_score(pars_info->pars_list[0],
                                      edge_list[0]->node_index,
                                      edge_list[0]->back->node_index);
//  printf("start cost: %u\n", cost);

  best_index = edge_count + 1;
  for (i = 0; i < edge_count; ++i)
  {
    pll_unode_t * regraft_edge = edge_list[i];
    pll_unode_t * regraft_back = regraft_edge->back;

    /* check constraint */
    if (constraint)
    {
      unsigned int s = constraint[inner_node->clv_index];
      unsigned int r1 = constraint[regraft_edge->clv_index];
      unsigned int r2 = constraint[regraft_back->clv_index];

      assert(s);
      if (s && s != r1 && s != r2)
        continue;
    }

    /* split the regraft edge and insert subtree rooted at inner_node */
    utree_edgesplit(regraft_edge, inner_node->next, inner_node->next->next);

    /* we only need to recompute one CLV vector at the inner node */
    utree_update_pars_vectors(inner_node, pars_info, PLL_TREE_TRAVERSE_NONE);

    total_ops += pars_info->ops_count;

    /* compute the cost for all parsimony partitions */
    cost = utree_pars_edge_score(inner_node, pars_info);

    /* if current cost is smaller than minimum cost save branch index */
    if (cost < min_cost)
    {
      min_cost = cost;
      best_index = i;
    }

    /* restore tree to its state before placing the tip (and inner) node */
    utree_link(regraft_edge, regraft_back);
    inner_node->next->back = NULL;
    inner_node->next->next->back = NULL;
  }

  /* perform the placement yielding the lowest cost */
  if (best_index < edge_count)
  {
//    printf("best cost: %u\n", min_cost);
    utree_edgesplit(edge_list[best_index], inner_node->next, inner_node->next->next);
  }
  else
  {
    // no valid placements found (only with constraint!) -> regraft to original edge
    assert(constraint && prune_edge);

    utree_edgesplit(prune_edge, inner_node->next, inner_node->next->next);

    utree_update_pars_vectors(inner_node, pars_info, PLL_TREE_TRAVERSE_NONE);

    total_ops += pars_info->ops_count;

    min_cost = utree_pars_edge_score(inner_node, pars_info);

//    printf("rollback cost: %u\n", min_cost);
  }

  /* add the two new edges to the end of the list */
  if (!prune_edge)
  {
    /* inner_node->next is linked to edge_list[best_index], so it is already in the list */
    empty_slot[0] = inner_node;
    empty_slot[1] = inner_node->next->next;
  }

  /* invalidate all CLVs */
  if (!pll_utree_traverse(edge_list[0],
                          PLL_TREE_TRAVERSE_POSTORDER,
                          cb_invalidate,
                          pars_info->travbuffer,
                          &pars_info->traversal_size))
    assert(0);

  /* re-validate CLVs that remain correct after new tip insertion -> not for SPR! */
  if (!prune_edge)
  {
    if (!pll_utree_traverse(inner_node,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_validate,
                            pars_info->travbuffer,
                            &pars_info->traversal_size))
      assert(0);
  }

  /* reset direction for the newly placed inner node */
  invalidate_node(inner_node);

  if (!utree_is_tip(inner_node->back))
    invalidate_node(inner_node->back);

  return min_cost;
}

PLL_EXPORT int pll_fastparsimony_stepwise_spr_round(pll_utree_t * tree,
                                                   pll_parsimony_t ** pars_list,
                                                   unsigned int pars_count,
                                                   const unsigned int * tip_msa_idmap,
                                                   unsigned int seed,
                                                   const int * clv_index_map,
                                                   unsigned int * cost)
{
  unsigned int i;
  unsigned int old_edge_count = tree->edge_count;
  unsigned int tip_count = tree->tip_count;
  unsigned int inner_count = tree->inner_count;
  unsigned int node_count = tip_count + inner_count;
  unsigned int edge_count = old_edge_count;
  unsigned int subtree_count = inner_count * 3;
  unsigned int new_tip_count = pars_list[0]->tips;
  unsigned int ext_tip_count = new_tip_count - tip_count;

  pars_info_t * pars_info = create_pars_info(pars_list, pars_count);

  pll_unode_t ** all_nodes = (pll_unode_t **) calloc(subtree_count,
                                                     sizeof(pll_unode_t *));

  pll_unode_t ** edge_list = (pll_unode_t **) calloc(old_edge_count,
                                                     sizeof(pll_unode_t *));

  unsigned int * constraint = (unsigned int *) calloc(node_count,
                                                      sizeof(unsigned int));

  unsigned int * orig_idmap = (unsigned int *) calloc(new_tip_count,
                                                      sizeof(unsigned int));
  for (i = 0; i < node_count; ++i)
  {
    unsigned int clv_id = tree->nodes[i]->clv_index;
    constraint[clv_id] = tree->nodes[i]->next ? clv_index_map[clv_id]+1 : 0;
  }

  /* special treatment for incomplete constraint trees:
   * tip indexing in the tree is different from MSA ordering */
  if (tip_msa_idmap)
  {
    /* remap node_index to be consistent with numbering in pll_parsimony_t ! */
    for (i = 0; i < tip_count; ++i)
    {
      unsigned int old_idx = tree->nodes[i]->node_index;
      unsigned int new_idx = tip_msa_idmap[old_idx];
      tree->nodes[i]->node_index = new_idx;
      orig_idmap[new_idx] = old_idx;
    }

    /* update node_index of inner nodes to correct for additional tips */
    for (i = tip_count; i < node_count; ++i)
    {
      pll_unode_t * node = tree->nodes[i];
      assert(node->next);
      node->node_index += ext_tip_count;
      node->next->node_index += ext_tip_count;
      node->next->next->node_index += ext_tip_count;
    }
  }

  unsigned int * order = create_shuffled(subtree_count, seed);
  if (!order)
    return PLL_FAILURE;

  /* collect all nodes */
  for (i = 0; i < inner_count; ++i)
  {
    pll_unode_t * node = tree->nodes[tip_count + i];
    assert(node->next);
    all_nodes[3*i] = node;
    all_nodes[3*i + 1] = node->next;
    all_nodes[3*i + 2] = node->next->next;
  }

  for (i = 0; i < subtree_count; ++i)
  {
    all_nodes[i]->data = (node_info_t *) calloc(1, sizeof(node_info_t));
  }

  /* prune and regraft subtrees in random order */
  for (i = 0; i < subtree_count; ++i)
  {
    pll_unode_t * new_inner = all_nodes[order[i]];
    pll_unode_t * new_root = NULL;
    pll_unode_t * prune_edge = NULL;

    assert(new_inner->next);

    /* if remaining pruned tree would only contain 2 taxa, skip this node */
    if (utree_is_tip(new_inner->next->back) &&
        utree_is_tip(new_inner->next->next->back))
      continue;

    /* prune a subtree */
    prune_edge = utree_prune(new_inner);

    new_root = prune_edge->next ? prune_edge : prune_edge->back;

    // collect remaining edges
    utree_collect_edges(new_root, edge_list, &edge_count);

    *cost = utree_insert_best(pars_info,
                              edge_list,
                              edge_count,
                              new_inner,
                              constraint,
                              prune_edge);
  }

  /* restore original node_index */
  if (tip_msa_idmap)
  {
    for (i = 0; i < tip_count; ++i)
    {
      unsigned int new_idx = tree->nodes[i]->node_index;
      unsigned int old_idx = orig_idmap[new_idx];
      tree->nodes[i]->node_index = old_idx;
    }

    /* update node_index of inner nodes */
    for (i = tip_count; i < node_count; ++i)
    {
      pll_unode_t * node = tree->nodes[i];
      assert(node->next);
      node->node_index -= ext_tip_count;
      node->next->node_index -= ext_tip_count;
      node->next->next->node_index -= ext_tip_count;
    }
  }

  destroy_pars_info(pars_info);

  /* delete data elements */
  for (i = 0; i < node_count; ++i)
    dealloc_data(tree->nodes[i]);

  free(order);
  free(edge_list);
  free(all_nodes);
  free(constraint);
  free(orig_idmap);

  return PLL_SUCCESS;
}

PLL_EXPORT int pll_fastparsimony_stepwise_extend(pll_utree_t * tree,
                                                 pll_parsimony_t ** pars_list,
                                                 unsigned int pars_count,
                                                 char * const * labels,
                                                 const unsigned int * tip_msa_idmap,
                                                 unsigned int seed,
                                                 unsigned int * cost)
{
  unsigned int i,j;
  unsigned int new_tip_count = pars_list[0]->tips;
  unsigned int new_inner_count = new_tip_count - 2;
  unsigned int new_node_count = new_tip_count + new_inner_count;
  unsigned int new_edge_count = 2 * new_tip_count - 3;
  unsigned int old_tip_count = tree->tip_count;
  unsigned int old_inner_count = tree->inner_count;
  unsigned int old_node_count = old_tip_count + old_inner_count;
  unsigned int old_edge_count = tree->edge_count;
  unsigned int ext_tip_count = new_tip_count - old_tip_count;
  unsigned int edge_count;

  pars_info_t * pars_info = create_pars_info(pars_list, pars_count);

  pll_unode_t ** old_nodes = tree->nodes;
  pll_unode_t ** new_nodes = (pll_unode_t **) calloc(new_node_count,
                                                     sizeof(pll_unode_t *));
  pll_unode_t ** edge_list = (pll_unode_t **) calloc(new_edge_count+1,
                                                     sizeof(pll_unode_t *));
  unsigned int * order = create_shuffled(ext_tip_count, seed);

  if (!new_nodes || !edge_list || !order)
  {
    free(new_nodes);
    free(edge_list);
    free(order);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Cannot allocate memory for nodes!");

    return PLL_FAILURE;
  }

  /* 1:1 mapping for old tips */
  for (i = 0; i < old_tip_count; ++i)
    new_nodes[i] = old_nodes[i];

  /* copy old inner nodes and adjust CLVs */
  for (i = old_tip_count; i < old_node_count; ++i)
  {
    unsigned int new_idx = i + ext_tip_count;
    new_nodes[new_idx] = old_nodes[i];
    pll_unode_t * snode = new_nodes[new_idx];
    assert(snode->next);
    do
    {
      snode->clv_index += ext_tip_count;
      snode->node_index += ext_tip_count;
      snode->data = (node_info_t *) calloc(1, sizeof(node_info_t));

      snode = snode->next;
    }
    while (snode != new_nodes[new_idx]);
  }

  /* create new tip and inner nodes */
  for (i = 0; i < ext_tip_count; ++i)
  {
    unsigned int tip_pos = old_tip_count + i;
    unsigned int inner_pos = new_tip_count + old_inner_count + i;
    unsigned int index = order[i] + old_tip_count;
    pll_unode_t * tip_node = utree_tip_create(index);
    pll_unode_t * inner_node = utree_inner_create(old_inner_count + i, new_tip_count);

    if (!tip_node || !inner_node)
    {
      free(new_nodes);
      free(edge_list);
      free(order);

      for (j = 0; j < i; ++j)
      {
        free(new_nodes[old_tip_count + j]);
        free(new_nodes[new_tip_count + old_inner_count + j]);
      }

      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Cannot allocate memory for nodes!");

      return PLL_FAILURE;
    }

    tip_node->label = xstrdup(labels[index - old_tip_count]);
    new_nodes[tip_pos] = tip_node;
    new_nodes[inner_pos] = inner_node;

    /* connect tip node with respective inner node */
    utree_link(inner_node, tip_node);
  }

  if (tip_msa_idmap)
  {
    /* remap node_index to be consistent with numbering in pll_parsimony_t !*/
    for (i = 0; i < new_tip_count; ++i)
    {
      unsigned int old_idx = new_nodes[i]->node_index;
      new_nodes[i]->node_index = tip_msa_idmap[old_idx];
    }
  }

  /* collect all edges */
  utree_collect_edges(tree->vroot, edge_list, &edge_count);

  assert(edge_count == old_edge_count);

  pll_unode_t ** new_inner_nodes = new_nodes + new_tip_count + old_inner_count;
  for (i = 0; i < ext_tip_count; ++i)
  {
//    printf("%d -- adding %u %s\n", i, new_tip_nodes[i]->clv_index, new_tip_nodes[i]->label);
    pll_unode_t * new_tip_node = new_inner_nodes[i]->back;
    assert(!new_tip_node->next);

    *cost = utree_insert_best(pars_info,
                          edge_list,
                          edge_count,
                          new_inner_nodes[i],
                          NULL,
                          NULL);

    /* after adding a leaf, we have two new edges */
    edge_count += 2;
  }

  assert(edge_count == new_edge_count);

  tree->nodes = new_nodes;
  tree->tip_count = new_tip_count;
  tree->inner_count = new_inner_count;
  tree->edge_count = edge_count;
  tree->vroot = tree->vroot->next ? tree->vroot : tree->vroot->back;

  destroy_pars_info(pars_info);

  /* delete data elements */
  for (i = 0; i < new_node_count; ++i)
    dealloc_data(new_nodes[i]);

  free(edge_list);
  free(old_nodes);
  free(order);

  return PLL_SUCCESS;
}

PLL_EXPORT pll_utree_t * pll_fastparsimony_stepwise(pll_parsimony_t ** list,
                                                    char * const * labels,
                                                    unsigned int * cost,
                                                    unsigned int count,
                                                    unsigned int seed)
{
  unsigned int i,j;

  unsigned int tips_count = list[0]->tips;
  unsigned int inner_nodes = list[0]->inner_nodes;

  if (tips_count < 3)
  {
    pll_errno = PLL_ERROR_STEPWISE_TIPS;
    snprintf(pll_errmsg, 200,
             "Stepwise parsimony requires at least three tips.");
    return NULL;
  }

  //if (tips_count != inner_nodes + 2)
  if (inner_nodes < tips_count-2)
  {
    pll_errno = PLL_ERROR_STEPWISE_UNSUPPORTED;
    snprintf(pll_errmsg, 200,
             "Stepwise parsimony currently supports only unrooted trees.");
    return NULL;
  }

  *cost = ~0u;

  pll_unode_t * root;

  /* check that all parsimony structures have the same number of tips and
     inner nodes */

  for (i = 1; i < count; ++i)
  {
    if ((list[i]->tips != tips_count) ||
        (list[i]->inner_nodes != inner_nodes))
    {
      pll_errno = PLL_ERROR_STEPWISE_STRUCT;
      snprintf(pll_errmsg, 200,
               "Parsimony structures tips/inner nodes not equal.");
      return NULL;
    }
  }
    

  /* 1. Make all allocations at the beginning and check everything was
        allocated, otherwise return an error */

  pars_info_t * pars_info = create_pars_info(list, count);

  root = utree_inner_create(tips_count-3, tips_count);

  /* create tip node list with a terminating NULL element */
  pll_unode_t ** tip_node_list = (pll_unode_t **)calloc(tips_count+1,
                                                        sizeof(pll_unode_t *));

  /* create inner node list for (tips_count - 3) inner nodes (root was already
     created, and leave the last slot NULL for termination */
  pll_unode_t ** inner_node_list = (pll_unode_t **)calloc(tips_count - 2,
                                                          sizeof(pll_unode_t *));

  if (!inner_node_list || !pars_info || !tip_node_list || !root)
  {
    pll_utree_graph_destroy(root,NULL);
    destroy_pars_info(pars_info);
    free(inner_node_list);
    free(tip_node_list);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  /* allocate all inner nodes */
  for (i=0; i<tips_count-3; ++i)
  {
    inner_node_list[i] = utree_inner_create(i, tips_count);
    if (!inner_node_list[i])
    {
      pll_utree_graph_destroy(root,NULL);
      destroy_pars_info(pars_info);
      free(tip_node_list);
      for (j = 0; j < i; ++j)
        pll_utree_graph_destroy(inner_node_list[j],NULL);
      free(inner_node_list);

      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
  }

  /* shuffle the order of iterating tip sequences */
  unsigned int * order = create_shuffled(tips_count, seed);
  if (!order) return NULL;

  /* allocate all tips */
  for (i=0; i<tips_count; ++i)
  {
    unsigned int index = order[i];
    tip_node_list[i] = utree_tip_create(index);
    if (tip_node_list[i])
    {
      tip_node_list[i]->label = xstrdup(labels[index]);
      if (i > 2)
        utree_link(inner_node_list[i-3], tip_node_list[i]);
    }

    if (!tip_node_list[i] || !tip_node_list[i]->label)
    {
      free(tip_node_list[i]); 

      pll_utree_graph_destroy(root,NULL);
      destroy_pars_info(pars_info);
      free(inner_node_list);
      for (j = 0; j < i; ++j)
        pll_utree_graph_destroy(tip_node_list[j],NULL);
      free(tip_node_list);

      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;

    }
  }
  free(order);

  /* 2. Create the following topology with three leaves
          
            *
           /
      *---*
           \
            *
  
  */

  /* place first three tips */
  utree_link(root, tip_node_list[0]);
  utree_link(root->next, tip_node_list[1]);
  utree_link(root->next->next, tip_node_list[2]);

  /* available placements */
  pll_unode_t ** edge_list = (pll_unode_t **)calloc(2*tips_count-3,
                                                    sizeof(pll_unode_t *));
  edge_list[0] = root;
  edge_list[1] = root->next;
  edge_list[2] = root->next->next;

  /* 3. The stepwise parsimony. Current topology is the tree with three leaves,
        and repeat the following steps for each remaining tip u:
         (i) compute the parsimony score of all possible topologies by placing
             u at every possible edge of the current tree topology.
        (ii) set current toplogy as the tree with the smallest parsimony score
  */
  if (tips_count > 3)
  {
    unsigned int edge_count = 3;
    
    for (i = 3; i < tips_count; ++i)
    {
      /* printf("%d -- adding %s\n", i, tip_node_list[i]->label); */
      *cost = utree_insert_best(pars_info,
                                edge_list,
                                edge_count,
                                inner_node_list[i-3],
                                NULL,
                                NULL);

      /* after adding a leaf, we have two new edges */
      edge_count += 2;
    }
  }
  else
  {
    *cost = 0;
    for (i = 0; i < count; ++i)
      *cost += list[i]->const_cost;
  }

  /* delete data elements */
  for (i = 0; i < tips_count-3; ++i)
    dealloc_data(inner_node_list[i]);
  dealloc_data(root);

  destroy_pars_info(pars_info);

  /* deallocate auxiliary arrays */
  free(inner_node_list);
  free(tip_node_list);
  free(edge_list);

  /* wrap tree */
  pll_utree_t * tree = pll_utree_wraptree(root,tips_count);

  return tree;
}
