/*
    Copyright (C) 2015-2018 Tomas Flouri, Alexey Kozlov

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
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static int indent_space = 4;

static void print_node_info(const pll_unode_t * node, int options)
{
  if (options & PLL_UTREE_SHOW_LABEL)
    printf (" %s", node->label);
  if (options & PLL_UTREE_SHOW_BRANCH_LENGTH)
    printf (" %f", node->length);
  if (options & PLL_UTREE_SHOW_CLV_INDEX)
    printf (" %d", node->clv_index);
  if (options & PLL_UTREE_SHOW_SCALER_INDEX)
    printf (" %d", node->scaler_index);
  if (options & PLL_UTREE_SHOW_PMATRIX_INDEX)
    printf (" %d", node->pmatrix_index);
  if (options & PLL_UTREE_SHOW_DATA)
    printf (" %p", node->data);
  printf("\n");
}

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

static void print_tree_recurse(pll_unode_t * node,
                               int indent_level,
                               int * active_node_order,
                               int options)
{
  int i,j;

  if (!node) return;

  for (i = 0; i < indent_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indent_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indent_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indent_space-1; ++j)
    printf ("-");
  if (node->next) printf("+");

  print_node_info(node, options);

  if (active_node_order[indent_level-1] == 2)
    active_node_order[indent_level-1] = 0;

  if (node->next)
  {
    pll_unode_t * snode = node->next;
    do {
      active_node_order[indent_level] = snode->next == node ? 2 : 1;
      print_tree_recurse(snode->back,
                         indent_level+1,
                         active_node_order,
                         options);
      snode = snode->next;
    }
    while (snode != node);
  }

}

static unsigned int tree_indent_level(const pll_unode_t * node, unsigned int indent)
{
  if (!node->next)
    return indent+1;

  unsigned int ind = 0;
  pll_unode_t * snode = node->next;
  do
  {
    unsigned int sind = tree_indent_level(snode->back, indent+1);
    ind = PLL_MAX(ind, sind);
    snode = snode->next;
  }
  while (snode && snode != node);

  return ind;
}

PLL_EXPORT void pll_utree_show_ascii(const pll_unode_t * root, int options)
{
  unsigned int a, b;

  if (!root->next) root=root->back;

  a = tree_indent_level(root->back,1);
  b = tree_indent_level(root,0);
  unsigned int max_indent_level = (a > b ? a : b);

  int * active_node_order = (int *)malloc((max_indent_level+1) * sizeof(int));
  if (!active_node_order)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return;
  }
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  const pll_unode_t * node = root;
  do {
    active_node_order[0] = node->next == root ? 2 : 1;
    print_tree_recurse(node->back, 1, active_node_order, options);
    node = node->next;
  }
  while (node != root);
  free(active_node_order);
}

static char * newick_utree_recurse(const pll_unode_t * root,
                                   char * (*cb_serialize)(const pll_unode_t *))
{
  char * newick;
  int size_alloced;
  assert(root != NULL);
  if (!root->next)
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = newick_utree_recurse(root->next->back,cb_serialize);
    if (subtree1 == NULL)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
    char * subtree2 = newick_utree_recurse(root->next->next->back,cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s:%f",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

char * utree_export_newick(const pll_unode_t * root,
                           int export_rooted,
                           double root_brlen,
                           char * (*cb_serialize)(const pll_unode_t *))
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (!root->next) root=root->back;

  char * subtree1 = newick_utree_recurse(root->back,cb_serialize);
  if (subtree1 == NULL)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }
  char * subtree2 = newick_utree_recurse(root->next->back,cb_serialize);
  if (subtree2 == NULL)
  {
    free(subtree1);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }
  char * subtree3 = newick_utree_recurse(root->next->next->back,cb_serialize);
  if (subtree3 == NULL)
  {
    free(subtree1);
    free(subtree2);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  if (export_rooted)
  {
    assert(!cb_serialize);
    size_alloced = asprintf(&newick,
                            "(%s,(%s,%s):%f)%s:0.0;",
                            subtree1,
                            subtree2,
                            subtree3,
                            root_brlen,
                            root->label ? root->label : "");
  }
  else
  {
    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = asprintf(&newick,
                              "(%s,%s,%s)%s",
                              subtree1,
                              subtree2,
                              subtree3,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = asprintf(&newick,
                              "(%s,%s,%s)%s:0.0;",
                              subtree1,
                              subtree2,
                              subtree3,
                              root->label ? root->label : "");
    }
  }
  free(subtree1);
  free(subtree2);
  free(subtree3);
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return (newick);
}

PLL_EXPORT char * pll_utree_export_newick(const pll_unode_t * root,
                                   char * (*cb_serialize)(const pll_unode_t *))
{
  return utree_export_newick(root, 0, 0, cb_serialize);
}

PLL_EXPORT char * pll_utree_export_newick_rooted(const pll_unode_t * root,
                                                 double root_brlen)
{
  return utree_export_newick(root, 1, root_brlen, NULL);
}

PLL_EXPORT void pll_utree_create_operations(pll_unode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count)
{
  const pll_unode_t * node;
  unsigned int i;

  *ops_count = 0;
  *matrix_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    /* if the current node is the second end-point of the edge
    shared with the root node, then do not add the edge to the
    list as it will be added in the end (avoid duplicate edges
    in the list) */
    if (node != trav_buffer[trav_buffer_size - 1]->back)
    {
      *branches++ = node->length;
      *pmatrix_indices++ = node->pmatrix_index;
      *matrix_count = *matrix_count + 1;
    }

    if (node->next)
    {
      ops[*ops_count].parent_clv_index = node->clv_index;
      ops[*ops_count].parent_scaler_index = node->scaler_index;

      ops[*ops_count].child1_clv_index = node->next->back->clv_index;
      ops[*ops_count].child1_scaler_index = node->next->back->scaler_index;
      ops[*ops_count].child1_matrix_index = node->next->back->pmatrix_index;

      ops[*ops_count].child2_clv_index = node->next->next->back->clv_index;
      ops[*ops_count].child2_scaler_index = node->next->next->back->scaler_index;
      ops[*ops_count].child2_matrix_index = node->next->next->back->pmatrix_index;

      *ops_count = *ops_count + 1;
    }
  }
}

PLL_EXPORT int pll_utree_every(pll_utree_t * tree,
                               int (*cb)(const pll_utree_t *,
                                         const pll_unode_t *))
{
  unsigned int i;
  int rc = 1;

  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
    rc &= cb(tree, tree->nodes[i]);

  return (rc ? PLL_SUCCESS : PLL_FAILURE);
}

PLL_EXPORT int pll_utree_every_const(const pll_utree_t * tree,
                                     int (*cb)(const pll_utree_t *,
                                               const pll_unode_t *))
{
  unsigned int i;
  int rc = 1;

  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
    rc &= cb(tree, tree->nodes[i]);

  return (rc ? PLL_SUCCESS : PLL_FAILURE);
}

static void utree_traverse_recursive(pll_unode_t * node,
                                     int traversal,
                                     int (*cbtrav)(pll_unode_t *),
                                     unsigned int * index,
                                     pll_unode_t ** outbuffer)
{
  if (!cbtrav(node))
    return;

  if (traversal == PLL_TREE_TRAVERSE_PREORDER)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }

  if (node->next)
  {
    pll_unode_t * snode = node->next;
    do
    {
      utree_traverse_recursive(snode->back, traversal, cbtrav, index, outbuffer);
      snode = snode->next;
    }
    while (snode && snode != node);
  }

  if (traversal == PLL_TREE_TRAVERSE_POSTORDER)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }
}

PLL_EXPORT int pll_utree_traverse(pll_unode_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_unode_t *),
                                  pll_unode_t ** outbuffer,
                                  unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->next) return PLL_FAILURE;

  if (traversal == PLL_TREE_TRAVERSE_POSTORDER ||
      traversal == PLL_TREE_TRAVERSE_PREORDER)
  {

    /* we will traverse an unrooted tree in the following way

                2
              /
        1  --*
              \
                3

       at each node the callback function is called to decide whether we
       are going to traversing the subtree rooted at the specific node */

    utree_traverse_recursive(root->back, traversal, cbtrav, trav_size, outbuffer);
    utree_traverse_recursive(root, traversal, cbtrav, trav_size, outbuffer);
  }
  else
  {
    snprintf(pll_errmsg, 200, "Invalid traversal value.");
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

/* a callback function for checking tree integrity */
static int cb_check_integrity_mult(const pll_utree_t * tree,
                                   const pll_unode_t * node)
{
  unsigned int clv_index = node->clv_index;
  int scaler_index = node->scaler_index;
  unsigned int pmatrix_index = node->pmatrix_index;
  char * label = node->label;
  double length = node->length;
  unsigned int subnodes = 1;

  /* edge attributes */
  if (node->back->length != length)
  {
    snprintf(pll_errmsg, 200, "Inconsistent branch lengths: %lf != %lf",
             length, node->back->length);
    return PLL_FAILURE;
  }

  if (node->back->pmatrix_index != pmatrix_index)
  {
    snprintf(pll_errmsg, 200, "Inconsistent pmatrix indices: %u != %u",
             pmatrix_index, node->back->pmatrix_index);
    return PLL_FAILURE;
  }

  if (node->next)
  {
    /* node attributes */
    pll_unode_t * snode = node->next;
    do
    {
      subnodes++;

      if (tree->binary && subnodes > 3)
      {
        snprintf(pll_errmsg, 200, "Multifurcation found in a binary tree "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (subnodes > tree->tip_count)
      {
        snprintf(pll_errmsg, 200, "Multifurcation exceeding the tree size found "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (snode->clv_index != clv_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent CLV indices: %u != %u",
                 clv_index, snode->clv_index);
        return PLL_FAILURE;
      }
      if (snode->scaler_index != scaler_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent scaler indices: %u != %u",
                 scaler_index, snode->scaler_index);
        return PLL_FAILURE;
      }
      if (snode->label != label)
      {
        snprintf(pll_errmsg, 200, "Inconsistent node labels: '%s' != '%s'",
                 label, snode->label);
        return PLL_FAILURE;
      }
      if (!snode->next)
      {
        snprintf(pll_errmsg, 200, "Open roundabout (node->next is NULL) "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }
      snode = snode->next;
    }
    while (snode != node);
  }

  return 1;
}


PLL_EXPORT int pll_utree_check_integrity(const pll_utree_t * tree)
{
  return pll_utree_every_const(tree, cb_check_integrity_mult);
}

/* TODO: Memory allocation checks were not implemented in this function!!! */
static pll_unode_t * clone_node(const pll_unode_t * node)
{
  pll_unode_t * new_node = (pll_unode_t *)malloc(sizeof(pll_unode_t));
  memcpy(new_node, node, sizeof(pll_unode_t));

  if (node->label)
  {
    new_node->label = (char *)malloc(strlen(node->label)+1);
    strcpy(new_node->label,node->label);
  }

  if (node->next)
  {
    pll_unode_t * snode = node->next;
    pll_unode_t * new_snode = new_node;
    do
    {
      new_snode->next = (pll_unode_t *)malloc(sizeof(pll_unode_t));
      memcpy(new_snode->next, snode, sizeof(pll_unode_t));
      new_snode->next->label = new_node->label;
      snode = snode->next;
      new_snode = new_snode->next;
    }
    while (snode != node);

    new_snode->next = new_node;
  }

  return new_node;
}

static void utree_recurse_clone(pll_unode_t * new_root, const pll_unode_t * root)
{
  const pll_unode_t * node = root->back;
  if (node)
  {
    new_root->back = clone_node(node);
    new_root->back->back = new_root;

    if (node->next)
    {
      pll_unode_t * snode = node->next;
      pll_unode_t * new_snode = new_root->back->next;
      do
      {
        utree_recurse_clone(new_snode, snode);
        snode = snode->next;
        new_snode = new_snode->next;
      }
      while (snode && snode != node);
    }
  }
}

PLL_EXPORT pll_unode_t * pll_utree_graph_clone(const pll_unode_t * root)
{
  pll_unode_t * new_root = clone_node(root);

  const pll_unode_t * snode = root;
  pll_unode_t * new_snode = new_root;
  do
  {
    utree_recurse_clone(new_snode, snode);
    snode = snode->next;
    new_snode = new_snode->next;
  }
  while (snode && snode != root);

  return new_root;
}

PLL_EXPORT pll_utree_t * pll_utree_clone(const pll_utree_t * tree)
{
   /* choose the last inner node as the starting point of the clone. It does not
     really matter which node to choose, but since the newick parser places the
     root node at the end of the list, we use the same notation here */
  pll_unode_t * root = pll_utree_graph_clone(tree->vroot);

  if (tree->binary)
    return pll_utree_wraptree(root, tree->tip_count);
  else
    return pll_utree_wraptree_multi(root, tree->tip_count, tree->inner_count);
}

static pll_unode_t * rtree_unroot(pll_rnode_t * root, pll_unode_t * back)
{
  pll_unode_t * uroot = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->back = back;
  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
  uroot->length = uroot->back->length;

  if (!root->left)
  {
    uroot->next = NULL;
    return uroot;
  }

  uroot->next = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot->next)
  {
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot->next->next)
  {
    free(uroot->next);
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next->next = uroot;

  uroot->next->length = root->left->length;
  uroot->next->back = rtree_unroot(root->left, uroot->next);
  uroot->next->next->length = root->right->length;
  uroot->next->next->back = rtree_unroot(root->right, uroot->next->next);

  return uroot;
}

PLL_EXPORT pll_utree_t * pll_rtree_unroot(pll_rtree_t * tree)
{
  pll_rnode_t * root = tree->root;

  if (!root->left->left && !root->right->left)
  {
    pll_errno = PLL_ERROR_TREE_CONVERSION;
    snprintf(pll_errmsg,
             200,
             "Tree requires at least three tips to be converted to unrooted");
    return NULL;
  }

  pll_rnode_t * new_root;

  pll_unode_t * uroot = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot->next)
  {
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next = (void *)calloc(1,sizeof(pll_unode_t));
  if (!uroot->next->next)
  {
    free(uroot->next);
    free(uroot);
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  uroot->next->next->next = uroot;
  uroot->length = root->left->length + root->right->length;

  /* get the first root child that has descendants and make  it the new root */
  if (root->left->left)
  {
    new_root = root->left;
    uroot->back = rtree_unroot(root->right,uroot);
    /* TODO: Need to clean uroot in case of error */
    if (!uroot->back) return NULL;
  }
  else
  {
    new_root = root->right;
    uroot->back = rtree_unroot(root->left,uroot);
    /* TODO: Need to clean uroot in case of error*/
    if (!uroot->back) return NULL;
  }

  uroot->label = (new_root->label) ? xstrdup(new_root->label) : NULL;

  uroot->next->label = uroot->label;
  uroot->next->length = new_root->left->length;
  uroot->next->back = rtree_unroot(new_root->left, uroot->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->back) return NULL;

  uroot->next->next->label = uroot->label;
  uroot->next->next->length = new_root->right->length;
  uroot->next->next->back = rtree_unroot(new_root->right, uroot->next->next);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->next->back) return NULL;

  return pll_utree_wraptree(uroot,0);
}

PLL_EXPORT void pll_utree_create_pars_buildops(pll_unode_t * const* trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count)
{
  const pll_unode_t * node;
  unsigned int i;

  *ops_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    if (node->next)
    {
      ops[*ops_count].parent_score_index = node->clv_index;
      ops[*ops_count].child1_score_index = node->next->back->clv_index;
      ops[*ops_count].child2_score_index = node->next->next->back->clv_index;

      *ops_count = *ops_count + 1;
    }
  }
}
