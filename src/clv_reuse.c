#include "pll.h"
#include <stdarg.h>
#include <search.h>
#include <time.h>

//const unsigned int EMPTY_ELEMENT = (unsigned int) -1;

time_t t;

/*
With this code we re-use (overwrite) the tip CLVs for likelihood
computations of full tree traversals.Most functions of the code
are simply copies of the original with the addition of overwriting
CLVs. The pll_clv_reuse function follows the motif of rooted.c example.
*/


/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_unode_t * node){
  return 1;
}

/* branch lengths not present in the newick file get a value of 0.000001 */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  unsigned int i;

  for (i = 0; i < tree->tip_count; ++i)
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;

  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
  {
    if (!tree->nodes[i]->length)
      tree->nodes[i]->length = length;
    if (!tree->nodes[i]->next->length)
      tree->nodes[i]->next->length = length;
    if (!tree->nodes[i]->next->next->length)
      tree->nodes[i]->next->next->length = length;
  }
}

/* ---------------------------------------------------- */

void print_clvs(pll_partition_t * partition, short a){
  unsigned i = 0;
  for (i = 0; i < 7; i++){ //partition -> nodes
    printf ("Tip %d: ",i);
    pll_show_clv(partition,i,PLL_SCALE_BUFFER_NONE, 4);
  }
}

static void pll_case_tiptip_clv_reuse(pll_partition_t * partition,
                        const pll_operation_t * op)
{
    const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
    const double * right_matrix = partition->pmatrix[op->child2_matrix_index];
    double * parent_clv = partition->clv[op->child1_clv_index];
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

static void pll_case_tipinner_clv_reuse(pll_partition_t * partition,
                          const pll_operation_t * op)
{
  unsigned x = partition->clv_map[op->parent_clv_index];
  if ( x > partition->tips)
    x = partition->clv_map[x];
  partition->clv_map[op->parent_clv_index] = x;
  double * parent_clv = partition->clv[x];
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

  pll_core_update_partial_ti(partition->states,
                             sites,
                             partition->rate_cats,
                             parent_clv,  //different that the original version
                             parent_scaler,
                             partition->tipchars[tip_clv_index],
                             partition->clv[inner_clv_index],
                             partition->pmatrix[tip_matrix_index],
                             partition->pmatrix[inner_matrix_index],
                             right_scaler,
                             partition->tipmap,
                             partition->maxstates,
                             partition->attributes);

}

static void pll_case_innerinner_clv_reuse(pll_partition_t * partition,
                            const pll_operation_t * op)
{
  const double * left_matrix = partition->pmatrix[op->child1_matrix_index];
  const double * right_matrix = partition->pmatrix[op->child2_matrix_index];

  /* ------------------- start of differences ------------------- */
  unsigned x = partition->clv_map[op->parent_clv_index];
  while ( x >= partition->tips)
    x = partition->clv_map[x];
  partition->clv_map[op->parent_clv_index] = x;

  double * parent_clv = partition->clv[partition->clv_map[op->parent_clv_index]]; //1
  double * left_clv = partition->clv[partition->clv_map[op->child1_clv_index]]; //2
  double * right_clv = partition->clv[partition->clv_map[op->child2_clv_index]]; //3
  /* ------------------------------------------------------------ */

  unsigned int * parent_scaler;
  unsigned int * left_scaler;
  unsigned int * right_scaler;
  unsigned int sites = partition->sites;
  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  /* get scalers */
  if (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[op->parent_scaler_index];
  if (op->child1_scaler_index != PLL_SCALE_BUFFER_NONE)
    left_scaler = partition->scale_buffer[op->child1_scaler_index];
  else
    left_scaler = NULL;
  if (op->child2_scaler_index != PLL_SCALE_BUFFER_NONE)
    right_scaler = partition->scale_buffer[op->child2_scaler_index];
  else
    right_scaler = NULL;

  pll_core_update_partial_ii(partition->states,
                             sites,
                             partition->rate_cats,
                             parent_clv,  //1
                             parent_scaler,
                             left_clv,  //2
                             right_clv, //3
                             left_matrix,
                             right_matrix,
                             left_scaler,
                             right_scaler,
                             partition->attributes);
}

void pll_update_partials_clv_reuse(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count,
                                    unsigned int update_repeats)
{
  unsigned int i;
  const pll_operation_t * op;

  for (i = 0; i < count; ++i)
  {
    op = &(operations[i]);
    if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    {
      if ((op->child1_clv_index < partition->tips) &&
          (op->child2_clv_index < partition->tips))
      {
        printf(" tiptip\n");
        pll_case_tiptip_clv_reuse(partition,op);
      }
      else if ((operations[i].child1_clv_index < partition->tips)
                || (operations[i].child2_clv_index < partition->tips))
      {
        printf(" tipinner 1\n");
        pll_case_tipinner_clv_reuse(partition,op);
      }
      else
      {
        printf(" innerinner\n");
        pll_case_innerinner_clv_reuse(partition,op);
      }
    }
    else
      pll_case_innerinner_clv_reuse(partition,op);
  }
}

void pll_utree_create_operations_clv_reuse(pll_unode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count,
                                            pll_partition_t * partition)
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

      partition->clv_map[node->clv_index] = node->next->back->clv_index; //map the left child

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

void pll_dealloc_partition_data_clv_reuse(pll_partition_t * partition)
{
  unsigned int i;

  if (!partition) return;

  free(partition->rates);
  free(partition->rate_weights);
  free(partition->eigen_decomp_valid);
  if (partition->prop_invar)
    free(partition->prop_invar);
  if (partition->invariant)
    free(partition->invariant);
  if (!partition->pattern_weights)
    free(partition->pattern_weights);

  if (partition->scale_buffer)
    for (i = 0; i < partition->scale_buffers; ++i)
      free(partition->scale_buffer[i]);
  free(partition->scale_buffer);

  if (partition->tipchars)
    for (i = 0; i < partition->tips; ++i)
      pll_aligned_free(partition->tipchars[i]);
  free(partition->tipchars);

  if (partition->ttlookup)
    pll_aligned_free(partition->ttlookup);

  if (partition->charmap)
    free(partition->charmap);

  if (partition->tipmap)
    free(partition->tipmap);

  if (partition->clv)
  {
    int start = (partition->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                    partition->tips : 0;
    for (i = start; i < partition->tips; ++i)
      pll_aligned_free(partition->clv[i]);
  }
  free(partition->clv);

  if (partition->pmatrix)
  {
    //for (i = 0; i < partition->prob_matrices; ++i)
      pll_aligned_free(partition->pmatrix[0]);
  }
  free(partition->pmatrix);

  if (partition->subst_params)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->subst_params[i]);
  free(partition->subst_params);

  if (partition->eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->eigenvecs[i]);
  free(partition->eigenvecs);

  if (partition->inv_eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->inv_eigenvecs[i]);
  free(partition->inv_eigenvecs);

  if (partition->eigenvals)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->eigenvals[i]);
  free(partition->eigenvals);

  if (partition->frequencies)
    for (i = 0; i < partition->rate_matrices; ++i)
      pll_aligned_free(partition->frequencies[i]);
  free(partition->frequencies);

  if (partition->pattern_weights)
    free(partition->pattern_weights);

  if (partition->repeats)
  {
    pll_repeats_t *repeats = partition->repeats;
    for (i = 0; i < partition->nodes; ++i)
    {
      free(repeats->pernode_site_id[i]);
      free(repeats->pernode_id_site[i]);
    }
    free(repeats->pernode_site_id);
    free(repeats->pernode_id_site);
    free(repeats->pernode_ids);
    free(repeats->perscale_ids);
    free(repeats->pernode_allocated_clvs);
    free(repeats->lookup_buffer);
    free(repeats->toclean_buffer);
    free(repeats->id_site_buffer);
    free(repeats->bclv_buffer);
    free(repeats);
  }

  free(partition -> clv_map);

  free(partition);
}

/* command line arguments: 1)newick_tree, 2)msa, 3)rooted, 4)version, 5)kernel */
int pll_clv_reuse(int argc, char * input_tree, char * input_msa,
            short int is_rooted, unsigned clv_re, unsigned input_ker)// char * argv[]) /* since this is just a test code, I assume arguments/files etc are correct */
{
  /* general variable initialization */
  unsigned int i = 0;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  unsigned int matrix_count, ops_count;
  unsigned int * matrix_indices;
  double * branch_lengths;
  pll_partition_t * partition;
  pll_operation_t * operations;
  pll_unode_t ** travbuffer;

  /* parse the unrooted binary tree in newick format, and store the number
     of tip nodes in tip_nodes_count */
  pll_utree_t * tree = pll_utree_parse_newick(input_tree);
  if (!tree)
    assert(0);

  tip_nodes_count = tree->tip_count;
  set_missing_branch_length(tree, 0.000001);
  /* compute and show node count information */
  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(tip_nodes_count *
                                               sizeof(unsigned int));
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = tree->nodes[i]->clv_index;
    ENTRY entry;
    entry.key = tree->nodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  /* read PHYLIP alignment */
  pll_phylip_t * fd = pll_phylip_open(input_msa, pll_map_phylip);
  pll_msa_t * msa = pll_phylip_parse_interleaved(fd);
  pll_phylip_close(fd);

  partition = pll_partition_create(tip_nodes_count,
                                   inner_nodes_count,
                                   4,
                                   (unsigned int)(msa->length),
                                   1,
                                   branch_count,
                                   4,
                                   inner_nodes_count,
                                   input_ker);

  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[6] = {1,1,1,1,1,1};
  double rate_cats[4] = {0};

  pll_compute_gamma_cats(1, 4, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst_params);
  pll_set_category_rates(partition, rate_cats);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = msa->label[i];
    ENTRY * found = NULL;
    found = hsearch(query,FIND);
    unsigned int tip_clv_index = *((unsigned int *)(found->data));
    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, msa->sequence[i]);
  }

  pll_msa_destroy(msa);
  hdestroy();
  free(data);

  travbuffer = (pll_unode_t **)malloc(nodes_count * sizeof(pll_unode_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations = (pll_operation_t *)malloc(inner_nodes_count *
                                                sizeof(pll_operation_t));

  /* perform a postorder traversal of the unrooted tree */
  pll_unode_t * root = tree->nodes[tip_nodes_count+inner_nodes_count-1];
  unsigned int traversal_size;

  pll_utree_traverse(root,
                     PLL_TREE_TRAVERSE_POSTORDER,
                     cb_full_traversal,
                     travbuffer,
                     &traversal_size);

  /* thus far there are no differences between the versions */
  if (clv_re){
    partition->clv_map = calloc(sizeof(unsigned), traversal_size); /* tracks in which clv we have actually stored the values */
    for (i = 0; i < traversal_size; i++){
      partition->clv_map[i] = i; /* initially each index points to its own clv */
    }
  }

  if (clv_re){
    pll_utree_create_operations_clv_reuse(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count,
                              partition);
    int start = (partition->attributes & PLL_ATTRIB_PATTERN_TIP) ?
                    partition->tips : 0;
    for (i = start + partition->tips; i < partition->clv_buffers + partition->tips; ++i)
      pll_aligned_free(partition->clv[i]);
  }
  else
    pll_utree_create_operations(travbuffer, traversal_size, branch_lengths, matrix_indices, operations,&matrix_count, &ops_count);

  unsigned int params_indices[4] = {0,0,0,0};

  pll_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);


  /* after creating the operations we need to update */
  if (clv_re)
    pll_update_partials_clv_reuse(partition, operations, ops_count, 0);
  else
    pll_update_partials_rep(partition, operations, ops_count, 0);

  /* now the root neeeds to be calculated */
  double logl = 0.0;
  if (clv_re){
      logl = pll_compute_edge_loglikelihood(partition,
                                                partition->clv_map[root->clv_index],
                                                root->scaler_index,
                                                partition->clv_map[root->back->clv_index],
                                                root->back->scaler_index,
                                                root->pmatrix_index,
                                                params_indices,
                                                NULL);
  }
  else{
      logl = pll_compute_edge_loglikelihood(partition,
                                                 root->clv_index,
                                                 root->scaler_index,
                                                 root->back->clv_index,
                                                 root->back->scaler_index,
                                                 root->pmatrix_index,
                                                 params_indices,
                                                 NULL);
  }
  fprintf(stderr, "%lf\n", logl);

  /* free everything */
  if (clv_re)
    pll_dealloc_partition_data_clv_reuse(partition);
  else
    pll_partition_destroy(partition);
  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);
  pll_utree_destroy(tree, NULL);

  return 0;
}
