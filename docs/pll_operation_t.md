Structures
================================================================================

A structure which simply keeps track of the "operations" required to calculate
the likelihood on a tree, during the partial likelihood calculation. In general,
these should not be produced by hand, but instead by
[`pll_utree_create_operations`](#Notable-Functions).

```
typedef struct pll_operation
{
  unsigned int parent_clv_index;
  int parent_scaler_index;
  unsigned int child1_clv_index;
  unsigned int child1_matrix_index;
  int child1_scaler_index;
  unsigned int child2_clv_index;
  unsigned int child2_matrix_index;
  int child2_scaler_index;
} pll_operation_t;
```

Notable Functions
================================================================================

```
PLL_EXPORT void pll_utree_create_operations(pll_unode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);
```

Creates a list of operations, given a list of nodes. These nodes should be
present in the `trav_buffer` in a post-order traversal. If `branches` or
`pmatrix_indices` are not `null`, then the nodes are assigned the values in
`branches` or `pmatrix_indices`. If `matrix_count` is not `null`, then it is an
out parameter which gives the number of matricies required to compute this
traversal. Finally, the operations are placed in `ops`, and the number of
operations is placed in `ops_count`.

----

```
PLL_EXPORT void pll_update_partials(pll_partition_t * partition,
                                    const pll_operation_t * operations,
                                    unsigned int count);
```

Computes the CLVS for all the trees in the nodes specified in the `operations`
list, where `operations` has length `count`.

----

```
PLL_EXPORT void pll_update_repeats(pll_partition_t * partition,
                    const pll_operation_t * op) ;
```

Calculates the repeats for `partition` for a particular operation.
