[This](This) will cover both the structures `pll_utree_t` and `pll_unode_t`, as well as
the concepts and design around the tree data structure in [`libpll`](libpll)

Concepts
================================================================================

Throughout this documentation, there is an intended difference between the terms
`pll_unode_t` and node. One refers to the actual datastructure that is used to
represent a vertex in a phylogenetic tree, while a node refers to that actual
vertex.  Each node in a phylogenetic tree is made up of 1 or more
`pll_unode_t`s. 

Each `pll_unode_t` is a collection of two pointers and some associated data. The
two pointers are the `next` and `back` pointers. If a `pll_unode_t` is an inner
node, then the `next` pointer points the the next `pll_unode_t` in the node. If
the `next` pointer is `null`, then that `pll_unode_t` represents a tip.

The `back` pointer represents edges. It points to a `pll_unode_t` accociated
with another node. Suppose that we have the tree `((a,b),c,d)`, then the libpll
representation of that would what is shown in the following figure.

![test](images/libpll_utree_figure.png)

Here, the dotted arcs indicate `next` pointers, and the solid lines represent
`back` pointers.

Structures
================================================================================

The datastructure is made up of two different structs. The first, `pll_utree_t`
wraps the tree. In general, when a tree is used for a function, it requires a
`pll_utree_t`. Some important things to know about this structure: the first
`inner_count` nodes in the `nodes` array are assumed to be "inner nodes". This
means that they have a non-null `next` pointer. Several functions that use
`pll_utree_t`s don't check for this, so they may fail when this assumption is
violated. To avoid this, use the [`pll_utree_wraptree`](#Notable-Functions)
function discussed below to create a `pll_utree_t`.

```
typedef struct pll_utree_s
{
  unsigned int tip_count;
  unsigned int inner_count;
  unsigned int edge_count;
  int binary;

  pll_unode_t ** nodes;
  pll_unode_t * vroot;
} pll_utree_t;
```

Fields:

- `tip_count`: Number of tips
- `inner_count`: number of inner nodes. Not the number of `pll_unode_t` that are
  part of inner nodes, but the number of interior nodes of a phylogenetic tree.
- `edge_count`: the number of edges in the tree.
- `nodes`: an array of pointers to nodes.
- `vroot`: A pointer to the virutal root. Is it aways an inner node?

----

It might be useful to at the above figure when looking at this structure, since
most of the complexity is in the behavior of the `next` and `back` pointers.

```
typedef struct pll_unode_s
{
  char * label;
  double length;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  struct pll_unode_s * next;
  struct pll_unode_s * back;

  void * data;
} pll_unode_t;
```

Fields:

- `label`: The label of the node. Optional.
- `length`: The length of the edge represented by the `back` pointer.
- `node_index`: Index of this node in the `nodes`.
- `clv_index`: Index of the CLVs to use when calculating a likelihood
- `scaler_index`: Index into the scaler array to represent the CLV scaler
- `pmatrix_index`: index into the array of pmatrices. These pmatrices need to be
  computed based on the length of the branch
- `next`, `back`: See the explaination in the concepts section.
- `data`: An extra pointer to store "user data". In practice, this can be used
  for any task, but exisiting functions might also use it, so be careful.


Notable Functions
================================================================================

```
PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename);
```

```
PLL_EXPORT pll_utree_t * pll_utree_parse_newick_unroot(const char * filename);
```

These functions will create a `pll_utree_t` from a newick _file_. If there is an
error, the function will return `PLL_ERROR` (which happens to be 0). If a rooted
tree is passed to the unroot version, then the tree is unrooted after parsing
(the root is suppressed).


-------

```
PLL_EXPORT void pll_utree_destroy(pll_utree_t * tree,
                                  void (*cb_destroy)(void *));
```

Deallocate the memory associated with a utree. `cb_destroy` is used to delete
the user data allocated in `data`.

-------

```
PLL_EXPORT pll_utree_t * pll_utree_wraptree(pll_unode_t * root,
                                            unsigned int tip_count);
```

Takes a tree, represented by a node, and optionally a tip count. Will produce a
`pll_utree_t` that contains that tree. The pointer to the original node is not
invalidated.

------

```
PLL_EXPORT int pll_utree_traverse(pll_unode_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_unode_t *),
                                  pll_unode_t ** outbuffer,
                                  unsigned int * trav_size);
```

Creates a list of nodes from a traversal, starting at `root`. The order of the
traversal can be controlled with the `traversal` arguement, which accepts either
`PLL_TREE_TRAVERSE_POSTORDER` or `PLL_TREE_TRAVERSE_PREORDER`. The callback
function controls which nodes are traversed by returning `PLL_SUCESS` for a
node node which should be added to the outbuffer. If `PLL_FAILURE` is returned
instead, then the traversal is halted for that subtree, and the node which
returned it is not added to the outbuffer. `trav_size` is an outparameter.
Returns `PLL_SUCCESS` on a traversal without errors, and `PLL_FAILURE` if there
was an error.

------

```
PLL_EXPORT void pll_utree_create_operations(pll_unode_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count);
```

Given the `pll_unode_t**` from a traversal, this will create a list of
[`pll_operation_t`](pll_operation_t) from that list. If `branches`,
`pmatrix_index` is not `null`, then the values of the branches and probability
matrix indices are stored in these buffers. Likewise if `matrix_count` is not
null, then it is an out paramter with the number of matricies needed for a
traversal. The remaining two parameters are both out parameters, the list of
operations, `ops` and the number of operations `ops_count`.

If a full traversal is desired, the callback should return 1 for all inputs.

------

```
PLL_EXPORT pll_utree_t * pll_utree_clone(const pll_utree_t * root);
```

Clones a tree. This is a semi-deep copy. The feilds `label` and a nodes `next`
and `back` pointers are deep copied, but the `data` field is just a shallow
copy.
