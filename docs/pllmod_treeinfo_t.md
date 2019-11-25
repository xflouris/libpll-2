pll_treeinfo_t
===============================================================================

Creating a treeinfo
-------------------------------------------------------------------------------

```
pllmod_teeinfo_t* pllmod_treeinfo_create(pll_unode_t* root,
                                         unsigned int tips,
                                         unsigned int partitions,
                                         int brlen_linkage);
```

- `root`: A pointer to the virtual root of the unrooted tree.
- `tips`: Number of tips in the tree. Almost always this will be the number of
    taxa under consideration.
- `partitions`: Number of partitions that are in the data set. These are
    initialize later.
- `brlen_linkage`: If the branch lengths are linked. This has the following
    constants assigned to it:
    - `PLLMOD_COMMON_BRLEN_UNLINKED`: The branch lengths have no relation to
        each other
    - `PLLMOD_COMMON_BRLEN_SCALED`: The branch lengths the product of a scalar
        and a global length.
    - `PLLMOD_COMMON_BRLEN_LINKED`: The branch lengths are linked together.

Initializing a partition
-------------------------------------------------------------------------------

```
int pllmod_treeinfo_init_partition(pllmod_treeinfo_t* treeinfo,
                                   unsigned int partition_index,
                                   pll_partition_t* partition,
                                   int params_to_optimize,
                                   int gamma_mode,
                                   double alpha,
                                   const unsigned int* param_indices,
                                   const int* subst_matrix_symmetries);
```

- `treeinfo`: the `pllmod_treeinfo_t` to use.
- `partition_index`: the index of this partition.
- `partition`: the partition itself. This partition should already be "setup"
    before it is added.
- `params_to_optimize`: Which parameters to optimize when performing
    optimizations. Can have the following values:
    - `PLLMOD_OPT_PARAM_ALL`: Optimize everything.
    - `PLLMOD_OPT_PARAM_SUBST_RATES`: Optimize the substitution rates.
    - `PLLMOD_OPT_PARAM_ALPHA`: Optimize the alpha parameter for a +G model.
    - `PLLMOD_OPT_PARAM_PINV`: Optimize the proportion of invariant sites.
    - `PLLMOD_OPT_PARAM_FREQUENCIES`: Optimize the initial frequencies.
    - `PLLMOD_OPT_PARAM_BRANCHES_SINGLE`: Optimize a single branch.
    - `PLLMOD_OPT_PARAM_BRANCHES_ALL`: Optimize all the branches, at once.
    - `PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE`: Optimize the branches, one by one.
    - `PLLMOD_OPT_PARAM_TOPOLOGY`: Optimize the topology of the tree.
    - `PLLMOD_OPT_PARAM_FREE_RATES`: No idea
    - `PLLMOD_OPT_PARAM_RATE_WEIGHTS`: Optimize the category rates, for a
        k-class gamma model.
    - `PLLMOD_OPT_PARAM_BRANCH_LEN_SCALAR`: Optimize the scalars for the branch
        lengths, in a scaled branch length scenario
    - `PLLMOD_OPT_PARAM_USER`: Run some user code
- `gamma_mode`: Controls the gamma rate discretion modes. Two options:
    - `PLL_GAMMA_RATES_MEAN`
    - `PLL_GAMMA_RATES_MEDIAN`
- `alpha`: the initial alpha to use for the model.
- `param_indices`: Specify the parameter indices to ... do stuff. Can be set to
    null, at which point the defaults are used.
- `subst_matrix_symmetries`: a list of symmetries in the model. Can be set to
    nullptr as well. If set to null, indicates that there are no symmetries
    in the model. If there are symmetries, the order is "left to right". For
    example, if our symmetry array was `{0,0,0,1,1,1}`, then the two symmetry
    groups would be the top row, and the reamaining triangle.

Performing a tree search
-------------------------------------------------------------------------------

To perform a tree search, it really is best to use the `pll_mod_spr_round`
function. This will look for good SPR moves to make, as well as optimize branch
lenghts. To optimize other model parameters, other functions should be used.
Here is the declaration of the function

```
double pllmod_algo_spr_round(pllmod_treeinfo_t* treeinfo,
                             int radius_min,
                             int radius_max,
                             int n_topologies,
                             int thorough,
                             int brlen_opt_method,
                             double bl_max,
                             double bl_min,
                             int smoothings,
                             double epsilon,
                             cutoff_info_t* cutoff_info,
                             double subtree_cutoff);
```

- `treeinfo`: A `treeinfo_t` struct that has been initalized using the methods
    above.
- `radius_min`: The minimum SPR radius to search.
- `radius_max`: The maximum SPR radius to search.
- `n_topologies`: The number of topologies to keep.
- `thorough`: Be thorough?
- `brlen_opt_method`: The optimization method to use when optimizing the branch
    lengths. Can be:
    - `PLLMOD_OPT_BLO_NEWTON_FAST`: This is the standard. Uses Newtons' method.
    - `PLLMOD_OPT_BLO_NEWTON_SAFE`: As fast, but with a per branch likelihood
        check.
    - `PLLMOD_OPT_BLO_NEWTON_FALLBACK`: Fast, but can fallback to safe
    - `PLLMOD_OPT_BLO_NEWTON_GLOBAL`: Newton, with additional searches for more
        optima
    - `PLLMOD_OPT_BLO_NEWTON_OLDFAST`: Old implementation of fast
    - `PLLMOD_OPT_BLO_NEWTON_OLDSAFE`: Old implementation of safe
- `bl_min`: The minimum branch length. Can be negative (i.e. the function
    works and won't behave badly in a operational sense. Branch lengths that
    are negative generally aren't meaningful in a biological sense).
- `bl_max`: The maximum branch length. Should be larger than `bl_min`.
- `smoothings`: Number of iterations for branch length optimization. Will
    operate if negative. I don't know what happens in that case.
- `epsilon`: The threshold to terminate optimization. If improvement falls
    below this, terminate. More commonly known as tolerance.
- `cutoff_info`: A struct that contains some cutoff information. It seems to be
  a return parameter.
- `subtree_cutoff`: Used to calculate a likelihood cutoff. A larger value measn
    that less subtrees are cut off.

Optimizing Model Parameters
-------------------------------------------------------------------------------

Model parameters need to be optimized outside of the SPR round searches. We can
use the following functions to perform the required operations. All of the
following functions find the parameter(s) which produce the highest likelihood
given a topology with branch lengths.

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_subst_rates_treeinfo(pllmod_treeinfo_t* treeinfo,
                                            unsigned int params_index,
                                            double min_rate,
                                            double max_rate,
                                            double bfgs_factor,
                                            double tolerance);
```

This function will optimize the substitution rates of the model. For example,
the $\alpha$ rate in JC69. The parameters are

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `params_index`: The index for the params. (For mixture models).
- `min_rate`: The minimum value for each substitution rate. Can be non-positive.
- `max_rate`: The maximum value for each substitution rate.
- `bfgs_factor`: Please see [bfgs_factor](bfgs_factor).
- `tolerance`: Please see [pg_tol](pg_tol).

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_frequencies_treeinfo(pllmod_treeinfo_t* treeinfo,
                                            unsigned int params_index,
                                            double min_freq,
                                            double max_freq,
                                            double bfgs_factor,
                                            double tolerance);
```

This function will optimize the initial frequencies, normally notated as $\pi$.

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `params_index`: The index for the params. (For mixture models).
- `min_freq`: The minimum value for each frequency.
- `max_freq`: The maximum value for each frequency.
- `bfgs_factor`: Please see [bfgs_factor](bfgs_factor).
- `tolerance`: Please see [pg_tol](pg_tol).

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_rates_weights_treeinfo(pllmod_treeinfo_t* treeinfo,
                                              double min_rate,
                                              double max_rate,
                                              double bfgs_factor,
                                              double tolerance);
```

TODO

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `min_rate`: The minimum value for each TODO rate.
- `max_rate`: The maximum value for each TODO rate.
- `bfgs_factor`: Please see [bfgs_factor](bfgs_factor).
- `tolerance`: Please see [pg_tol](pg_tol).

Seems to cause an assertion failure for me. I don't know if I am doing
something wrong here. Here is the assertion error:

```
hellopll: algo_search.c:1268: pllmod_algo_spr_round: Assertion `fabs(loglh - best_lh) < 1e-6' failed.
```

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_alpha_pinv_treeinfo(pllmod_treeinfo_t* treeinfo,
                                           unsigned int params_index,
                                           double min_alpha,
                                           double max_alpha,
                                           double min_pinv,
                                           double max_pinv,
                                           double bfgs_factor,
                                           double tolerance);
```

Optimize the $\alpha$ (i.e. the gamma parameter) and the proportion of
invariant sites.

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `params_index`: The index for the params. (For mixture models).
- `min_alpha`: The minimum value for alpha.
- `max_alpha`: The maximum value for alpha.
- `min_pinv`: The minimum value for the invariant site proportion.
- `max_pinv`: The maximum value for the invariant site proportion.
- `bfgs_factor`: Please see [bfgs_factor](bfgs_factor).
- `tolerance`: Please see [pg_tol](pg_tol).

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_brlen_scalers_treeinfo(pllmod_treeinfo_t* treeinfo,
                                              double min_scaler,
                                              double max_scaler,
                                              double min_brlen,
                                              double max_brlen,
                                              double lh_epsilon);
```

Optimize the branch lengths when using `PLLMOD_COMMON_BRLEN_SCALED`.

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `params_index`: The index for the params. (For mixture models).
- `min_scaler`: The minimum value for the scalar parameter.
- `max_scaler`: The maximum value for the scalar parameter.
- `min_brlen`: The minimum value for branch lengths.
- `max_brlen`: The maximum value for branch lengths.
- `lh_epsilon`: Tolerance for the optimization

-------------------------------------------------------------------------------

```
double pllmod_algo_opt_brlen_treeinfo(pllmod_treeinfo_t* treeinfo,
                                              double min_brlen,
                                              double max_brlen,
                                              double lh_epsilon,
                                              int max_iters);
```

Optimize the branch lengths of the tree. Normally, this is done during
`pllmod_algo_spr_round`, so this is not necessary. Regardless, it is here for
completion.

- `treeinfo`: a `pllmod_treeinfo_t` strut that contains the model and tree.
- `params_index`: The index for the params. (For mixture models).
- `min_brlen`: The minimum value for each substitution rate. Can be non-positive.
- `max_brlen`: The maximum value for each substitution rate.
- `lh_epsilon`: The tolerance for the optimization algorithm.
- `max_iters`: Maximum iterations to go through when iterating through branches
    and optimizing.

-------------------------------------------------------------------------------
