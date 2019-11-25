Checklist for creating a new partition
================================================================================

- Create the partition
- Set the substitution parameters
- Set the tip states
  - Compress the MSA first
  - Set the pattern weights
- Set the frequencies
- Set Pinvar
- Set category rates
- Set category weights

Concepts
================================================================================

A partition is a section of the genome for which all sites evolved under the
same model. Informally, one can think of a partition being a gene, but
understand that this is not always the case. In particular partitions of a
genome might not code for something, and there are no requirements that
partitions are even _contiguous_.

Structure
================================================================================

```
typedef struct pll_partition
{
  unsigned int tips;
  unsigned int clv_buffers;
  unsigned int nodes; // tips + clv_buffer
  unsigned int states;
  unsigned int sites;
  unsigned int pattern_weight_sum;
  unsigned int rate_matrices;
  unsigned int prob_matrices;
  unsigned int rate_cats;
  unsigned int scale_buffers;
  unsigned int attributes;

  /* vectorization options */
  size_t alignment;
  unsigned int states_padded;

  double ** clv;
  double ** pmatrix;
  double * rates;
  double * rate_weights;
  double ** subst_params;
  unsigned int ** scale_buffer;
  double ** frequencies;
  double * prop_invar;
  int * invariant;
  unsigned int * pattern_weights;

  int * eigen_decomp_valid;
  double ** eigenvecs;
  double ** inv_eigenvecs;
  double ** eigenvals;

  /* tip-tip precomputation data */
  unsigned int maxstates;
  unsigned char ** tipchars;
  unsigned char * charmap;
  double * ttlookup;
  pll_state_t * tipmap;

  /* ascertainment bias correction */
  int asc_bias_alloc;
  int asc_additional_sites; // partition->asc_bias_alloc ? states : 0

  /* site repeats */
  struct pll_repeats *repeats;
} pll_partition_t;
```

This is generally organized into the following sections:

- Buffer counts and sizes
- Memory alignment state variables
- Model parameter buffers
- Buffers for the eigen[values|vectors]
- Sequence meta information
- Runtime flags
- Site repeat buffers

## Buffer Counts and Sizes

- `tips`: The number of tips present in this partition.
- `clv_buffers`: The number of CLVs, generally the number of edges in the tree.
- `nodes`: Node count. This includes the tips count.
- `states`: The number of states of the underlying data.
- `rate_matrices`: How many rate matrices present in the partition.
- `prob_matrices`: The number of edges in the tree. Practically this is the
  number of edges
- `scale_buffers`: The number of scale buffers.

## Memory Alignment State Variables

- `alignment`: One of three constants depending on what architecture is being
  used. At the time of writing these are:
    - `PLL_ALIGNMENT_CPU`
    - `PLL_ALIGNMENT_SSE`
    - `PLL_ALIGNMENT_AVX`
- `states_padded`: How many states are used, padded up to be aligned for the
  architecture flag that is used.

## Attributes

The `attributes` field is a bitset that has the is combination of the following
flags.

- Architecture attributes: Only one may be set
  - `PLL_ATTRIB_ARCH_CPU`
  - `PLL_ATTRIB_ARCH_SSE`
  - `PLL_ATTRIB_ARCH_AVX`
  - `PLL_ATTRIB_ARCH_AVX2`
- Ascertainment Bias: Only one of the "types" may be set, and if they are,
  `PLL_ATTRIB_AB_FLAG` must also be set
  - `PLL_ATTRIB_AB_LEWIS`
  - `PLL_ATTRIB_AB_FELSENSTEIN`
  - `PLL_ATTRIB_AB_STAMATAKIS`
  - `PLL_ATTRIB_AB_FLAG`
- Scalers
  - `PLL_ATTRIB_RATE_SCALERS`
- Optimizations: Only one may be set
  - `PLL_ATTRIB_PATTERN_TIP`
  - `PLL_ATTRIB_SITE_REPEATS`

## `clv`

In any given run, this will probably be the biggest allocation of memory, and
will be involved in almost all of the computation that `libpll` performs. As
such, it can be useful to understand what CLVs are, and how the `clv` buffer
relates to them.

A conditional likelihood vector (CLV) is an intermediate calculation which
informally represents the likelihood of a subtree. Conceptually, every node (not
a `pll_unode_t`) has a CLV, which is oriented with respect to the virutal root.

Notable Functions
================================================================================

    pll_partition_t* pll_partition_create(unsigned int tips,
                                          unsigned int clv_buffers,
                                          unsigned int states,
                                          unsigned int sites,
                                          unsigned int rate_matrices,
                                          unsigned int prob_matrices,
                                          unsigned int rate_cats,
                                          unsigned int scale_buffers,
                                          unsigned int attributes);

- `tips`: The number of tips of the tree. In phylogenetic terms, this is the
  number of taxa.
- `clv_buffers`: This is the number of CLVs that will be required to compute the
  tree. Practically, this is the number of edges. The number of rate categories
  is automatically accounted for, so no need to add it.
- `states`: The number of states that the model has, ie the type of sequence
  data that is being worked on. Practically, this is:
    - 4 for nucleotide data,
    - 20 for nucleotide data,
    - and 61 for codon data.
- `sites`: how long is the alignment. Note that, this is going to be the post
  compressed length of the sequence, i.e. the length that is from
  `pll_compress_site_patterns`.
- `rate_matrices`: The number of rate matrices that are allocated. In a
  simple and standard model, this is 1. In the case of a mixture model, this
  should be equal to the number of classes.
- `prob_matrices`: The number of probabilty matrices need for calculation. This
  will almost always be equal to the number of branches in the tree. In
  particular, the number of rate categories is automatically accounted for.
- `rate_cats`: Number of different rate categories to consider.
- `scale_buffers`: Number of scalling buffers to allocate. Practially, this is
  equal to the number of inner nodes in the tree.
- `attributes`: A bitflag set which controls several features of runtime.
    - Architecture Flags: The first category of attributes are the architecture
    flags, ie whether to use avx/2, sse3, or none.
      - `PLL_ATTRIB_ARCH_CPU`
      - `PLL_ATTRIB_ARCH_SSE`
      - `PLL_ATTRIB_ARCH_AVX`
      - `PLL_ATTRIB_ARCH_AVX2`
    - Optimization Flags: If these are set, then some optimizations, which are
    not always free, are used.
      - `PLL_ATTRIB_PATTERN_TIP`
      - `PLL_ATTRIB_SITE_REPEATS`
    - Bias correction flags: Flags that control the different methods of
    compute the bias correction
      - `PLL_ATTRIB_AB_LEWIS`
      - `PLL_ATTRIB_AB_FELSENSTEIN`
      - `PLL_ATTRIB_AB_STAMATAKIS`
      - `PLL_ATTRIB_AB_FLAG`
    - Misc: Everything else.
      - `PLL_ATTRIB_RATE_SCALERS`

----

    PLL_EXPORT void pll_partition_destroy(pll_partition_t * partition);

Destroys the given partition, deallocating the memory.

----

Together, these two functions are used to compress the MSA

    unsigned int* pll_compress_site_patterns(char** sequence,
                                       const pll_state_t* map,
                                       int count,
                                       int* length);

- `sequence`: the sequences, as represented in a `pll_msa_t`.
- `map`: A predefined map from `char` to `int`. The choices are:
    - `pll_map_bin`: For an alphabet of ${0,1}$.
    - `pll_map_nt`: For nucleotide data.
    - `pll_map_aa`: For amino acid data.
- `count`: The number of sequences in the alignment.
- `length`: This is both a in parameter, as well as an out parameter.
    Initially, this is the length of the sequences. After compression, it
    returns the new length via this parameter.


    void pll_set_pattern_weights(pll_partition_t* partition,
                                 const unsigned int* pattern_weights);

- `partition`: The partition that we want to set the weights for
- `pattern_weights`: The array of weights. While you could set this yourself,
    in general it should be the output of `pll_compress_site_patterns`.

----

    int pll_set_tip_states(pll_partition_t* partition,
                           unsigned int tip_index,
                           const pll_state_t map,
                           const char* sequence);

Where the arguements are

- `partition`: The partition to set the tip states for.
- `tip_index`: The index of the tip to set.
- `map`: A predefined map from `char` to `int`. The choices are:
    - `pll_map_bin`: For an alphabet of ${0,1}$.
    - `pll_map_nt`: For nucleotide data.
    - `pll_map_aa`: For amino acid data.
- `sequence`: The sequence associated with that tip.

-----

    void pll_set_frequencies(pll_partition_t* partition,
                             unsigned int params_index,
                             const double* params);

- `partition`: The partition to set the frequencies for.
- `params_index`: The model index to set the frequencies for.
- `params`: The values to set the frequencies to. This needs to be as long as
    the number of states in the model. These should add up to 1.

If you are using a model variant where we want to use emperical frequencies, we
can use

    double* pllmod_msa_empirical_frequencies(pll_partition_t* partition);

which will calculate the frequencies for the given partition.

----

    pll_set_subst_params(pll_partition_t* partition,
                         unsigned int params_index,
                         const double* params);

If we wanted to set substitution parameters of the model to be equal to the
following matrix

$$
\begin{bmatrix}
 * & a & b & c \\
 a & * & d & e \\
 b & d & * & f \\
 c & e & f & *
\end{bmatrix}
$$

then the following list will do the trick:


    double subst_params[] = {a, b, c, d, e, f}
