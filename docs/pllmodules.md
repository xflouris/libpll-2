pllmodules
===============================================================================

How to run an SPR round on a tree.
-------------------------------------------------------------------------------

1. Read a MSA in
2. Set up the `pll_partition_t`.
    1. (Optional) compress the site pattters
    2. Create partition with appropriate variables.
    3. Set the tip states.
    4. Set initial model parameters
        1. Frequencies
        2. Substitution parameters
        3. $p_{inv}$
        4. Rate categories, and associated rates
3. Make an initial tree
4. Make a `pllmod_treeinfo_t`
5. Initalize the partitions for the treeinfo struct

###1.0: Read MSA

There is a datastructure to hold an MSA. It has the following definition

```
struct pll_msa_s{
    int count;
    int length;

    char** sequence;
    char** label;
}pll_msa_t;
```

Where the parameters are as follows:

- `count`: The number of taxa in the MSA.
- `length`: the length of each sequence in the alignment. More litterally, the
    length of the character arrays indexed by `sequence`.
- `sequence`: The sequences. Not null terminated.
- `label`: sequence label. Is null terminated.

###2.0: Create a Partition

The information about a sequence that evolves along the same line is stored  in
a partition structure, as  well  as  the  model  parameters  relating  to  that
process.  Roughly, this means that each partition is representative of a  gene,
though that language is slightly problematic, since we  can  analyze  sequences
that   don't   code   for   anything,   and    thus    are    not    a    gene.

To accomplish steps 2.1 to 2.4, I would recommend reading
[pll_partition_t](pll_partition_t) for this.

####3.0: Make an initial tree

There are two easy methods to create an initial tree. The first is to make a
random tree using

```
pll_utree_t* pllmod_utree_create_random(unsigned int taxa_count,
                           const char* const* names,
                           unsigned int random_seed);
```

This will create a bifurcated tree at random. It does this by creating a
minimal 3 tip tree, selecting a branch at random, and then hanging a new tip
off of that. The algorithm repeats this until all the tips have been added to
the tree.

- `taxa_count`: The number of tips that the resulting tree will have.
- `names`: The labels that will be assigned to the tips. Represented as null
    terminated `char*`.
- `random_seed`: The seed to pass to the random number generator.

The other method is to use a parsimony tree. This can be created using

```
pll_utree_t* pllmod_utree_create_parsimony(unsigned int taxon_count,
                                           unsigned int seq_length,
                                           char* const* names,
                                           char* const* sequences,
                                           const unsigned int* site_weights,
                                           const pll_state_t* map,
                                           unsigned int states,
                                           unsigned int attributes,
                                           unsigned int random_seed,
                                           unsigned int* score);
```

This will create a tree under a parsimony model. The arguments are

- `taxon_count`: The number of tips on the tree.
- `seq_length`: The length of the sequences in the alignment.
- `names`: An array with length `taxon_count` of `char*` null terminated
    strings.
- `sequences`: An array with length `taxon_count` of `char*` that all have
    length `seq_length`.
- `site_weights`: An array that is `seq_length` long of site weights. Can be
    `nullptr` which means that all sites have equal weight.
- `map`: A predefined map from `char` to `int`. The choices are:
    - `pll_map_bin`: For an alphabet of ${0,1}$.
    - `pll_map_nt`: For nucleotide data.
    - `pll_map_aa`: For amino acid data.
- `states`: Number of states in the data.
- `attributes`: Exactly one of the following.
    - `PLL_ATTRIB_ARCH_CPU`: No special extensions,
    - `PLL_ATTRIB_ARCH_SSE`: Use the SSE3 extensions,
    - `PLL_ATTRIB_ARCH_AVX`: Use the AVX extensions,
    - `PLL_ATTRIB_ARCH_AVX2`: Use the AVX2 extensions,
    - `PLL_ATTRIB_ARCH_AVX512`: Use the AVX512 extensions,
    In addition, set the `PLL_ATTRIB_PATTERN_TIP` to enable site repeats.
- `random_seed`: Seed to pass to the random number generator
- `score`: output parameter that contains the score of the tree.

####4.0: Make a `pllmod_treeinfo_t`

Please see the relevant section in [pllmod_treeinfo_t](pllmod_treeinfo_t)

####5.0: Initialize the partitions for the treeinfo struct

Please see the relevant section in [pllmod_treeinfo_t](pllmod_treeinfo_t)
