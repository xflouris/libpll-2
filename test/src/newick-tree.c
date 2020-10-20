// for asprintf
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include "common.h"

#include <stdarg.h>
#include <search.h>

#define NEWICK_BL_PREC 4

#define TREE_COUNT     5
#define TREEFILE_COUNT 5

const char * newick_tree_list[TREE_COUNT] =
    {
    "(A,B,(C,D));",
    "((A,B),C,(D,E));",
    "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5)0:0;",
    "((A:0.1,B:0.2,C:0.3):1,(D:0.4,(E:0.5,F):0.1):0.6);",
    "((taxon1:0.100,2:2)100,  (taxon3:0.2,\n(t4:0.5,t5:0.3)95:0.22):0.6);"
    };

const char * newick_file_list[TREEFILE_COUNT] =
    {
        "testdata/small.tree",
        "testdata/small.rooted.tree",
        "testdata/small.rooted.tip.tree",
        "testdata/medium.tree",
        "testdata/ribosomal_l5_pf00673.tree"
    };

char * newick_print_cb(const pll_unode_t * node)
{
  char * newick;
  size_t size_alloced = asprintf(&newick, "%s[%u]:%.*lf",
                                 node->label ? node->label : "" ,
                                 node->clv_index, NEWICK_BL_PREC, node->length);

  assert(size_alloced >= 0);

  return newick;
}

int test_tree(pll_utree_t * tree)
{
  int res;
  char * out_newick;

  printf("Tree: %s, virtual root at clv_id: %u\n",
         tree->binary ? "BINARY" : "MULTIFURCATING", tree->vroot->clv_index);
  printf("Number of tips/inner nodes/edges in tree: %d / %d / %d\n",
         tree->tip_count, tree->inner_count, tree->edge_count);

  res = pll_utree_check_integrity(tree);
  if (!res)
  {
    printf("ERROR in tree validation: %s\n", pll_errmsg);
    return PLL_FAILURE;
  }

  pll_utree_show_ascii(tree->vroot,
                       PLL_UTREE_SHOW_LABEL |
                       PLL_UTREE_SHOW_BRANCH_LENGTH |
                       PLL_UTREE_SHOW_CLV_INDEX);

  out_newick = pll_utree_export_newick(tree->vroot, NULL);
  printf("Newick export (default): %s\n", out_newick);
  free(out_newick);

  out_newick = pll_utree_export_newick(tree->vroot, newick_print_cb);
  printf("Newick export (custom): %s\n", out_newick);
  free(out_newick);

  out_newick = pll_utree_export_newick_rooted(tree->vroot, 6.13);
  printf("Newick export (rooted): %s\n", out_newick);
  free(out_newick);

  pll_utree_destroy(tree, NULL);

  return PLL_SUCCESS;
}

int test_newick_string(const char * newick)
{
  pll_utree_t * tree = pll_utree_parse_newick_string(newick);

  if (!tree && pll_errno == PLL_ERROR_TREE_INVALID)
  {
    pll_errno = 0;

    // tree must be rooted -> parse as rooted tree
    tree = pll_utree_parse_newick_string_rooted(newick);
    if (!tree || !pll_utree_is_rooted(tree))
    {
      printf("ERROR parsing newick string as ROOTED tree: %s\n", pll_errmsg);
      return PLL_FAILURE;
    }
    pll_utree_destroy(tree, NULL);

    printf("NOTE: tree was automatically unrooted!\n");
    tree = pll_utree_parse_newick_string_unroot(newick);
  }

  assert(!pll_utree_is_rooted(tree));

  if (!tree)
  {
    printf("ERROR parsing tree string: %s\n", pll_errmsg);
    exit(1);
  }

  return test_tree(tree);
}

int test_newick_file(const char * fname)
{
  pll_utree_t * tree = pll_utree_parse_newick(fname);

  if (!tree && pll_errno == PLL_ERROR_TREE_INVALID)
  {
    pll_errno = 0;
    
    // tree must be rooted -> parse as rooted tree
    tree = pll_utree_parse_newick_rooted(fname);
    if (!tree || !pll_utree_is_rooted(tree))
    {
      printf("ERROR parsing newick file as ROOTED tree: %s\n", pll_errmsg);
      return PLL_FAILURE;
    }
    pll_utree_destroy(tree, NULL);

    printf("NOTE: tree was automatically unrooted!\n");
    tree = pll_utree_parse_newick_unroot(fname);
  }

  assert(!pll_utree_is_rooted(tree));

  if (!tree)
  {
    printf("ERROR parsing tree file: %s\n", pll_errmsg);
    exit(1);
  }

  return test_tree(tree);
}


int main(int argc, char * argv[])
{
  unsigned int i;
  unsigned int attributes = get_attributes(argc, argv);

  if (attributes != PLL_ATTRIB_ARCH_CPU)
    skip_test();

  for (i = 0; i < TREE_COUNT; ++i)
  {
    printf("*** TREE # S%d\n", i+1);
    test_newick_string(newick_tree_list[i]);
    printf("\n");
  }

  for (i = 0; i < TREEFILE_COUNT; ++i)
  {
    printf("*** TREE # F%d\n", i+1);
    test_newick_file(newick_file_list[i]);
    printf("\n");
  }

  return (EXIT_SUCCESS);
}
