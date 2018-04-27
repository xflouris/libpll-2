#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <stdarg.h>
#include <time.h>

const pll_state_t odd5_map[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1f, 0, 0, 0x1f, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1f, 0, 0x01, 0x02, 0x04,
    0x08, 0x0c, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static unsigned int get_attribute(char *arg)
{
  if (!strcmp (arg, "tv"))
  {
    /* tipvector */
    return PLL_ATTRIB_PATTERN_TIP;
  }
  else if (!strcmp (arg, "sr"))
  {
    /* site repeats */
    return PLL_ATTRIB_SITE_REPEATS;
  }
  else if (!strcmp (arg, "sml"))
  {
    /* SIMD memory layout */
    return PLL_ATTRIB_SIMD_MEM_LAYOUT;
  }
  else if (!strcmp (arg, "avx"))
  {
    /* avx vectorization */
    return PLL_ATTRIB_ARCH_AVX;
  }
  else if (!strcmp (arg, "sse"))
  {
    /* sse3 vectorization */
    return PLL_ATTRIB_ARCH_SSE;
  }
  else if (!strcmp (arg, "avx2"))
  {
    /* avx2 vectorization */
    return PLL_ATTRIB_ARCH_AVX2;
  }
  else if (!strcmp (arg, "avx512f"))
  {
    /* avx512f vectorization */
    return PLL_ATTRIB_ARCH_AVX512F;
  }
  else
  {
    printf("Unrecognised attribute: %s\n", arg);
    exit(1);
  }
}

unsigned int get_attributes(int argc, char **argv)
{
  unsigned int attributes = PLL_ATTRIB_ARCH_CPU;

  for (int i=1; i<argc; ++i)
  {
    attributes |= get_attribute(argv[i]);
  }
  return attributes;
}

static unsigned int count_values(const char* str) {
  unsigned int count = 0;
  char* values_copy = xstrdup(str);
  char* token = strtok(values_copy, ",");
  while (token != NULL) {
    count++;
    token = strtok(NULL, ",");
  }
  return count;
}

pll_common_args_t* get_common_args(int argc, char **argv) {
  pll_common_args_t* common_args = xmalloc(sizeof(pll_common_args_t));
  common_args->attributes = PLL_ATTRIB_ARCH_CPU;
  common_args->alpha_values = xmalloc(sizeof(double));
  common_args->n_alpha_values = 1;
  common_args->alpha_values[0] = 0.1;
  common_args->n_pmatrix_itr = 1;
  common_args->n_categories = 4;
  common_args->seed = (unsigned int) time(NULL);
  common_args->print_seq = 0;
  common_args->n_sites = 1000000;
  common_args->n_states = 20;
  common_args->pinvar = 0.0;

  for (int i=1; i<argc; ++i)
  {
    if (strstr (argv[i], "-alpha=") != NULL)
    {
      char* values = strstr (argv[i], "=") + 1;
      if(*values == '\0' || *values == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }

      common_args->n_alpha_values = count_values(values);
      common_args->alpha_values = xmalloc(sizeof(double)*common_args->n_alpha_values);

      unsigned int idx = 0;
      char* values_copy = xstrdup(values);
      char* token = strtok(values_copy, ",");
      while (token != NULL) {
        common_args->alpha_values[idx] = strtod(token, NULL);
        if(common_args->alpha_values[idx] <= 0) {
          printf("alpha value must not be smaller 0: %s\n", token);
          exit(1);
        }
        idx++;
        token = strtok(NULL, ",");
      }
    }
    else if (strstr (argv[i], "-n-states=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->n_states = (unsigned int) strtol(value, NULL, 10);
      if(common_args->n_states != 4 && common_args->n_states != 20) {
        printf("Only specific number of states are allowed (4 / DNA or 20 / AA): %s\n", argv[i]);
        exit(1);
      }
    }
    else if (strstr (argv[i], "-n-categories=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->n_categories = (unsigned int) strtol(value, NULL, 10);
      if(common_args->n_categories == 0) {
        printf("Number of categories must not be 0: %s\n", argv[i]);
        exit(1);
      }
    }
    else if (strstr (argv[i], "-seed=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->seed = (unsigned int) strtol(value, NULL, 10);
    }
    else if (strstr (argv[i], "-n-sites=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->n_sites = (unsigned int) strtol(value, NULL, 10);
      if(common_args->n_sites == 0) {
        printf("Number of sites must not be 0: %s\n", argv[i]);
        exit(1);
      }
    }
    else if (strstr (argv[i], "-print-seq") != NULL)
    {
      common_args->print_seq = 1;
    }
    else if (strstr (argv[i], "-n-pmatrix-itr=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->n_pmatrix_itr = (unsigned int) strtol(value, NULL, 10);
      if(common_args->n_pmatrix_itr == 0) {
        printf("Number of sites must not be 0: %s\n", argv[i]);
        exit(1);
      }
    }
    else if (strstr (argv[i], "-p-invar=") != NULL)
    {
      char* value = strstr (argv[i], "=") + 1;
      if(*value == '\0' || *value == ',') {
        printf("Unable to read value from: %s\n", argv[i]);
        exit(1);
      }
      common_args->pinvar = strtod(value, NULL);
    }
    else
    {
      common_args->attributes |= get_attribute(argv[i]);
    }
  }

  return common_args;
}

void destroy_common_args(pll_common_args_t** args) {
  if(args == NULL) {
    return;
  }
  if((*args)->alpha_values != NULL) {
    free((*args)->alpha_values);
  }
  free(*args);
  *args = NULL;
}

void skip_test ()
{
  printf ("Skip\n");
  exit (0);
}

/* note: There is no exhaustive error checking in this function.
         It is intended to use with the test datasets that were
         validated in advance. */
pll_partition_t * parse_msa(const char * filename,
                           unsigned int states,
                           unsigned int rate_cats,
                           unsigned int rate_matrices,
                           pll_utree_t * tree,
                           unsigned int attributes)
{
  return parse_msa_reduced(filename,
                          states,
                          rate_cats,
                          rate_matrices,
                          tree,
                          attributes,
                          -1);
}

pll_partition_t * parse_msa_reduced(const char * filename,
                            unsigned int states,
                            unsigned int rate_cats,
                            unsigned int rate_matrices,
                            pll_utree_t * tree,
                            unsigned int attributes,
                            unsigned int max_sites)
{
  unsigned int i;
  unsigned int taxa_count = tree->tip_count;
  pll_partition_t * partition;
  long hdrlen, seqlen, seqno;
  char * seq = NULL,
       * hdr = NULL;

  /* open FASTA file */
  pll_fasta_t * fp = pll_fasta_open(filename, pll_map_fasta);
  if (!fp)
  {
    printf("Error opening file %s", filename);
    return NULL;
  }

  /* allocate arrays to store FASTA headers and sequences */
  char ** headers = (char **)calloc(taxa_count, sizeof(char *));
  char ** seqdata = (char **)calloc(taxa_count, sizeof(char *));

  /* read FASTA sequences and make sure they are all of the same length */
  int sites = -1;
  for (i = 0; pll_fasta_getnext(fp,&hdr,&hdrlen,&seq,&seqlen,&seqno); ++i)
  {
    if (sites == -1) sites = seqlen;

    headers[i] = hdr;
    seqdata[i] = seq;
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
  {
    printf("Error while reading file %s", filename);
    return NULL;
  }

  /* close FASTA file */
  pll_fasta_close(fp);

  if (sites == -1)
  {
    printf("Unable to read alignment");
    return NULL;
  }

  if (max_sites != -1)
    sites = max_sites;

  partition = pll_partition_create(taxa_count,           /* tip nodes */
                                   taxa_count - 2,       /* inner nodes */
                                   states,
                                   (unsigned int)sites,
                                   rate_matrices,        /* rate matrices */
                                   2*taxa_count - 3,     /* prob matrices */
                                   rate_cats,            /* rate categories */
                                   taxa_count - 2,       /* scale buffers */
                                   attributes);

  /* create a libc hash table of size tip_nodes_count */
  hcreate(taxa_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *)malloc(taxa_count *
                                               sizeof(unsigned int));
  for (i = 0; i < taxa_count; ++i)
  {
    data[i] = tree->nodes[i]->clv_index;
    ENTRY entry;
#ifdef __APPLE__
    entry.key = xstrdup(tree->nodes[i]->label);
#else
    entry.key = tree->nodes[i]->label;
#endif
    entry.data = (void *)(data+i);

    hsearch(entry, ENTER);
  }

  for (i = 0; i < taxa_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
    {
      printf("Sequence with header %s does not appear in the tree", headers[i]);
      return NULL;
    }

    unsigned int tip_clv_index = *((unsigned int *)(found->data));

    const pll_state_t * map = states == 4? pll_map_nt : pll_map_aa;
    pll_set_tip_states(partition, tip_clv_index, map, seqdata[i]);

    free(headers[i]);
    free(seqdata[i]);
  }

  /* clean */
  hdestroy();
  free(data);
  free(headers);
  free(seqdata);

  return partition;
}

int cb_full_traversal(pll_unode_t * node)
{
  return 1;
}

int cb_rfull_traversal(pll_rnode_t * node)
{
  return 1;
}

__attribute__((format(printf, 1, 2)))
void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

void * xmalloc(size_t size)
{
  void * t;
  t = malloc(size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}

char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
}  
