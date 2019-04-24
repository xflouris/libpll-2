/*
 Copyright (C) 2015 Diego Darriba, Tomas Flouri

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#include "common.h"
#include <string.h>
#include <assert.h>

#define ALPHA        0.5
#define N_STATES     4
#define N_RATE_CATS  4

#define N_TAXA_BIG    43
#define N_TAXA_SMALL   5
#define N_SITES_BIG   492
#define N_SITES_SMALL 1964

#define SMALL_FASTA "testdata/small.fas"
#define BIG_FASTA "testdata/worms16s.fas"
#define BAD_FASTA "unexistent-file"

static void load_fasta_and_set_tips(const char * fasta_fname, pll_partition_t * partition)
{
  unsigned int i;

  pll_msa_t * msa = pll_fasta_load(fasta_fname);
  if (!msa)
  {
    printf (" ERROR loading MSA from FASTA file (%d): %s\n", pll_errno, pll_errmsg);
    exit (PLL_FAILURE);
  }

  assert(msa->count == partition->tips);
  assert(msa->length == partition->sites);

  for (i = 0; i < msa->count; ++i)
  {
    if (!pll_set_tip_states (partition, i, pll_map_nt, msa->sequence[i]))
    {
      printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
      exit (PLL_FAILURE);
    }
  }

  pll_msa_destroy(msa);
}


static int failtest(unsigned int attributes, pll_bool_t oneliner)
{
  if (oneliner)
  {
    pll_msa_t * msa = pll_fasta_load(BAD_FASTA);
    assert(!msa && pll_errno == PLL_ERROR_FILE_OPEN);
  }
  else
  {
    pll_fasta_t * fp;
    fp = pll_fasta_open (BAD_FASTA, pll_map_fasta);
    assert(!fp && pll_errno == PLL_ERROR_FILE_OPEN);
  }

  return PLL_SUCCESS;
}

static int bigtest(unsigned int attributes, pll_bool_t oneliner)
{
  unsigned int i;
  pll_partition_t * partition;

  printf ("Creating PLL partition\n");

  partition = pll_partition_create(N_TAXA_BIG, /* tips */
                                    4, /* clv buffers */
                                    N_STATES, /* states */
                                    N_SITES_BIG, /* sites */
                                    1, /* different rate parameters */
                                    8, /* probability matrices */
                                    N_RATE_CATS, /* rate categories */
                                    1, 
                                    attributes
                                    );

  if (oneliner)
  {
    load_fasta_and_set_tips(BIG_FASTA, partition);
  }
  else
  {
    char * seq, *header;
    long seq_len, header_len, seqno;
    pll_fasta_t * fp = pll_fasta_open (BIG_FASTA, pll_map_fasta);
    if (!fp)
    {
      printf (" ERROR opening file (%d): %s\n", pll_errno, pll_errmsg);
      exit (PLL_FAILURE);
    }

    i = 0;
    while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
    {
      if (seq_len != N_SITES_BIG)
      {
        printf (
            " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
            i, seq_len, N_SITES_BIG);
        exit (PLL_FAILURE);
      }
      if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
      {
        printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
        exit (PLL_FAILURE);
      }
      printf ("Header of sequence %d(%ld) %s (%ld sites)\n", i, seqno, header,
              seq_len);
      free (header);
      free (seq);
      ++i;
    }

    if (pll_errno != PLL_ERROR_FILE_EOF)
    {
      printf (" ERROR at the end (%d): %s\n", pll_errno, pll_errmsg);
      exit (PLL_FAILURE);
    }

    pll_fasta_close (fp);

    assert(i == N_TAXA_BIG);
  }

  pll_partition_destroy(partition);

  return PLL_SUCCESS;
}

static int smalltest (unsigned int attributes, pll_bool_t oneliner)
{
  unsigned int i;
  pll_partition_t * partition;
  pll_operation_t * operations;
  double rate_cats[4];
  unsigned int params_indices[N_RATE_CATS] = {0,0,0,0};

  double branch_lengths[4] =
    { 0.1, 0.2, 1, 1 };
  double frequencies[4] =
    { 0.1, 0.2, 0.3, 0.4 };
  unsigned int matrix_indices[4] = { 0, 1, 2, 3 };
  double subst_params[6] = { 1, 5, 1, 1, 5, 1 };

  partition = pll_partition_create(N_TAXA_SMALL, 4,
                                   N_STATES,
                                   N_SITES_SMALL,
                                   1, 2 * N_TAXA_SMALL - 3,
                                   N_RATE_CATS,
                                   1, 
                                   attributes);


  if (oneliner)
  {
    load_fasta_and_set_tips(SMALL_FASTA, partition);
  }
  else
  {
    char * seq, *header;
    long seq_len, header_len, seqno;
    pll_fasta_t * fp = pll_fasta_open (SMALL_FASTA, pll_map_fasta);
    i = 0;
    while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
    {
      if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
        exit (PLL_FAILURE);
      free (header);
      free (seq);
      ++i;
    }
    pll_fasta_close (fp);
    assert(i == (N_TAXA_SMALL));
  }

  operations = (pll_operation_t *) malloc (4 * sizeof(pll_operation_t));

  operations[0].parent_clv_index = 5;
  operations[0].child1_clv_index = 0;
  operations[0].child2_clv_index = 1;
  operations[0].child1_matrix_index = 1;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index = 6;
  operations[1].child1_clv_index = 5;
  operations[1].child2_clv_index = 2;
  operations[1].child1_matrix_index = 0;
  operations[1].child2_matrix_index = 1;
  operations[1].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  operations[2].parent_clv_index = 7;
  operations[2].child1_clv_index = 3;
  operations[2].child2_clv_index = 4;
  operations[2].child1_matrix_index = 1;
  operations[2].child2_matrix_index = 1;
  operations[2].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[2].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  pll_compute_gamma_cats (ALPHA, N_RATE_CATS, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_subst_params (partition, 0, subst_params);
  pll_set_frequencies (partition, 0, frequencies);
  pll_set_category_rates (partition, rate_cats);
  pll_update_prob_matrices (partition, params_indices, matrix_indices, branch_lengths, 4);
  pll_update_partials (partition, operations, 3);

  printf ("logL: %17.6f\n", 
  pll_compute_edge_loglikelihood(partition,
                                 6,
                                 PLL_SCALE_BUFFER_NONE,
                                 7,
                                 PLL_SCALE_BUFFER_NONE,
                                 0,
                                 params_indices,
                                 NULL));

  free (operations);
  pll_partition_destroy(partition);

  return PLL_SUCCESS;
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);

  if (bigtest (attributes, PLL_FALSE))
    printf ("Big test (low-level): OK\n\n");

  if (bigtest (attributes, PLL_TRUE))
    printf ("Big test (one-liner): OK\n\n");

  if (smalltest (attributes, PLL_FALSE))
    printf ("Small test (low-level): OK\n\n");

  if (smalltest (attributes, PLL_TRUE))
    printf ("Small test (one-liner): OK\n\n");


  if (failtest (attributes, PLL_FALSE))
    printf ("Fail test (low-level): OK\n");

  if (failtest (attributes, PLL_TRUE))
    printf ("Fail test (one-liner): OK\n");

  return PLL_SUCCESS;
}
