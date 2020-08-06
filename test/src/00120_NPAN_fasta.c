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
#define N_STATES    20
#define N_RATE_CATS  4

#define N_TAXA   21
#define N_SITES 113

#define FASTA_FILE "testdata/ribosomal_l5_pf00673.fas"


static int load_fasta_and_set_tips(const char * fasta_fname, pll_partition_t * partition,
                                   const pll_state_t * map)
{
  unsigned int i;
  int retval = PLL_FAILURE;

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
    if (!pll_set_tip_states (partition, i, map, msa->sequence[i]))
    {
      printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
      retval = i+1;
      break;
    }
  }

  pll_msa_destroy(msa);

  return retval;
}

static int failtest(unsigned int attributes, pll_bool_t oneliner)
{
  unsigned int i;
  pll_partition_t * partition;
  int retval = PLL_FAILURE;

  partition = pll_partition_create(N_TAXA, /* tips */
                                    4, /* clv buffers */
                                    N_STATES, /* states */
                                    N_SITES, /* sites */
                                    1, /* different rate parameters */
                                    8, /* probability matrices */
                                    N_RATE_CATS, /* rate categories */
                                    1, 
                                    attributes
                                    );

  if (oneliner)
  {
    retval = load_fasta_and_set_tips(FASTA_FILE, partition, pll_map_nt);
  }
  else
  {
    char * seq, *header;
    long seq_len, header_len, seqno;
    pll_fasta_t * fp = pll_fasta_open (FASTA_FILE, pll_map_fasta);

    i = 0;
    while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
    {
      if (!pll_set_tip_states (partition, i, pll_map_nt, seq))
      {
        free (header);
        free (seq);
        retval = i+1;
        printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
        break;
      }
      free (header);
      free (seq);
      ++i;
    }
    pll_fasta_close (fp);
  }

  pll_partition_destroy(partition);

  return retval;
}

static int proteintest(unsigned int attributes, pll_bool_t oneliner)
{
  unsigned int i;
  pll_partition_t * partition;

  printf ("Creating PLL partition\n");

  partition = pll_partition_create(N_TAXA, /* tips */
                                    4, /* clv buffers */
                                    N_STATES, /* states */
                                    N_SITES, /* sites */
                                    1, /* different rate parameters */
                                    8, /* probability matrices */
                                    N_RATE_CATS, /* rate categories */
                                    1,
                                    attributes);

  if (oneliner)
  {
    int retval = load_fasta_and_set_tips(FASTA_FILE, partition, pll_map_aa);
    assert(!retval);
  }
  else
  {
    char * seq, *header;
    long seq_len, header_len, seqno;

    pll_fasta_t * fp = pll_fasta_open (FASTA_FILE, pll_map_fasta);
    if (!fp)
    {
      printf (" ERROR opening file (%d): %s\n", pll_errno, pll_errmsg);
      return (PLL_FAILURE);
    }

    i = 0;
    while (pll_fasta_getnext (fp, &header, &header_len, &seq, &seq_len, &seqno))
    {
      if (seq_len != N_SITES)
      {
        printf (
            " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
            i, seq_len, N_SITES);
        return (PLL_FAILURE);
      }
      if (!pll_set_tip_states (partition, i, pll_map_aa, seq))
      {
        printf (" ERROR setting states (%d): %s\n", pll_errno, pll_errmsg);
        return (PLL_FAILURE);
      }
      printf ("Header of sequence %d(%ld) %s (%ld sites)\n", i, seqno, header,
              seq_len);
      printf ("   %s\n", seq);
      free (header);
      free (seq);
      ++i;
    }

    if (pll_errno != PLL_ERROR_FILE_EOF)
    {
      printf (" ERROR at the end (%d): %s\n", pll_errno, pll_errmsg);
      return (PLL_FAILURE);
    }

    if (i != N_TAXA)
    {
      printf (" ERROR: Number of taxa mismatch (%d): %d\n", i, N_TAXA);
      return (PLL_FAILURE);
    }

    pll_fasta_close (fp);
  }

  pll_partition_destroy(partition);

  return PLL_SUCCESS;
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);
  int fail_retval;

  if (proteintest (attributes, PLL_FALSE))
    printf ("Test (low-level): OK\n\n");

  if (proteintest (attributes, PLL_TRUE))
    printf ("Test (one-liner): OK\n\n");

  fail_retval = failtest (attributes, PLL_FALSE);
  if (fail_retval)
    printf ("Fail test (low-level): OK (sequence %d)\n\n", fail_retval);

  fail_retval = failtest (attributes, PLL_TRUE);
  if (fail_retval)
    printf ("Fail test (one-liner): OK (sequence %d)\n\n", fail_retval);

  return PLL_SUCCESS;
}
