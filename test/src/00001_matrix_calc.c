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

#define N_STATES_AA 20
#define N_CAT_GAMMA 8
#define FLOAT_PRECISION 4

static double alpha = 0.5;
static unsigned int n_cat_gamma = N_CAT_GAMMA;
unsigned int params_indices[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

int main(int argc, char * argv[])
{
  unsigned int j;
  unsigned int n_sites = 0;
  unsigned int n_tips = 8;
  double rate_cats[N_CAT_GAMMA];

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  pll_partition_t * partition;

  partition = pll_partition_create(
          n_tips,      /* numer of tips */
          5,           /* clv buffers */
          N_STATES_AA, /* number of states */
          n_sites,     /* sequence length */
          1,           /* different rate parameters */
          2*n_tips-3,  /* probability matrices */
          n_cat_gamma, /* gamma categories */
          0,           /* scale buffers */
          attributes
  );          /* attributes */

  if (!partition)
  {
    printf("Error %d: %s\n", pll_errno, pll_errmsg);
    fatal("Fail creating partition");
  }

  double branch_lengths[8] = { 0.1, 0.2, 1, 1, 1.1, 1.2, 1.3, 1.4};
  unsigned int matrix_indices[8] = { 0, 1, 2, 3, 4, 5, 6, 7};

  if (pll_compute_gamma_cats(alpha, n_cat_gamma, rate_cats, PLL_GAMMA_RATES_MEAN) == PLL_FAILURE)
  {
    printf("Error %d: %s\n", pll_errno, pll_errmsg);
    fatal("Fail computing gamma cats");
  }

  pll_set_frequencies(partition, 0, pll_aa_freqs_dayhoff);
  pll_set_subst_params(partition, 0, pll_aa_rates_dayhoff);

  pll_set_category_rates(partition, rate_cats);

  for(int i = 0; i < 1000; i++)
  pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 8);

  for (j = 0; j < 8; ++j)
  {
    printf ("[%d] P-matrix for branch length %f\n", j+1, branch_lengths[j]);
    pll_show_pmatrix(partition, j, FLOAT_PRECISION);
    printf ("\n");
  }

  assert(n_sites == partition->sites);

  pll_partition_destroy(partition);

  return (0);
}
