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

/*
    derivatives-aa-benchmark.c

    Benchmark version of derivatives-aa.c
 */
#include "common.h"
#include <time.h>

#define NUM_ALPHAS   3
#define NUM_BRANCHES 9
#define NUM_CATS     3
#define NUM_PINV     4
#define N_STATES_AA 20

#define FLOAT_PRECISION 4

static double alpha[NUM_ALPHAS] = {0.1, 0.75, 1.5};
static double pinvar[NUM_PINV] = {0.0, 0.3, 0.5, 0.6};
static unsigned int n_cat_gamma[NUM_CATS] = {1, 2, 4};
unsigned int params_indices[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int main(int argc, char *argv[]) {
  unsigned int i, j, k, p;
  unsigned int n_sites = 1000;
  unsigned int n_tips = 5;
  pll_operation_t *operations;
  double *sumtable;

  operations = (pll_operation_t *) malloc(4 * sizeof(pll_operation_t));

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

  /* additional operation for moving the root into the tip branch */

  operations[3].parent_clv_index = 7;
  operations[3].child1_clv_index = 6;
  operations[3].child2_clv_index = 3;
  operations[3].child1_matrix_index = 0;
  operations[3].child2_matrix_index = 0;
  operations[3].parent_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[3].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
  operations[3].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  for (k = 0; k < NUM_CATS; ++k) {
    pll_partition_t *partition;
    partition = pll_partition_create(
            n_tips,      /* numer of tips */
            4,           /* clv buffers */
            N_STATES_AA, /* number of states */
            n_sites,     /* sequence length */
            1,           /* different rate parameters */
            2 * n_tips - 3,  /* probability matrices */
            n_cat_gamma[k], /* gamma categories */
            0,           /* scale buffers */
            attributes
    );          /* attributes */

    if (!partition) {
      printf("Fail creating partition");
      return (-1);
    }

    sumtable = pll_aligned_alloc(
            partition->sites * partition->rate_cats * partition->states_padded *
            sizeof(double), partition->alignment);

    if (!sumtable) {
      printf("Fail creating sumtable");
      pll_partition_destroy(partition);
      return (-1);
    }

    double branch_lengths[4] = {0.1, 0.2, 0.3, 0.4};
    unsigned int matrix_indices[4] = {0, 1, 2, 3};

    pll_set_frequencies(partition, 0, pll_aa_freqs_dayhoff);
    pll_set_subst_params(partition, 0, pll_aa_rates_dayhoff);

    int return_val;

    const char *ref_seq = "PIGLRVTLRRDRMWI";

    size_t align_size = n_sites * strlen(ref_seq);
    size_t align_seqs = 5;

    char **align = calloc(align_seqs, sizeof(char *));

    srand(time(NULL));

    return_val = PLL_SUCCESS;
    for (size_t i = 0; i < align_seqs; i++) {
      align[i] = calloc(align_size + 1, sizeof(char));
      for (size_t k = 0; k < align_size; k++) {
        align[i][k] = ref_seq[rand() % strlen(ref_seq)];
      }
      return_val &= pll_set_tip_states(partition, i, pll_map_aa,
                                       align[i]);
    }

    if (!return_val)
      fatal("Error setting tip states");

    for (i = 0; i < NUM_ALPHAS; ++i) {
      for (p = 0; p < NUM_PINV; ++p) {

        double *rate_cats = (double *) malloc(n_cat_gamma[k] * sizeof(double));

        if (pll_compute_gamma_cats(alpha[i], n_cat_gamma[k], rate_cats, PLL_GAMMA_RATES_MEAN) == PLL_FAILURE) {
          printf("Fail computing the gamma rates\n");
          continue;
        }

        pll_set_category_rates(partition, rate_cats);
        free(rate_cats);

        for (j = 0; j < partition->rate_matrices; ++j) {
          pll_update_invariant_sites_proportion(partition, j, pinvar[p]);
        }

        pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);
        pll_update_partials(partition, operations, 3);


        for (int i = 0; i < 1000; i++)
          pll_compute_edge_loglikelihood(partition,
                                         6,
                                         PLL_SCALE_BUFFER_NONE,
                                         7,
                                         PLL_SCALE_BUFFER_NONE,
                                         0,
                                         params_indices,
                                         NULL);
      }
    }

    pll_aligned_free(sumtable);
    pll_partition_destroy(partition);

    for (size_t i = 0; i < align_seqs; i++) {
      free(align[i]);
    }
    free(align);
  }


  free(operations);
  return (0);
}
