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

#define NUM_BRANCHES 9
#define N_STATES_AA 20
#define N_CATEGORIES 4

#define FLOAT_PRECISION 4

unsigned int params_indices[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static double testbranches[NUM_BRANCHES] = {0.1, 0.2, 0.5, 0.9, 1.5, 5, 10, 50, 90};

int main(int argc, char *argv[]) {
  unsigned int j, b;
  double d_f, dd_f;
  unsigned int n_sites = 10000;
  unsigned int n_tips = 5;
  pll_operation_t *operations;
  double *sumtable;

  clock_t begin = clock();

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

  pll_partition_t *partition;
  partition = pll_partition_create(
          n_tips,      /* numer of tips */
          4,           /* clv buffers */
          N_STATES_AA, /* number of states */
          n_sites,     /* sequence length */
          1,           /* different rate parameters */
          2 * n_tips - 3,  /* probability matrices */
          4,           /* gamma categories */
          0,           /* scale buffers */
          attributes
  );          /* attributes */

  if (!partition) {
    printf("Fail creating partition");
    return (-1);
  }

  int elems_per_reg = partition->alignment / sizeof(double);
  int sites_padded = (partition->sites + elems_per_reg - 1) & (0xFFFFFFFF - elems_per_reg + 1);

  sumtable = pll_aligned_alloc(
          sites_padded * partition->rate_cats * partition->states *
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

  const char *ref_seq = "PIGLRVTLRRDRMWI";
  size_t align_seqs = 5;

  char **align = calloc(align_seqs, sizeof(char *));

  srand(time(NULL));

  int return_val = PLL_SUCCESS;
  for (size_t i = 0; i < align_seqs; i++) {
    align[i] = calloc(n_sites + 1, sizeof(char));
    for (size_t k = 0; k < n_sites; k++) {
      align[i][k] = ref_seq[rand() % strlen(ref_seq)];
    }
    return_val &= pll_set_tip_states(partition, i, pll_map_aa,
                                     align[i]);
  }

  if (!return_val)
    fatal("Error setting tip states");

  double *rate_cats = (double *) malloc(N_CATEGORIES * sizeof(double));

  if (pll_compute_gamma_cats(0.1, N_CATEGORIES, rate_cats, PLL_GAMMA_RATES_MEAN) == PLL_FAILURE) {
    printf("Fail computing the gamma rates\n");
    exit(0);
  }

  pll_set_category_rates(partition, rate_cats);
  free(rate_cats);

  for (j = 0; j < partition->rate_matrices; ++j) {
    pll_update_invariant_sites_proportion(partition, j, 0.0);
  }

  clock_t partials_time = clock();
  float init_secs = (float) (partials_time - begin) / CLOCKS_PER_SEC;
  printf("Init time:                     %f\n", init_secs);

  pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);
  pll_update_partials(partition, operations, 3);

  clock_t edge_loglikelihood_time = clock();
  float partials_secs = (float) (edge_loglikelihood_time - partials_time) / CLOCKS_PER_SEC;
  printf("Prob matrices & partials time: %f\n", partials_secs);

  pll_compute_edge_loglikelihood(partition,
                                 6,
                                 PLL_SCALE_BUFFER_NONE,
                                 7,
                                 PLL_SCALE_BUFFER_NONE,
                                 0,
                                 params_indices,
                                 NULL);

  clock_t sumtable_time = clock();
  float edge_loglikelihood_secs = (float) (sumtable_time - edge_loglikelihood_time) / CLOCKS_PER_SEC;
  printf("Edge loglikelihood time:       %f\n", edge_loglikelihood_secs);

  pll_update_sumtable(partition, 6, 7,
                      PLL_SCALE_BUFFER_NONE, PLL_SCALE_BUFFER_NONE,
                      params_indices, sumtable);

  clock_t derivetive_time = clock();
  float sumtable_secs = (float) (derivetive_time - sumtable_time) / CLOCKS_PER_SEC;
  printf("Sumtable time:                 %f\n", sumtable_secs);

  for (b = 0; b < NUM_BRANCHES; ++b) {
    pll_compute_likelihood_derivatives(partition,
                                       PLL_SCALE_BUFFER_NONE,
                                       PLL_SCALE_BUFFER_NONE,
                                       testbranches[b],
                                       params_indices,
                                       sumtable,
                                       &d_f, &dd_f);
  }

  clock_t cleanup_time = clock();
  float derivetive_secs = (float) (cleanup_time - derivetive_time) / CLOCKS_PER_SEC;
  printf("Derivetive time:               %f\n", derivetive_secs);

  pll_aligned_free(sumtable);
  pll_partition_destroy(partition);

  for (size_t i = 0; i < align_seqs; i++) {
    free(align[i]);
  }
  free(align);

  clock_t end = clock();
  float total_secs = (float) (end - begin) / CLOCKS_PER_SEC;

  printf("\n");
  printf("Total exec time:               %f\n", total_secs);

  free(operations);
  return (0);
}
