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
#include <locale.h>

#define NUM_BRANCHES 9
#define N_STATES_AA 20

#define FLOAT_PRECISION 4

unsigned int params_indices[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static double testbranches[NUM_BRANCHES] = {0.1, 0.2, 0.5, 0.9, 1.5, 5, 10, 50, 90};

int main(int argc, char *argv[]) {
  double d_f, dd_f;
  unsigned int n_tips = 5;

  clock_t begin = clock();

  pll_operation_t *operations;
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

  float pll_partition_create_secs = 0;
  float pll_update_prob_matrices_secs = 0;
  float pll_update_partials_secs = 0;
  float pll_set_tip_states_secs = 0;
  float pll_update_sumtable_secs = 0;
  float pll_compute_edge_loglikelihood_secs = 0;
  float pll_compute_likelihood_derivatives_secs = 0;

  /* check attributes */
  pll_common_args_t *common_args = get_common_args(argc, argv);
  unsigned int n_sites = common_args->n_sites;
  unsigned int n_categories = common_args->n_categories;

  clock_t pll_partition_create_begin_time = clock();

  pll_partition_t *partition;
  partition = pll_partition_create(
          n_tips,                 /* numer of tips */
          4,                      /* clv buffers */
          N_STATES_AA,            /* number of states */
          n_sites,                /* sequence length */
          1,                      /* different rate parameters */
          2 * n_tips - 3,         /* probability matrices */
          n_categories,           /* gamma categories */
          0,                      /* scale buffers */
          common_args->attributes /* attributes */
  );

  clock_t pll_partition_create_end_time = clock();
  pll_partition_create_secs +=
          (float) (pll_partition_create_end_time - pll_partition_create_begin_time) / CLOCKS_PER_SEC;

  if (!partition) {
    printf("Fail creating partition");
    return (-1);
  }

  double *sumtable = pll_allocate_sumtable(partition);

  if (!sumtable) {
    printf("Fail creating sumtable");
    pll_partition_destroy(partition);
    return (-1);
  }

  double branch_lengths[4] = {0.1, 0.2, 0.3, 0.4};
  unsigned int matrix_indices[4] = {0, 1, 2, 3};

  pll_set_frequencies(partition, 0, pll_aa_freqs_dayhoff);
  pll_set_subst_params(partition, 0, pll_aa_rates_dayhoff);

  const char *ref_seq = "ARNDCQEGHILKMFPSTWYV";
  size_t align_seqs = 5;

  char **align = calloc(align_seqs, sizeof(char *));

  srand(common_args->seed);

  clock_t pll_set_tip_states_before_time = clock();

  int return_val = PLL_SUCCESS;
  for (unsigned int i = 0; i < align_seqs; i++) {
    align[i] = calloc(n_sites + 1, sizeof(char));
    for (unsigned int j = 0; j < n_sites; j++) {
      align[i][j] = ref_seq[rand() % strlen(ref_seq)];
    }
    return_val &= pll_set_tip_states(partition, i, pll_map_aa,
                                     align[i]);
  }

  clock_t pll_set_tip_states_end_time = clock();
  pll_set_tip_states_secs +=
          (float) (pll_set_tip_states_end_time - pll_set_tip_states_before_time) / CLOCKS_PER_SEC;


  if (!return_val)
    fatal("Error setting tip states");

  for (unsigned int n = 0; n < common_args->n_alpha_values; ++n) {

    double *rate_cats = (double *) malloc(n_categories * sizeof(double));

    if (pll_compute_gamma_cats(common_args->alpha_values[n], n_categories, rate_cats, PLL_GAMMA_RATES_MEAN)
        == PLL_FAILURE) {
      printf("Fail computing the gamma rates\n");
      exit(0);
    }

    pll_set_category_rates(partition, rate_cats);
    free(rate_cats);

    for (unsigned int j = 0; j < partition->rate_matrices; ++j) {
      pll_update_invariant_sites_proportion(partition, j, common_args->pinvar);
    }

    clock_t pll_update_prob_matrices_begin_time = clock();

    for (int i = 0; i < common_args->n_pmatrix_itr; i++) {
      pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);
    }

    clock_t pll_update_prob_matrices_end_time = clock();
    pll_update_prob_matrices_secs +=
            (float) (pll_update_prob_matrices_end_time - pll_update_prob_matrices_begin_time) / CLOCKS_PER_SEC;

    clock_t pll_update_partials_begin_time = clock();

    pll_update_partials(partition, operations, 3);

    clock_t pll_update_partials_end_time = clock();
    pll_update_partials_secs +=
            (float) (pll_update_partials_end_time - pll_update_partials_begin_time) / CLOCKS_PER_SEC;

    clock_t edge_loglikelihood_begin_time = clock();

    pll_compute_edge_loglikelihood(partition,
                                   6,
                                   PLL_SCALE_BUFFER_NONE,
                                   7,
                                   PLL_SCALE_BUFFER_NONE,
                                   0,
                                   params_indices,
                                   NULL);

    clock_t edge_loglikelihood_end_time = clock();
    pll_compute_edge_loglikelihood_secs +=
            (float) (edge_loglikelihood_end_time - edge_loglikelihood_begin_time) / CLOCKS_PER_SEC;

    clock_t pll_update_sumtable_begin_time = clock();

    pll_update_sumtable(partition, 6, 7,
                        PLL_SCALE_BUFFER_NONE, PLL_SCALE_BUFFER_NONE,
                        params_indices, sumtable);

    clock_t pll_update_sumtable_end_time = clock();
    pll_update_sumtable_secs +=
            (float) (pll_update_sumtable_end_time - pll_update_sumtable_begin_time) / CLOCKS_PER_SEC;


    clock_t pll_compute_likelihood_derivatives_begin_time = clock();

    for (unsigned int b = 0; b < NUM_BRANCHES; ++b) {
      pll_compute_likelihood_derivatives(partition,
                                         PLL_SCALE_BUFFER_NONE,
                                         PLL_SCALE_BUFFER_NONE,
                                         testbranches[b],
                                         params_indices,
                                         sumtable,
                                         &d_f, &dd_f);
    }

    clock_t pll_compute_likelihood_derivatives_end_time = clock();
    pll_compute_likelihood_derivatives_secs +=
            (float) (pll_compute_likelihood_derivatives_end_time - pll_compute_likelihood_derivatives_begin_time) /
            CLOCKS_PER_SEC;
  }

  pll_aligned_free(sumtable);
  pll_partition_destroy(partition);

  for (unsigned int i = 0; i < align_seqs; i++) {
    free(align[i]);
  }
  free(align);

  clock_t end = clock();
  float total_secs = (float) (end - begin) / CLOCKS_PER_SEC;

  printf("Execution Mode:                     ");
  if (partition->attributes & PLL_ATTRIB_ARCH_AVX)
  {
    printf("avx");
  }
  else if (partition->attributes & PLL_ATTRIB_ARCH_SSE)
  {
    printf("sse");
  }
  else if (partition->attributes & PLL_ATTRIB_ARCH_AVX2)
  {
    printf("avx2");
  }
  else if (partition->attributes & PLL_ATTRIB_ARCH_AVX512F)
  {
    printf("avx512f");
  }
  else
  {
    printf("cpu");
  }

  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    printf(" tv");
  }
  if (partition->attributes & PLL_ATTRIB_SITE_REPEATS)
  {
    printf(" sr");
  }
  if (partition->attributes & PLL_ATTRIB_SIMD_MEM_LAYOUT)
  {
    printf(" sml");
  }
  printf("\n");
  printf("Seed:                               %d\n", common_args->seed);
  setlocale(LC_NUMERIC, "");
  printf("No. Sites:                          %'d\n", n_sites);
  printf("No. Categories:                     %d\n", n_categories);
  printf("Proportion of invariant Sites:      %f\n", common_args->pinvar);
  printf("Alpha Values:                       {");
  printf("%f", common_args->alpha_values[0]);
  for(int i = 1; i < common_args->n_alpha_values; i++) {
    printf(", %f", common_args->alpha_values[i]);
  }
  printf("}\n");
  printf("No. pll_update_prob_matrices calls: %'d\n", common_args->n_pmatrix_itr);
  printf("pll_partition_create:               %f\n", pll_partition_create_secs);
  printf("pll_set_tip_states:                 %f\n", pll_set_tip_states_secs);
  printf("pll_update_prob_matrices:           %f\n", pll_update_prob_matrices_secs);
  printf("pll_update_partials:                %f\n", pll_update_partials_secs);
  printf("pll_compute_edge_loglikelihood:     %f\n", pll_compute_edge_loglikelihood_secs);
  printf("pll_update_sumtable:                %f\n", pll_update_sumtable_secs);
  printf("pll_compute_likelihood_derivatives: %f\n", pll_compute_likelihood_derivatives_secs);
  printf("\n");
  printf("Total exec time:                    %f\n", total_secs);

  free(operations);
  destroy_common_args(&common_args);
  return (0);
}
