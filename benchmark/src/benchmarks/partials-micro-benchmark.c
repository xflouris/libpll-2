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
    This minimal benchmark executes a fixed set of partial kernels with different memory layouts and/or instruction sets
 */

#include "common.h"
#include <time.h>
#include <locale.h>

#include <xmmintrin.h>
#include <pmmintrin.h>

#define NUM_BRANCHES 9
#define ALIGN_SEQS 5

#define FLOAT_PRECISION 4

unsigned int params_indices[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

typedef void (*kernel_ptr)(unsigned int, unsigned int, double *, unsigned int *,
                           const double *, const double *, const double *, const double *,
                           const unsigned int *, const unsigned int *, unsigned int);

void pll_core_update_partial_ii_20x20_avx2(unsigned int sites,
                                           unsigned int rate_cats,
                                           double *parent_clv,
                                           unsigned int *parent_scaler,
                                           const double *left_clv,
                                           const double *right_clv,
                                           const double *left_matrix,
                                           const double *right_matrix,
                                           const unsigned int *left_scaler,
                                           const unsigned int *right_scaler,
                                           unsigned int attrib);


void pll_core_update_partial_ii_20x20_avx2_sml(unsigned int sites,
                                               unsigned int rate_cats,
                                               double *parent_clv,
                                               unsigned int *parent_scaler,
                                               const double *left_clv,
                                               const double *right_clv,
                                               const double *left_matrix,
                                               const double *right_matrix,
                                               const unsigned int *left_scaler,
                                               const unsigned int *right_scaler,
                                               unsigned int attrib);

void pll_core_update_partial_ii_4x4_avx(unsigned int sites,
                                        unsigned int rate_cats,
                                        double *parent_clv,
                                        unsigned int *parent_scaler,
                                        const double *left_clv,
                                        const double *right_clv,
                                        const double *left_matrix,
                                        const double *right_matrix,
                                        const unsigned int *left_scaler,
                                        const unsigned int *right_scaler,
                                        unsigned int attrib);

void pll_core_update_partial_ii_20x20_avx512f_sml(unsigned int sites,
                                                  unsigned int rate_cats,
                                                  double *parent_clv,
                                                  unsigned int *parent_scaler,
                                                  const double *left_clv,
                                                  const double *right_clv,
                                                  const double *left_matrix,
                                                  const double *right_matrix,
                                                  const unsigned int *left_scaler,
                                                  const unsigned int *right_scaler,
                                                  unsigned int attrib);

void pll_core_update_partial_ii_20x20_avx512f(unsigned int sites,
                                              unsigned int rate_cats,
                                              double *parent_clv,
                                              unsigned int *parent_scaler,
                                              const double *left_clv,
                                              const double *right_clv,
                                              const double *left_matrix,
                                              const double *right_matrix,
                                              const unsigned int *left_scaler,
                                              const unsigned int *right_scaler,
                                              unsigned int attrib);

void pll_core_update_partial_ii_4x4_avx2_sml(unsigned int sites,
                                             unsigned int rate_cats,
                                             double *parent_clv,
                                             unsigned int *parent_scaler,
                                             const double *left_clv,
                                             const double *right_clv,
                                             const double *left_matrix,
                                             const double *right_matrix,
                                             const unsigned int *left_scaler,
                                             const unsigned int *right_scaler,
                                             unsigned int attrib);

void pll_core_update_partial_ii_4x4_avx512f_sml(unsigned int sites,
                                                unsigned int rate_cats,
                                                double *parent_clv,
                                                unsigned int *parent_scaler,
                                                const double *left_clv,
                                                const double *right_clv,
                                                const double *left_matrix,
                                                const double *right_matrix,
                                                const unsigned int *left_scaler,
                                                const unsigned int *right_scaler,
                                                unsigned int attrib);

#ifdef SIMD_COUNTING_MODE
static void print_and_reset_instr_counts() {
  printf("   _mm256_setzero_pd_counter:      %'lu\n", _mm256_setzero_pd_counter);
  printf("   _mm256_store_pd_counter:        %'lu\n", _mm256_store_pd_counter);
  printf("   _mm256_fmadd_pd_counter:        %'lu\n", _mm256_fmadd_pd_counter);
  printf("   _mm256_set1_pd_counter:         %'lu\n", _mm256_set1_pd_counter);
  printf("   _mm256_load_pd_counter:         %'lu\n", _mm256_load_pd_counter);
  printf("   _mm256_unpackhi_pd_counter:     %'lu\n", _mm256_unpackhi_pd_counter);
  printf("   _mm256_unpacklo_pd_counter:     %'lu\n", _mm256_unpacklo_pd_counter);
  printf("   _mm256_add_pd_counter:          %'lu\n", _mm256_add_pd_counter);
  printf("   _mm256_mul_pd_counter:          %'lu\n", _mm256_mul_pd_counter);
  printf("   _mm256_cmp_pd_counter:          %'lu\n", _mm256_cmp_pd_counter);
  printf("   _mm256_permute2f128_pd_counter: %'lu\n", _mm256_permute2f128_pd_counter);
  printf("   _mm256_movemask_pd_counter:     %'lu\n", _mm256_movemask_pd_counter);
  printf("   _mm256_blend_pd_counter:        %'lu\n", _mm256_blend_pd_counter);

  pll_reset_simd_counter();
}
#endif /* SIMD_COUNTING_MODE */

static void count_instrs_of_kernel(unsigned int n_tips,
                                    char **align,
                                    const double *pll_freqs,
                                    const double *pll_rates,
                                    const pll_state_t *pll_map,
                                    const unsigned int n_states,
                                    const pll_common_args_t *common_args,
                                    pll_operation_t *operations,
                                    kernel_ptr kernel_fun,
                                    unsigned int attributes) {
  unsigned int n_sites = common_args->n_sites;
  unsigned int n_categories = common_args->n_categories;

  pll_partition_t *partition;
  partition = pll_partition_create(
          n_tips,                 /* numer of tips */
          4,                      /* clv buffers */
          n_states,               /* number of states */
          n_sites,                /* sequence length */
          1,                      /* different rate parameters */
          2 * n_tips - 3,         /* probability matrices */
          n_categories,           /* gamma categories */
          0,                      /* scale buffers */
          attributes /* attributes */
  );

  if (!partition) {
    fatal("Fail creating partition");
  }

  double branch_lengths[4] = {0.1, 0.2, 0.3, 0.4};
  unsigned int matrix_indices[4] = {0, 1, 2, 3};

  pll_set_frequencies(partition, 0, pll_freqs);
  pll_set_subst_params(partition, 0, pll_rates);

  int return_val = PLL_SUCCESS;
  for (unsigned int i = 0; i < ALIGN_SEQS; i++) {
    return_val &= pll_set_tip_states(partition, i, pll_map,
                                     align[i]);
  }

  if (!return_val)
    fatal("Error setting tip states");

  double *rate_cats = (double *) malloc(n_categories * sizeof(double));

  if (pll_compute_gamma_cats(0.1, n_categories, rate_cats, PLL_GAMMA_RATES_MEAN)
      == PLL_FAILURE) {
    printf("Fail computing the gamma rates\n");
    exit(0);
  }

  pll_set_category_rates(partition, rate_cats);
  free(rate_cats);

  for (unsigned int j = 0; j < partition->rate_matrices; ++j) {
    pll_update_invariant_sites_proportion(partition, j, common_args->pinvar);
  }

  pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);

  pll_reset_simd_counter();

  const pll_operation_t *op = &(operations[0]);
  kernel_fun(partition->sites,
             partition->rate_cats,
             partition->clv[op->parent_clv_index],
             NULL,
             partition->clv[op->child1_clv_index],
             partition->clv[op->child2_clv_index],
             partition->pmatrix[op->child1_matrix_index],
             partition->pmatrix[op->child2_matrix_index],
             NULL,
             NULL,
             partition->attributes);

  print_and_reset_instr_counts();

  pll_partition_destroy(partition);
}

static float benchmark_kernel(unsigned int n_tips,
                              char **align,
                              const double *pll_freqs,
                              const double *pll_rates,
                              const pll_state_t *pll_map,
                              const unsigned int n_states,
                              const pll_common_args_t *common_args,
                              pll_operation_t *operations,
                              kernel_ptr kernel_fun,
                              unsigned int attributes) {
  unsigned int n_sites = common_args->n_sites;
  unsigned int n_categories = common_args->n_categories;

  pll_partition_t *partition;
  partition = pll_partition_create(
          n_tips,                 /* numer of tips */
          4,                      /* clv buffers */
          n_states,               /* number of states */
          n_sites,                /* sequence length */
          1,                      /* different rate parameters */
          2 * n_tips - 3,         /* probability matrices */
          n_categories,           /* gamma categories */
          0,                      /* scale buffers */
          attributes /* attributes */
  );

  if (!partition) {
    fatal("Fail creating partition");
  }

  double branch_lengths[4] = {0.1, 0.2, 0.3, 0.4};
  unsigned int matrix_indices[4] = {0, 1, 2, 3};

  pll_set_frequencies(partition, 0, pll_freqs);
  pll_set_subst_params(partition, 0, pll_rates);

  int return_val = PLL_SUCCESS;
  for (unsigned int i = 0; i < ALIGN_SEQS; i++) {
    return_val &= pll_set_tip_states(partition, i, pll_map,
                                     align[i]);
  }

  if (!return_val)
    fatal("Error setting tip states");

  double *rate_cats = (double *) malloc(n_categories * sizeof(double));

  if (pll_compute_gamma_cats(0.1, n_categories, rate_cats, PLL_GAMMA_RATES_MEAN)
      == PLL_FAILURE) {
    printf("Fail computing the gamma rates\n");
    exit(0);
  }

  pll_set_category_rates(partition, rate_cats);
  free(rate_cats);

  for (unsigned int j = 0; j < partition->rate_matrices; ++j) {
    pll_update_invariant_sites_proportion(partition, j, common_args->pinvar);
  }

  pll_update_prob_matrices(partition, params_indices, matrix_indices, branch_lengths, 4);

  clock_t pll_update_partials_begin_time = clock();

  for (int i = 0; i < common_args->n_itr; i++) {
    for (int o = 0; o < 3; o++) {
      const pll_operation_t *op = &(operations[o]);
      kernel_fun(partition->sites,
                 partition->rate_cats,
                 partition->clv[op->parent_clv_index],
                 NULL,
                 partition->clv[op->child1_clv_index],
                 partition->clv[op->child2_clv_index],
                 partition->pmatrix[op->child1_matrix_index],
                 partition->pmatrix[op->child2_matrix_index],
                 NULL,
                 NULL,
                 partition->attributes);
    }
  }

  clock_t pll_update_partials_end_time = clock();
  float secs = (float) (pll_update_partials_end_time - pll_update_partials_begin_time) / CLOCKS_PER_SEC;

  pll_partition_destroy(partition);
  return secs;
}

// Creates a MSA by choosing randomly characters from the given comm ref_seq
static char **create_random_align(const char *ref_seq, unsigned int n_sites, unsigned int print_seq) {

  char **align = calloc(ALIGN_SEQS, sizeof(char *));

  if (print_seq) {
    printf("Random alignment:\n");
  }
  for (unsigned int i = 0; i < ALIGN_SEQS; i++) {
    align[i] = calloc(n_sites + 1, sizeof(char));
    for (unsigned int j = 0; j < n_sites; j++) {
      align[i][j] = ref_seq[rand() % strlen(ref_seq)];
    }
    if (print_seq) {
      printf("  [%d]=%s\n", i, align[i]);
    }
  }
  return align;
}

// deletes a randomly allocated MSA
static void delete_random_align(char **align) {
  for (unsigned int i = 0; i < ALIGN_SEQS; i++) {
    free(align[i]);
  }
  free(align);
}

int main(int argc, char *argv[]) {
  unsigned int n_tips = 5;

  const char *ref_seq_aa = "*-?ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz";
  const char *ref_seq_nt = "-?ABCDGHKMNORSTUVWXYabcdghkmnorstuvwxy";

  double titv = 2.5;
  double pll_nt_freqs[4] = {0.3, 0.4, 0.1, 0.2};
  double pll_nt_rates[6] = {1, titv, 1, 1, titv, 1};

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

  /* check attributes */
  pll_common_args_t *common_args = get_common_args(argc, argv);
  unsigned int n_sites = common_args->n_sites;
  unsigned int n_categories = common_args->n_categories;

  srand(common_args->seed);

  char **align_aa = create_random_align(ref_seq_aa, n_sites, common_args->print_seq);
  char **align_nt = create_random_align(ref_seq_nt, n_sites, common_args->print_seq);

  #define MM_DAZ_ON    0x0040
  _mm_setcsr(_mm_getcsr() | (_MM_FLUSH_ZERO_ON | MM_DAZ_ON));

  float avx2_aa_secs = benchmark_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                                        20, common_args, operations,
                                        &pll_core_update_partial_ii_20x20_avx2,
                                        PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2);

  float avx2sml_aa_secs = benchmark_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                                           20, common_args, operations,
                                           &pll_core_update_partial_ii_20x20_avx2_sml,
                                           PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2 |
                                           PLL_ATTRIB_SIMD_MEM_LAYOUT);

  float avx512_aa_secs = benchmark_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                                          20, common_args, operations,
                                          &pll_core_update_partial_ii_20x20_avx512f,
                                          PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX512F);

  float avx512sml_aa_secs = benchmark_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                                             20, common_args, operations,
                                             &pll_core_update_partial_ii_20x20_avx512f_sml,
                                             PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX512F |
                                             PLL_ATTRIB_SIMD_MEM_LAYOUT);

  float avx2_nt_secs = benchmark_kernel(n_tips, align_nt, pll_nt_freqs, pll_nt_rates, pll_map_nt, 4, common_args,
                                        operations,
                                        &pll_core_update_partial_ii_4x4_avx,
                                        PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX);

  float avx2sml_nt_secs = benchmark_kernel(n_tips, align_nt, pll_nt_freqs, pll_nt_rates, pll_map_nt, 4, common_args,
                                           operations,
                                           &pll_core_update_partial_ii_4x4_avx2_sml,
                                           PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2 |
                                           PLL_ATTRIB_SIMD_MEM_LAYOUT);

  float avx512sml_nt_secs = benchmark_kernel(n_tips, align_nt, pll_nt_freqs, pll_nt_rates, pll_map_nt, 4, common_args,
                                             operations,
                                             &pll_core_update_partial_ii_4x4_avx512f_sml,
                                             PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX512F |
                                             PLL_ATTRIB_SIMD_MEM_LAYOUT);

  setlocale(LC_NUMERIC, "");
  printf("Seed:                                  %d\n", common_args->seed);
  printf("No. Sites:                             %'d\n", n_sites);
  printf("No. Categories:                        %d\n", n_categories);
  printf("No. Iterations:                        %d\n", common_args->n_itr);
  printf("Proportion of invariant Sites:         %f\n", common_args->pinvar);
  printf("CLV buffer size (AA):                  %.1f KB\n", common_args->n_sites
                                                             * common_args->n_categories * 20 *sizeof(double)/1024.0);
  printf("CLV buffer size (NT):                  %.1f KB\n", common_args->n_sites
                                                           * common_args->n_categories * 4 *sizeof(double)/1024.0);
  printf("\n");
  printf("pll_update_partials AA (AVX2):         %f\n", avx2_aa_secs);
  printf("pll_update_partials AA (AVX2 SML):     %f %f (vs AVX2)\n", avx2sml_aa_secs, avx2_aa_secs/avx2sml_aa_secs);
  printf("pll_update_partials AA (AVX512F):      %f %f (vs AVX2)\n", avx512_aa_secs, avx2_aa_secs/avx512_aa_secs);
  printf("pll_update_partials AA (AVX512F SML):  %f %f (vs AVX2)\n", avx512sml_aa_secs, avx2_aa_secs/avx512sml_aa_secs);
  printf("\n");
  printf("pll_update_partials NT (AVX):          %f\n", avx2_nt_secs);
  printf("pll_update_partials NT (AVX2 SML):     %f %f (vs AVX)\n", avx2sml_nt_secs, avx2_nt_secs/avx2sml_nt_secs);
  printf("pll_update_partials NT (AVX512F SML):  %f %f (vs AVX)\n", avx512sml_nt_secs, avx2_nt_secs/avx512sml_nt_secs);

#ifdef SIMD_COUNTING_MODE

  printf("Instruction counts for pll_core_update_partial_ii_20x20_avx2:\n");
  count_instrs_of_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                         20, common_args, operations,
                         &pll_core_update_partial_ii_20x20_avx2,
                         PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2);

  printf("Instruction counts for pll_core_update_partial_ii_20x20_avx2_sml:\n");
  count_instrs_of_kernel(n_tips, align_aa, pll_aa_freqs_dayhoff, pll_aa_rates_dayhoff, pll_map_aa,
                         20, common_args, operations,
                         &pll_core_update_partial_ii_20x20_avx2_sml,
                         PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2 |
                         PLL_ATTRIB_SIMD_MEM_LAYOUT);

  printf("Instruction counts for pll_core_update_partial_ii_4x4_avx:\n");
  count_instrs_of_kernel(n_tips, align_nt, pll_nt_freqs, pll_nt_rates, pll_map_nt, 4, common_args,
                         operations,
                         &pll_core_update_partial_ii_4x4_avx,
                         PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX);

  printf("Instruction counts for pll_core_update_partial_ii_4x4_avx2_sml:\n");
  count_instrs_of_kernel(n_tips, align_nt, pll_nt_freqs, pll_nt_rates, pll_map_nt, 4, common_args,
                         operations,
                         &pll_core_update_partial_ii_4x4_avx2_sml,
                         PLL_ATTRIB_ARCH_CPU | PLL_ATTRIB_ARCH_AVX2 |
                         PLL_ATTRIB_SIMD_MEM_LAYOUT);

#endif /* SIMD_COUNTING_MODE */

  delete_random_align(align_aa);
  delete_random_align(align_nt);

  free(operations);
  destroy_common_args(&common_args);
  return (0);
}
