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
#include <x86intrin.h>

double test_avx2(const double *A, const double *B, double *C, size_t elem ,size_t max_repeat);
double test_avx512(const double *A, const double *B, double *C, size_t elem, size_t max_repeat);

typedef double (*test_func)(const double *, const double *, double *, size_t, size_t);

double *createArray(size_t dim, size_t align) {
  double *p;
  if (posix_memalign((void **) &p, align, dim * sizeof(double))) {
    fprintf(stderr, "Unable to allocate memory block");
    exit(-1);
  }

  for (size_t i = 0; i < dim; i++) {
    p[i] = (rand() % 100) / 100000.0;
  }

  return p;
}

void benchmark_test_func(size_t align, unsigned int seed, test_func func, size_t elems, size_t repeats) {
  srand(seed);

  double *A = createArray(elems, align);
  double *B = createArray(elems, align);
  double *C = createArray(elems, align);

//  for (size_t i = 0; i < ARRAY_SIZE; i++) {
//    printf("% .3e,", A[i]);
//  }
//  printf("\n");
//  for (size_t i = 0; i < ARRAY_SIZE; i++) {
//    printf("% .3e,", B[i]);
//  }
//  printf("\n");
//  for (size_t i = 0; i < ARRAY_SIZE; i++) {
//    printf("% .3e,", C[i]);
//  }
//  printf("\n");

  clock_t begin_time = clock();

  double r = func(A, B, C, elems, repeats);

  clock_t end_time = clock();
  float secs = (float) (end_time - begin_time) / CLOCKS_PER_SEC;

  printf("Elapsed time: % .6f,\n", secs);

//  for (size_t i = 0; i < ARRAY_SIZE; i++) {
//    printf("% .3e,", C[i]);
//  }
//  printf("\n");

  printf("Check: %f\n", r);

  free(A);
  free(B);
  free(C);
}

int main(int argc, char *argv[]) {
  unsigned int seed = (unsigned int) time(NULL);

  if(argc < 3) {
    printf("Usage: %s <array_factor> <n_repeats>\n",argv[0]);
    return -1;
  }

  size_t elems = (unsigned int)abs(atoi(argv[1]));
  size_t repeats = (unsigned int)abs(atoi(argv[2]));

  printf("Array elements %lu*8 (= %f KB)\n", elems, elems*8*sizeof(double)/1024.0);
  printf("Repeated calculation %lu\n", repeats);

  benchmark_test_func(32, seed, test_avx2, elems*8, repeats);
  benchmark_test_func(64, seed, test_avx512, elems*8, repeats);

  return (0);
}
