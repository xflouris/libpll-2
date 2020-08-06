/*
    Copyright (C) 2020 Alexey Kozlov

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

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

// for asprintf
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include "common.h"

#include <stdarg.h>
#include <search.h>

#define DATATYPE_NT 0
#define DATATYPE_AA  1
#define DATATYPE_ODD 2

#define N_STATES_NT 4
#define N_STATES_AA 20
#define N_STATES_ODD 7

/* odd map with 7 states: A..G */
const pll_state_t odd7_map[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x3f, 0, 0, 0x3f, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x3f, 0, 0x01, 0x02, 0x04,
    0x08, 0x0c, 0x10, 0x20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0x01, 0x02, 0x04, 0x08, 0x0c, 0x10, 0x20, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


#define MSA_SEQCOUNT 5

const char* nt_msa[]  = {"WWAC-CTA-ATACTT",
                         "CCCC-TTA-ATAGTT",
                         "AA-C-TAG-CTGCTT",
                         "CCTCTTAA-A-ACGG",
                         "CCAC-TCA-A-ATGG"};

const char* aa_msa[]  = {"PIGL-RVTLGRRGD-RMWI",
                         "IQGM-DITIGVTG------",
                         "--AF-ALLQAKIAG-MPFE",
                         "MDIS-IVT-I-------TA",
                         "GLSE-QTVFSHESI-DQDK"};

const char* odd_msa[] = {"AAB-CCD-EAFAA",
                         "ACC-FBA-ACBGG",
                         "A-C-GAG-G-CCF",
                         "ADCFCAA-AD-CG",
                         "ABC-BCA-AB-BG"};

pll_msa_t * copy_msa(const char* sequence[], unsigned int count, unsigned int len)
{
  unsigned int i;
  pll_msa_t * msa = (pll_msa_t *) calloc(1, sizeof(pll_msa_t));

  msa->count = count;
  msa->length = len ? len : strlen(sequence[0]);

  msa->sequence = (char**) malloc(count * sizeof(char*));
  for (i = 0; i < count; ++i)
  {
    msa->sequence[i] = (char*) malloc((msa->length+1) * sizeof(char));
    memcpy(msa->sequence[i], sequence[i], (msa->length+1)*sizeof(char));
  }

  return msa;
}

static void print_msa(pll_msa_t * msa)
{
  unsigned int i,j;

  for (i = 0; i < msa->count; ++i)
  {
    for (j = 0; j < msa->length; ++j)
      printf("%c", msa->sequence[i][j]);
    printf("\n");
  }
}

static void print_uncompressed_msa(pll_msa_t * msa, unsigned int * site_pattern_map,
                                   unsigned int uncomp_len)
{
  unsigned int i,j,pat;

  for (i = 0; i < msa->count; ++i)
  {
    for (j = 0; j < uncomp_len; ++j)
    {
      pat = site_pattern_map[j];
      assert(pat < msa->length);
      printf("%c", msa->sequence[i][pat]);
    }
    printf("\n");
  }
}

static void test_compress(int datatype, pll_bool_t backmap)
{
  unsigned int i;

  unsigned int seqlen = 0;
  const pll_state_t * map = NULL;
  pll_msa_t * msa = NULL;
  unsigned int * site_pattern_map = NULL;
  unsigned int * w = NULL;
  const char * dt_name = NULL;
  unsigned int sumw = 0;

  switch(datatype)
  {
    case DATATYPE_NT:
      map = pll_map_nt;
      msa = copy_msa(nt_msa, MSA_SEQCOUNT, 0);
      dt_name = "DNA";
      break;
    case DATATYPE_AA:
      map = pll_map_aa;
      msa = copy_msa(aa_msa, MSA_SEQCOUNT, 0);
      dt_name = "AA";
      break;
    case DATATYPE_ODD:
      map = odd7_map;
      msa = copy_msa(odd_msa, MSA_SEQCOUNT, 0);
      dt_name = "ODD7";
      break;
    default:
      assert(0);
  }

  printf("* TEST: DATATYPE = %s, BACKMAP = %s\n\n",
         dt_name, backmap ? "YES" : "NO");

  seqlen = msa->length;

  if (backmap)
    site_pattern_map = (unsigned int *) calloc(seqlen, sizeof(unsigned int));

  printf("ORIGINAL MSA (%u):\n", msa->length);
  print_msa(msa);

  printf("\n");

  w = pll_compress_site_patterns_msa(msa, map, site_pattern_map);

  if (w)
  {
    printf("COMPRESSED MSA (%u):\n", msa->length);
    print_msa(msa);

    printf("\n");

    printf("PATTERN WEIGHTS: ");
    for (i = 0; i < msa->length; ++i)
    {
      sumw += w[i];
      printf("%u ", w[i]);
    }
    printf("\n");

    assert(sumw == seqlen);

    if (site_pattern_map)
    {
      printf("SITE-TO-PATTERN MAP: ");
      for (i = 0; i < seqlen; ++i)
        printf("%u ", site_pattern_map[i]);
      printf("\n\n");

      printf("UNCOMPRESSED MSA (%u):\n", sumw);
      print_uncompressed_msa(msa, site_pattern_map, sumw);
      printf("\n");
    }
  }
  else
   printf("Pattern compression failed: ERR-%d  %s\n ", pll_errno, pll_errmsg);

  printf("\n\n");

  pll_msa_destroy(msa);
  free(w);
  free(site_pattern_map);
}

int main(int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);

  if (attributes != PLL_ATTRIB_ARCH_CPU)
    skip_test();

  test_compress(DATATYPE_NT, PLL_FALSE);
  test_compress(DATATYPE_NT, PLL_TRUE);

  test_compress(DATATYPE_AA, PLL_FALSE);
  test_compress(DATATYPE_AA, PLL_TRUE);

  test_compress(DATATYPE_ODD, PLL_FALSE);
  test_compress(DATATYPE_ODD, PLL_TRUE);

  return (EXIT_SUCCESS);
}
