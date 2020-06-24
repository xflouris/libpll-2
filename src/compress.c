/*
    Copyright (C) 2016-2020 Tomas Flouri, Alexey Kozlov

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static void tracked_swap(int i, int j, char ** x, unsigned int * sort_backmap)
{
  PLL_SWAP(x[i],x[j]);
  if (sort_backmap)
    PLL_SWAP(sort_backmap[i],sort_backmap[j]);
}

static void vecswap(int i, int j, int n, char ** x, unsigned int * sort_backmap)
{
  while (n--)
  {
    tracked_swap(i, j, x, sort_backmap);
    ++i; ++j;
  }
}

static void ssort1(char ** x, int n, int depth, unsigned int * sort_backmap,
                   pll_random_state * rstate)
{
  int a, b, c, d, r, v;

  if (n <= 1) return;

  a = pll_random_getint(rstate, n);

  tracked_swap(0, a, x, sort_backmap);

  v = x[0][depth];

  a = b = 1;
  c = d = n-1;

  while (1)
  {
    while (b <= c && (r = x[b][depth]-v) <= 0)
    {
      if (r == 0)
      {
        tracked_swap(a, b, x, sort_backmap);
        ++a;
      }
      ++b;
    }
    while (b <= c && (r = x[c][depth]-v) >= 0)
    {
      if (r == 0)
      {
        tracked_swap(c, d, x, sort_backmap);
        --d;
      }
      --c;
    }
    if (b > c) break;
    tracked_swap(b, c, x, sort_backmap);
    ++b; --c;
  }

  r = PLL_MIN(a, b-a);
  vecswap(0, b-r, r, x, sort_backmap);
  r = PLL_MIN(d-c, n-d-1);
  vecswap(b, n-r, r, x, sort_backmap);
  r = b-a;
  ssort1(x, r, depth, sort_backmap, rstate);

  if (x[r][depth] != 0)
  {
    ssort1 (x + r, a + n - d - 1, depth + 1,
            sort_backmap ? sort_backmap + r : NULL, rstate);
  }

  r = d - c;
  ssort1(x + n - r, r, depth,
         sort_backmap ? sort_backmap + n - r : NULL, rstate);
}

static void remap_range(const pll_state_t * map,
                        unsigned char * charmap)
{
  pll_state_t oldmap[PLL_ASCII_SIZE];
  unsigned int i,j;
  unsigned char k = 1;

  memcpy(oldmap, map, PLL_ASCII_SIZE * sizeof(pll_state_t));
  memset(charmap, 0, PLL_ASCII_SIZE * sizeof(unsigned char));

  for (i = 0; i < PLL_ASCII_SIZE; ++i)
    if (oldmap[i])
    {
      charmap[i] = k;

      for (j = i+1; j < PLL_ASCII_SIZE; ++j)
        if (oldmap[i] == oldmap[j])
        {
          charmap[j] = k;
          oldmap[j] = 0;
        }

      ++k;
    }
}

static pll_state_t findmax(const pll_state_t * map)
{
  int i;
  pll_state_t max = 0;

  for (i = 0; i < PLL_ASCII_SIZE; ++i)
    if (map[i] > max)
      max = map[i];

  return max;
}

static int encode(char ** sequence, const unsigned char * map, int count, int len)
{
  int i,j;
  unsigned char * p;
  unsigned char c;

  /* reset error */
  pll_errno = 0;

  for (i = 0; i < count; ++i)
  {
    p = (unsigned char*)sequence[i];
    j = len;
    while (j--)
    {
      c = map[(int)(*p)];
      if (!c)
      {
        pll_errno = PLL_ERROR_TIPDATA_ILLEGALSTATE;
        snprintf (pll_errmsg, 200,
                  "Cannot encode character %c at sequence %d position %d.",
                  *p,
                  i+1,
                  len-j);
        return PLL_FAILURE;
      }
      *p = c;
      ++p;
    }
  }

  return PLL_SUCCESS;
}

static unsigned int * compress_site_patterns(char ** sequence,
                                             const pll_state_t * map,
                                             int count,
                                             int * length,
                                             unsigned int * site_pattern_map)
{
  int i,j;
  char * memptr;
  char ** column;
  unsigned int * weight;
  unsigned int * sort_backmap;
  pll_random_state * rnd_state;

  unsigned char charmap[PLL_ASCII_SIZE];
  unsigned char inv_charmap[PLL_ASCII_SIZE];

  /* check that at least one sequence is given */
  if (!count)
  {
    pll_errno = PLL_ERROR_MSA_EMPTY;
    snprintf (pll_errmsg, 200,
              "Number of sequences must be greater than 0.");
    return NULL;
  }

  /* a map must be given */
  if (!map)
  {
    pll_errno = PLL_ERROR_MSA_MAP_INVALID;
    snprintf (pll_errmsg, 200,
              "Map is undefined.");
    return NULL;
  }

  /* a zero can never be used as a state */
  if (map[0])
  {
    pll_errno = PLL_ERROR_MSA_MAP_INVALID;
    snprintf (pll_errmsg, 200,
              "'0' cannot be used as a state.");
    return NULL;
  }

  /* if map states are out of the BYTE range, remap */
  if (findmax(map) >= PLL_ASCII_SIZE)
  {
    remap_range(map,charmap);
  }
  else
  {
    for (i = 0; i < PLL_ASCII_SIZE; ++i)
      charmap[i] = (unsigned char)(map[i]);
  }

  /* create inverse charmap to decode states back to characters when
     compression is finished */
  memset(inv_charmap, 0, PLL_ASCII_SIZE);
  for (i = 0; i < PLL_ASCII_SIZE; ++i)
  {
    /* always use '-' to represent gap, otherwise if multiple chars code for the
     * same state, use char with lowest ASCII code (= capital AA/DNA letters) */
    if (map[i] && (!inv_charmap[charmap[i]] || i == '-'))
      inv_charmap[charmap[i]] = (unsigned char)i;
  }

  /* encode sequences using charmap */
  if (!encode(sequence,charmap,count,*length))
  {
    /* spread errno and message */
    assert(pll_errno);
    return NULL;
  }

  if (site_pattern_map)
  {
    sort_backmap = (unsigned int *) malloc((size_t)(*length)*sizeof(unsigned int));
    if (!sort_backmap)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf (pll_errmsg, 200,
                "Cannot allocate space for sort backmap.");
      return NULL;
    }
    for (i = 0; i < *length; ++i)
      sort_backmap[i] = i;
  }
  else
    sort_backmap = NULL;

  /* allocate memory for columns */
  column = (char **)malloc((size_t)(*length)*sizeof(char *));
  if (!column)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for matrix columns.");
    free(sort_backmap);
    return NULL;
  }

  /* allocate memory for the alignment */
  memptr = column[0] = (char *)malloc((size_t)(*length) *
                                      (size_t)(count+1) *
                                      sizeof(char));
  if (!memptr)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf (pll_errmsg, 200,
              "Cannot allocate space for matrix data.");
    free(sort_backmap);
    free(column);
    return NULL;
  }

  /* map memory to each column */
  for (i = 1; i < *length; ++i)
  {
    column[i] = column[i-1] + (count+1);
  }

  /* allocate space for weight vector */
  weight = (unsigned int *)malloc((size_t)(*length)*sizeof(unsigned int));
  if (!weight)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Cannot allocate space for storing site weights.");
    free(sort_backmap);
    free(column);
    free(memptr);
    return NULL;
  }

  /* split alignment into columns instead of rows */
  for (i = 0; i < (*length); ++i)
  {
    for (j = 0; j < count; ++j)
      column[i][j] = sequence[j][i];
    column[i][j] = 0;
  }

  rnd_state = pll_random_create(0);

  if (!rnd_state)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200,
             "Cannot allocate space for storing RNG state.");
    free(weight);
    free(sort_backmap);
    free(column);
    free(memptr);
    return NULL;
  }

  /* sort the columns */
  ssort1(column, *length, 0, sort_backmap, rnd_state);

  /* we have at least one unique site with weight 1 (the first site) */
  int compressed_length = 1;
  size_t ref = 0;
  weight[ref] = 1;
  if (site_pattern_map)
    site_pattern_map[sort_backmap[0]] = ref;

  /* find all unique columns and set their weights */
  for (i = 1; i < *length; ++i)
  {
    if (strcmp(column[i],column[i-1]))
    {
      column[ref+1] = column[i];
      ++ref;
      ++compressed_length;
      weight[ref] = 1;
    }
    else
      weight[ref]++;
    if (site_pattern_map)
      site_pattern_map[sort_backmap[i]] = ref;
  }

  /* copy the unique columns over the original sequences */
  for (i = 0; i < compressed_length; ++i)
    for (j = 0; j < count; ++j)
      sequence[j][i] = column[i][j];

  /* add terminating zero */
  for (j = 0; j < count; ++j)
    sequence[j][compressed_length] = 0;

  /* deallocate memory */
  free(memptr);
  free(column);
  free(sort_backmap);
  pll_random_destroy(rnd_state);

  /* adjust weight vector size to compressed length */
  unsigned int * mem = (unsigned int *)malloc((size_t)compressed_length *
                                              sizeof(unsigned int));
  if (mem)
  {
    /* copy weights */
    for (i = 0; i < compressed_length; ++i)
      mem[i] = weight[i];

    /* free and re-point */
    free(weight);
    weight = mem;
  }


  /* update length */
  *length = compressed_length;

  /* decode sequences using inv_charmap */
  if (!encode(sequence,inv_charmap,count,compressed_length))
  {
    /* if previous enconding was successful, decoding should not fail */
    assert(0);
  }

  return weight;
}

PLL_EXPORT unsigned int * pll_compress_site_patterns(char ** sequence,
                                                     const pll_state_t * map,
                                                     int count,
                                                     int * length)
{
  return compress_site_patterns(sequence, map, count, length, NULL);
}

PLL_EXPORT
unsigned int * pll_compress_site_patterns_msa(pll_msa_t * msa,
                                              const pll_state_t * map,
                                              unsigned int * site_pattern_map)
{
  return compress_site_patterns(msa->sequence, map, msa->count, &msa->length,
                                site_pattern_map);
}
