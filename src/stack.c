/*
    Copyright (C) 2020 Pierre Barbera, HITS gGmbH

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

    Contact: Alexandros Stamatakis <alexandros.stamatakis@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

pll_uint_stack_t* pll_uint_stack_create(const size_t size)
{
  assert(size);

  // alloc the struct
  pll_uint_stack_t* stack = (pll_uint_stack_t*)malloc(sizeof(pll_uint_stack_t));

  // alloc the data
  stack->data = (unsigned int *)malloc(size * sizeof(unsigned int));
  if (!stack->data)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg,
             200,
             "Unable to allocate enough memory for pll_stack.");
    return PLL_FAILURE;
  }

  // set other vars
  stack->size = size;
  stack->top  = stack->data - 1;
  stack->empty= true;

  return stack;
}

void pll_uint_stack_destroy(pll_uint_stack_t* stack)
{
  if (stack)
  {
    free(stack->data);
    free(stack);
  }
}

unsigned int pll_uint_stack_push(pll_uint_stack_t* stack,
                                 const unsigned int val)
{
  assert(stack);

  if (stack->top + 1 >= stack->data + stack->size)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg,
             200,
             "Trying to push to a full pll_stack.");
    return PLL_FAILURE;
  }

  ++stack->top;

  *(stack->top) = val;

  // we are pushing, so the stack can't be empty afterwards
  stack->empty = false;
  return PLL_SUCCESS;
}

unsigned int pll_uint_stack_pop(pll_uint_stack_t* stack)
{
  assert(stack);

  if (stack->top < stack->data)
  {
    printf("Trying to pop from an empty pll_stack.");
    exit(PLL_FAILURE);
  }

  // if the element being popped is the last remaining, the stack will be empty
  stack->empty = (stack->top == stack->data);

  unsigned int ret = *stack->top;
  stack->top--;

  return ret;
}
