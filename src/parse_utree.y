/*
    Copyright (C) 2015-2018 Tomas Flouri, Alexey Kozlov

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
%{
#include "pll.h"

extern int pll_utree_lex();
extern FILE * pll_utree_in;
extern void pll_utree_lex_destroy();
extern int pll_utree_lineno;
extern int pll_utree_colstart;
extern int pll_utree_colend;

extern int pll_utree_parse();
extern struct pll_utree_buffer_state * pll_utree__scan_string(const char * str);
extern void pll_utree__delete_buffer(struct pll_utree_buffer_state * buffer);

static unsigned int tip_cnt = 0;

static pll_unode_t * alloc_node()
{
  pll_unode_t * node = (pll_unode_t *)calloc(1, sizeof(pll_unode_t)); 
  return node;
}

static void dealloc_data(pll_unode_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

static void close_roundabout(pll_unode_t * first)
{
  pll_unode_t * last = first;
  while(last->next != NULL && last->next != first)
  {
  	if (!last->next->label)
  	  last->next->label = last->label;
  	last = last->next;
  }
  last->next = first;
}

static void dealloc_graph_recursive(pll_unode_t * node,
                                   void (*cb_destroy)(void *), 
                                   int level)
{
  if (!node->next)
  {
    /* tip node */
    dealloc_data(node, cb_destroy);
    free(node->label);
    free(node);
  }
  else
  {
    /* inner node */
    if (node->label)
      free(node->label);
    
    pll_unode_t * snode = node;
	do
    {
      if (node != snode || level == 0)
        dealloc_graph_recursive(snode->back, cb_destroy, level+1);
      pll_unode_t * next = snode->next;
      dealloc_data(snode, cb_destroy);
      free(snode);
      snode = next;
    }
    while(snode && snode != node);
  }
}

PLL_EXPORT void pll_utree_graph_destroy(pll_unode_t * root,
                                        void (*cb_destroy)(void *))
{
  if (!root) return;
  
  dealloc_graph_recursive(root, cb_destroy, 0);
}

PLL_EXPORT void pll_utree_destroy(pll_utree_t * tree,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;

  /* deallocate tip nodes */
  for (i = 0; i < tree->tip_count; ++i)
  {
    dealloc_data(tree->nodes[i], cb_destroy);
    if (tree->nodes[i]->label)
      free(tree->nodes[i]->label);
    free(tree->nodes[i]);
  }

  /* deallocate inner nodes */
  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
  {
    pll_unode_t * first = tree->nodes[i];
    
    assert(first);

    if (first->label)
      free(first->label);

    pll_unode_t * node = first;
	do
    {
      pll_unode_t * next = node->next;
      dealloc_data(node, cb_destroy);
      free(node);
      node = next;
    }
    while(node && node != first);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void pll_utree_error(pll_unode_t * node, const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pll_utree_colstart == pll_utree_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pll_utree_lineno, pll_utree_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n",
             s, pll_utree_lineno, pll_utree_colstart, pll_utree_colend);
}

%}

%union
{
  char * s;
  char * d;
  struct pll_unode_s * tree;
}

%error-verbose
%parse-param {struct pll_unode_s * tree}
%destructor { pll_utree_graph_destroy($$,NULL); } subtree
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label


%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<tree> subtree descendant_list_item descendant_list
%start input
%%

input: descendant_list optional_label optional_length SEMICOLON
{
 tree->back = $1->back;
 $1->back->back = tree;
 tree->next = $1->next;
 tree->node_index = $1->node_index;
 tree->length = $1->length;
 tree->label = $2;
 close_roundabout(tree);
 free($1);
 /* ignore root length if specified -> we create an unrooted tree structure! */
 if ($3)
   free($3);
};
	
descendant_list: OPAR  descendant_list_item CPAR
{
  $$=$2;
};
	
descendant_list_item: subtree
{
  /* create inner node (1st subtree)  */
  $$ = alloc_node();
  $$->back = $1;
  $1->back = $$;
  $$->length = $1->length;
}
	| descendant_list_item COMMA subtree
{
  $$=$1;
  pll_unode_t * last = $1;
  while(last->next != NULL)
  {
    last = last->next;
  }
  last->next = alloc_node();
  last->next->label = last->label;
  last->next->length = $3->length;
  last->next->back = $3;
  $3->back = last->next;
};

subtree : descendant_list optional_label optional_length
{
  /* create internal node */
  $$ = alloc_node();
  $$->next = $1;
  $$->label=$2;
  if ($3)
  {
    $$->length = atof($3);
    free($3);
  }
  else
    $$->length = 0;
  
  close_roundabout($$);
}
         | label optional_length
{
  /* create tip node */
  $$ = alloc_node();
  $$->label = $1;
  if ($2)
  {
    $$->length = atof($2);
    free($2);
  }
  else
    $$->length = 0;
  
  tip_cnt++;
};

optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;} | NUMBER {$$=$1;};
number: NUMBER   { $$=$1;};

%%

static void recursive_assign_indices(pll_unode_t * node,
                                    unsigned int * tip_clv_index,
                                    unsigned int * inner_clv_index,
                                    int * inner_scaler_index,
                                    unsigned int * inner_node_index,
                                    unsigned int level)
{
  if (!node->next)
  {
    /* tip node */
    node->node_index = *tip_clv_index;
    node->clv_index = *tip_clv_index;
    node->pmatrix_index = *tip_clv_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    *tip_clv_index = *tip_clv_index + 1;
  }
  else
  {
    /* inner node */
    pll_unode_t * snode = level ? node->next : node;
    do 
    {
      recursive_assign_indices(snode->back,
                               tip_clv_index,
                               inner_clv_index,
                               inner_scaler_index,
                               inner_node_index,
                               level+1);
      snode = snode->next;
    }
    while (snode != node);

    snode = node;
    do 
    {
      snode->node_index = (*inner_node_index)++;
      snode->clv_index = *inner_clv_index;
      snode->scaler_index = *inner_scaler_index;
      if (snode == node && level > 0)
      	snode->pmatrix_index = *inner_clv_index; 
      else
      	snode->pmatrix_index = snode->back->pmatrix_index;
      snode = snode->next;
    }
    while (snode != node);
    
    *inner_clv_index += 1;
    *inner_scaler_index += 1;
  }
}

PLL_EXPORT void pll_utree_reset_template_indices(pll_unode_t * root,
                                                 unsigned int tip_count)
{
  unsigned int tip_clv_index = 0;
  unsigned int inner_clv_index = tip_count;
  unsigned int inner_node_index = tip_count;
  int inner_scaler_index = 0;

  if (!root->next)
    root = root->back;

  recursive_assign_indices(root,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index,
                           0);
}

static void fill_nodes_recursive(pll_unode_t * node,
                                 pll_unode_t ** array,
                                 unsigned int array_size,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index,
                                 unsigned int level)
{
  unsigned int index;
  if (!node->next)
  {
    /* tip node */
    index = *tip_index;
    *tip_index += 1;
  }
  else
  {
    /* inner node */
    pll_unode_t * snode = level ? node->next : node;
    do 
    {
      fill_nodes_recursive(snode->back, array, array_size, tip_index, 
                           inner_index, level+1);
      snode = snode->next;
    }
    while (snode != node);

    index = *inner_index;
    *inner_index += 1;
  }

  assert(index < array_size);
  array[index] = node;
}

static unsigned int utree_count_nodes_recursive(pll_unode_t * node, 
                                                unsigned int * tip_count,
                                                unsigned int * inner_count,
                                                unsigned int level)
{
  if (!node->next)
  {
    *tip_count += 1;
    return 1;
  }
  else
  {
    unsigned int count = 0;

    pll_unode_t * snode = level ? node->next : node;
	do 
	{
	  count += utree_count_nodes_recursive(snode->back, tip_count, inner_count, level+1);
	  snode = snode->next;
	}
	while (snode != node);

    *inner_count += 1;
	
	return count + 1;
  }
}

static unsigned int utree_count_nodes(pll_unode_t * root, unsigned int * tip_count,
                                      unsigned int * inner_count)
{
  unsigned int count = 0;
  
  if (tip_count)
    *tip_count = 0; 
  
  if (inner_count)
    *inner_count = 0; 

  if (!root->next && !root->back->next)
    return 0;

  if (!root->next)
    root = root->back;
    
  count = utree_count_nodes_recursive(root, tip_count, inner_count, 0);
  
  if (tip_count && inner_count)
    assert(count == *tip_count + *inner_count); 

  return count;
}

static pll_utree_t * utree_wraptree(pll_unode_t * root,
                                    unsigned int tip_count,
                                    unsigned int inner_count,
                                    int binary)
{
  unsigned int node_count;
  
  pll_utree_t * tree = (pll_utree_t *)malloc(sizeof(pll_utree_t));
  if (!tree)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  if (tip_count < 3 && tip_count != 0)
  {
    snprintf(pll_errmsg, 200, "Invalid tip_count value (%u).", tip_count);
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }
  
  if (!root->next)
    root = root->back;

  if (binary)
  {
    if (tip_count == 0)
    {
      node_count = utree_count_nodes(root, &tip_count, &inner_count);
      if (inner_count != tip_count - 2)
      {
        snprintf(pll_errmsg, 200, "Input tree is not strictly bifurcating.");
        pll_errno = PLL_ERROR_PARAM_INVALID;
        return PLL_FAILURE;
      }
    }
    else
    {
      inner_count = tip_count - 2;
      node_count = tip_count + inner_count;
    }
  }
  else
  {
    if (tip_count == 0 || inner_count == 0)
      node_count = utree_count_nodes(root, &tip_count, &inner_count);
    else
      node_count = tip_count + inner_count;
  }

  if (!tip_count)
  {
    snprintf(pll_errmsg, 200, "Input tree contains no inner nodes.");
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  tree->nodes = (pll_unode_t **)malloc(node_count*sizeof(pll_unode_t *));
  if (!tree->nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  fill_nodes_recursive(root, tree->nodes, node_count, &tip_index, &inner_index, 0);
 
  assert(tip_index == tip_count);
  assert(inner_index == tip_count + inner_count);

  tree->tip_count = tip_count;
  tree->inner_count = inner_count;
  tree->edge_count = node_count - 1;
  tree->binary = (inner_count == tip_count-2);
  tree->vroot = root;

  return tree;
}

/* wraps/encalupsates the unrooted tree graph into a tree structure
   that contains a list of nodes, number of tips and number of inner
   nodes. If 0 is passed as tip_count, then an additional recrursion
   of the tree structure is done to detect the number of tips */
PLL_EXPORT pll_utree_t * pll_utree_wraptree(pll_unode_t * root,
                                            unsigned int tip_count)
{
  return utree_wraptree(root, tip_count, 0, 1);
}

PLL_EXPORT pll_utree_t * pll_utree_wraptree_multi(pll_unode_t * root,
                                                  unsigned int tip_count,
                                                  unsigned int inner_count)
{
  return utree_wraptree(root, tip_count, inner_count, 0);
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick(const char * filename)
{
  pll_utree_t * tree;

  struct pll_unode_s * root;

  /* reset tip count */
  tip_cnt = 0;

  pll_utree_in = fopen(filename, "r");
  if (!pll_utree_in)
  {
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    return PLL_FAILURE;
  }

  if (!(root = (pll_unode_t *)calloc(1, sizeof(pll_unode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  if (pll_utree_parse(root))
  {
    pll_utree_graph_destroy(root,NULL);
    root = NULL;
    fclose(pll_utree_in);
    pll_utree_lex_destroy();
    return PLL_FAILURE;
  }

  if (pll_utree_in) fclose(pll_utree_in);

  pll_utree_lex_destroy();

  /* initialize clv and scaler indices to the default template */
  pll_utree_reset_template_indices(root, tip_cnt);

  /* wrap tree */
  tree = utree_wraptree(root, 0, 0, 0);

  return tree;
}

PLL_EXPORT pll_utree_t * pll_utree_parse_newick_string(const char * s)
{
  int rc; 
  struct pll_unode_s * root;
  pll_utree_t * tree = NULL;

  /* reset tip count */
  tip_cnt = 0;

  if (!(root = (pll_unode_t *)calloc(1, sizeof(pll_unode_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  struct pll_utree_buffer_state * buffer = pll_utree__scan_string(s);
  rc = pll_utree_parse(root);
  pll_utree__delete_buffer(buffer);

  pll_utree_lex_destroy();

  if (!rc)
  {
    /* initialize clv and scaler indices */
    pll_utree_reset_template_indices(root, tip_cnt);
    
    tree = utree_wraptree(root, 0, 0, 0);
  }
  else
    free(root);

  return tree;
}
