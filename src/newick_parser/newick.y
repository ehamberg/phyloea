/**
   Copyright James H. Bullard 2009

   This program is free software: you can redistribute it and/or modify
   it under the terms of the Lesser GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   
   You should have received a copy of the Lesser GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "newick.h"

#ifndef PRINT
#define PRINT printf
#endif

_newick_node_t* make_node(char* name); 
void yyerror(const char *str); 
int yywrap();
int yylex(void);  
void _newick_tree_print(_newick_node_t*);
void _newick_add_child(_newick_node_t*, _newick_node_t*);
 

%}

%union {
       double   numberVal;
       char*    nameVal;
       _newick_node_t*  nPtr;
}

/* Bison declarations.  */
%token <numberVal>   NUMBER
%token <nameVal>     WORD
%type <nPtr>         tree subtree leaf internal branch branchlist 
%type <nameVal>      name
%type <numberVal>    length 
     
%%  /* The grammar follows.  */
tree: 
     subtree ';'			        	{ 
	                                                  /* printf("tree\n"); */
	 
		                                          _newick_tree_print($1); 
							  $$ = $1; 
							  _newick_tree_ = $1;
							  return(1); 
							}
              ;

subtree: leaf						{ 
                                                          /* printf("subtree-leaf\n"); */
	 						  $$ = $1; 
							}
       | internal					{ 
	                                                  /* printf("subtree-internal\n"); */
       	 						  $$ = $1; 
							}
       ;

leaf: name						{ 
                                                          /* printf("leaf\n"); */
							  _newick_node_t* tmp = make_node($1); 
							  $$ = tmp;
							}
    ;

internal: '(' branchlist ')' name       		{ 
                                                          /* printf("internal\n"); */
							  $2->name = $4;
                                                          $$ = $2;
							}
        ;

branchlist: branch					{ 
                                                          /* printf("branchlist-branch\n");*/
							  _newick_node_t* tmp = make_node("");
							  _newick_add_child(tmp, $1);
	    						  $$ = tmp; 
                                                        }
          | branchlist  ','  branch  		        { 
	                                                  /* printf("branchlist-branch,branchlist\n"); */
							  _newick_add_child($1, $3);
	    	       					  $$ = $1; 
							 }
          ;

branch: subtree length					{ 
                                                          /* printf("branch\n"); */
							  $1->dist = $2; 
							  $$ = $1; 
							}
      ;

name: /* empty */					{ $$ = ""; }
    | WORD						{ 
	                                                   /* printf("WORD\n");*/
							   $$ = $1; 
							}
    ;

length: /* empty */					{ $$ = 0; }
      | ':' NUMBER					{ 
	                                                  /* printf("NUMBER\n"); */
							  $$ = $2; 
							}
      ;
%%


_newick_node_t* make_node (char* name) {
	_newick_node_t* ret = (_newick_node_t*) malloc(sizeof(_newick_node_t));
	ret->name                = name;
	ret->dist                = 0;
	ret->_n_children         = 0;
	ret->_max_children       = 0;
	ret->children            = NULL; 
	ret->parent              = NULL;
	return ret;
}
 
void _newick_add_child(_newick_node_t* parent, _newick_node_t* child) {
    if (parent->children == NULL) {
	parent->children = (_newick_node_t**) malloc (sizeof(_newick_node_t*) * 4);
	parent->_max_children = 4;
    }
    else if (parent->_n_children == parent->_max_children) {
	_newick_node_t** tmp = (_newick_node_t**) malloc (sizeof(_newick_node_t*) * 2 * parent->_n_children);
	int i; 
	for (i = 0; i < parent->_n_children; i++) {
	    tmp[i] = parent->children[i];
	}
	parent->children = tmp; 
	parent->_max_children = (2 * parent->_n_children);
	free(parent->children);
    }

    child->parent = parent; 
    parent->children[parent->_n_children++] = child;

}



void _print_node(_newick_node_t* node) {
    if (node->parent == NULL) {
	PRINT("%s", node->name);
    }
    else {
	if (node->dist == 0) 
	    PRINT("%s", node->name);
	else
	    PRINT("%s:%g", node->name, node->dist);
    }
}

void _newick_tree_print_i (_newick_node_t* node) {
    if (node != NULL) {
	if (node->_n_children == 0) {
	    _print_node(node);
	}
	else {
	    PRINT("(");
	    int i; 
	    for (i = 0; i < (node->_n_children - 1); i++) {
		_newick_tree_print_i(node->children[i]);
		PRINT(", ");
	    }
	    _newick_tree_print_i(node->children[i]);
	    PRINT(")"); 
	    _print_node(node);
	}
    }
}

void _newick_tree_print(_newick_node_t* node) {
     printf("\n");
     _newick_tree_print_i(node);
     printf(";\n");
}
  
main()
{
    printf(" ---- Welcome To the Newick Parser ---- \n");
    while (1) yyparse();
} 

void yyerror(const char *str)
{
    fprintf(stderr,"error: %s\n",str);
}
 
int yywrap()
{
        return 1;
} 
