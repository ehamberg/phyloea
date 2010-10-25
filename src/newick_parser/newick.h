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

/**
 * The data structure which holds nodes in the tree. 
 *
 */
struct _newick_node {
    double       dist;
    char*        name;              
    int          _n_children;
    int          _max_children;
    struct _newick_node*      parent;
    struct _newick_node**     children;
};

typedef struct _newick_node _newick_node_t; 


_newick_node_t* _newick_tree_;
