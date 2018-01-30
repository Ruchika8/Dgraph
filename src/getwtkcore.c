/* -*- mode: C -*-  */
/* 
   DGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph_community.h"
#include "igraph_memory.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "config.h"
#include "R.h"

/**
 * \function igraph_getkcore 
 * \brief Finding the induced subgraph for core k in a network.
 *
 * The k-core of a graph is a maximal subgraph in which each vertex
 * has at least degree k. (Degree here means the degree in the
 * subgraph of course.). The coreness of a vertex is the highest order 
 * of a k-core containing the vertex.
 * 
 * </para><para>
 * This function implements the algorithm presented in Vladimir
 * Batagelj, Matjaz Zaversnik: An O(m) Algorithm for Cores
 * Decomposition of Networks. 
 * \param graph The input graph.
 * \param cores Pointer to an initialized vector, the result of the
 *        computation will be stored here. It will be resized as
 *        needed. For each vertex it contains the highest order of a
 *        core containing the vertex.
 * \param mode For directed graph it specifies whether to calculate
 *        in-cores, out-cores or the undirected version. It is ignored
 *        for undirected graphs. Possible values: \c IGRAPH_ALL
 *        undirected version, \c IGRAPH_IN in-cores, \c IGRAPH_OUT
 *        out-cores. 
 * \return Error code.
 *
 * Time complexity: O(|E|), the number of edges.
 */

int dgraph_getwtkcore(const igraph_t *graph, igraph_vector_t *res, igraph_neimode_t mode, igraph_integer_t k,
		      igraph_vector_t *weights, igraph_integer_t type, igraph_integer_t alpha,
	 	      igraph_integer_t beta, igraph_real_t lambda, igraph_integer_t wtdeg) {
  
  igraph_vector_t cores;
  long int maxk, i, j, c;
  //long int *vids;
 // igraph_vector_t vids;

  igraph_vector_init(&cores, 0);
 // igraph_vector_init(&vids, 0);
  long int no_of_nodes=igraph_vcount(graph);
 // Rprintf("\nhi");

  if(wtdeg==1)
  	dgraph_wtcoreness1(graph, &cores, mode, weights, type, alpha, beta, lambda);
  else
        dgraph_wtcoreness2(graph, &cores, mode, weights, type, alpha, beta, lambda);

 // Rprintf("\nm");
  maxk = (long int) igraph_vector_max(&cores);
//  Rprintf("\nmaxk:%li",maxk);
 // IGRAPH_VECTOR_INIT_FINALLY(&vids,100);
  if(k <= maxk){

	j=0;
        for(i=0; i<no_of_nodes; i++){
		c = (long int)VECTOR(cores)[i];
	//	Rprintf("\ncore:%li",c);
		if(c >= k){	
			j++;
		}
	}
	IGRAPH_CHECK(igraph_vector_resize(res, j));
  	igraph_vector_null(res);

	/*vids=igraph_Calloc(no_of_nodes, long int);
  	if (vids==0) {
    		IGRAPH_ERROR("Cannot calculate", IGRAPH_ENOMEM);
  	}
  	IGRAPH_FINALLY(igraph_free, vids);
	*/
	//Rprintf("\nno: %li",no_of_nodes);
	j=0;
        for(i=0; i<no_of_nodes; i++){
		c = (long int)VECTOR(cores)[i];
		if(c >= k){
			VECTOR(*res)[j] = i+1;
			//Rprintf("\nvid:%li,%li",(long int)VECTOR(vids)[j],i);
			//Rprintf("\nsizevectinloop: %li", igraph_vector_size(&vids));		
			j++;
		}
	}
//	Rprintf("\nbeforeinduced");
//	igraph_vs_t vs;
//		Rprintf("\nballe");
 //	igraph_vs_vector_copy(&vs,&vids);
  //	Rprintf("\ninbetween");
	//Rprintf("\ncores: %li", igraph_vector_size(&cores));
//	Rprintf("\nsizevect: %li", igraph_vector_size(&vids));
//	igraph_induced_subgraph(graph, res, igraph_vss_vector(&vids), IGRAPH_SUBGRAPH_AUTO);
 //	Rprintf("\nafterinduced"); 
 }

  igraph_vector_destroy(&cores);
 // igraph_vector_destroy(&vids);
  IGRAPH_FINALLY_CLEAN(1);

//igraph_vector_destroy(&vids);
 // IGRAPH_FINALLY_CLEAN(1);
  //igraph_free(vids);
  //igraph_vs_destroy(&vs);
  //IGRAPH_FINALLY_CLEAN(1);
 
 return 0;
}
