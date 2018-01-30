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
 * \function dgraph_getkcore 
 * \brief Finding the induced subgraph for core value k in a network.
 *
 * The k-core of a graph is a maximal subgraph in which each vertex
 * has at least degree k. (Degree here means the degree in the
 * subgraph of course.). The coreness of a vertex is the highest order 
 * of a k-core containing the vertex. getkcore computes this maximal 
 * subgraph containing nodes with degree greater than equal to k.
 * 
 * </para><para>
 * \param graph The input graph.
 * \param mode For directed graph it specifies whether to calculate
 *        in-cores, out-cores or the undirected version. It is ignored
 *        for undirected graphs. Possible values: \c IGRAPH_ALL
 *        undirected version, \c IGRAPH_IN in-cores, \c IGRAPH_OUT
 *        out-cores. 
 * \param k Nodes with core value greater than equal to k will be included in the subgraph.
 * \return Error code.
 *
 * Time complexity: O(|E|), the number of edges.
 */

int dgraph_getkcore(const igraph_t *graph, igraph_vector_t *res,
		    igraph_neimode_t mode, igraph_integer_t k) {
  
  igraph_vector_t cores;
  long int maxk, i, j, c;
  //long int *vids;
 // igraph_vector_t vids;

  igraph_vector_init(&cores, 0);
 // igraph_vector_init(&vids, 0);
  long int no_of_nodes=igraph_vcount(graph);
 // Rprintf("\nhi");
  igraph_coreness(graph, &cores, mode);
 // Rprintf("\nm");
  maxk = (long int) igraph_vector_max(&cores);
 // Rprintf("\nmaxk:%li",maxk);
 // IGRAPH_VECTOR_INIT_FINALLY(&vids,100);
  if(k <= maxk){

	j=0;
        for(i=0; i<no_of_nodes; i++){
		c = (long int)VECTOR(cores)[i];
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
//	Rprintf("\nno: %li",no_of_nodes);
	j=0;
        for(i=0; i<no_of_nodes; i++){
		c = (long int)VECTOR(cores)[i];
		if(c >= k){
			VECTOR(*res)[j] =i+1;
			//VECTOR(vids)[j] = i;
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
// 	Rprintf("\nafterinduced"); 
 }

  igraph_vector_destroy(&cores);
//  igraph_vector_destroy(&vids);
  IGRAPH_FINALLY_CLEAN(1);
//Rprintf("\nafter"); 
//igraph_vector_destroy(&vids);
 // IGRAPH_FINALLY_CLEAN(1);
  //igraph_free(vids);
  //igraph_vs_destroy(&vs);
  //IGRAPH_FINALLY_CLEAN(1);
 
 return 0;
}
