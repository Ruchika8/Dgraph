/* -*- mode: C -*-  */
/* 
   IGraph library.
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
#include <math.h>
#include "igraph_memory.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "config.h"
#include "R.h"
#include <time.h>

/**
 * \function igraph_wtcoreness 
 * \brief Finding the coreness of the vertices in a network (both
 *  directed as well as undirected).
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
 * \param weights A vector containing the weights of the edges
 *        in the same order as the simple edge iterator visits them
 *        (i.e. in increasing order of edge IDs).
 * \return Error code.
 *

 * Time complexity: O(|E|), the number of edges.
 */

int dgraph_wtcoreness1(const igraph_t *graph, igraph_vector_t *cores, igraph_neimode_t mode, igraph_vector_t *weights, 
		      igraph_integer_t type, igraph_integer_t alpha, igraph_integer_t beta, igraph_real_t lambda) {
/*
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Rprintf ( "Current local time and date: %s", asctime (timeinfo) );
*/
  clock_t start_time,end_time;
  start_time = clock();
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);


  long double a=alpha, b=beta, l=lambda;
  long int *bin, *vert, *pos;
  long int maxdeg;
  long double minwt,sum=0,meanwt;
  long int i, j=0;
  long double *sum_wt;// sum of adjacent nodes
  igraph_vector_t neis;
  igraph_neimode_t omode;
  //igraph_bool_t isdirected=FALSE;
  
  //isdirected=igraph_is_directed(graph);
 /* if(igraph_is_directed(graph))
  	IGRAPH_ERROR("Directed graph is not allowed in weighted k-cores", IGRAPH_EINVAL);
  */
  if (mode != IGRAPH_ALL && mode != IGRAPH_OUT && mode != IGRAPH_IN) {
    IGRAPH_ERROR("Invalid mode in k-cores", IGRAPH_EINVAL);
  }
  if (!igraph_is_directed(graph) || mode==IGRAPH_ALL) {
    mode=omode=IGRAPH_ALL;
  } else if (mode==IGRAPH_IN) {
    omode=IGRAPH_OUT;
  } else {
    omode=IGRAPH_IN;
  }
	
//	mode=omode=IGRAPH_ALL;
  vert=igraph_Calloc(no_of_nodes, long int);
  if (vert==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vert);
  pos=igraph_Calloc(no_of_nodes, long int);
  if (pos==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, pos);

  /* maximum degree + degree of vertices */
  IGRAPH_CHECK(igraph_degree(graph, cores, igraph_vss_all(), mode, 
			     IGRAPH_LOOPS));
  maxdeg = (long int) igraph_vector_max(cores);

  //Rprintf("\nwt: %Lf",(long double) igraph_vector_max(weights));

  if (weights != 0){
	maxdeg = -1;

	if (igraph_vector_size(weights) != igraph_ecount(graph)) {
    		IGRAPH_ERROR("Invalid weights length", IGRAPH_EINVAL);
  	}

/*	for(i=0;i<no_of_edges;i++){
		sum = sum + (long double)VECTOR(*weights)[i];	    
	}

		  
	meanwt = sum/no_of_edges;
	for(i=0;i<no_of_edges;i++){
		
		VECTOR(*weights)[i] = ((long double)(VECTOR(*weights)[i]))/(long double)meanwt;	    
		//sum = sum/minwt;
		 //= roundl(sum);//confusion
	}
*/	minwt =(long double) igraph_vector_min(weights);
	for(i=0;i<no_of_edges;i++){
		
		sum = ((long double)(VECTOR(*weights)[i]))/(long double)minwt;	    
		//sum = sum/minwt;
		VECTOR(*weights)[i] = roundl(sum);//confusion
	}
	sum_wt=igraph_Calloc(no_of_nodes, long double);
	if (sum_wt==0) {
    		IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  	}
  	IGRAPH_FINALLY(igraph_free, sum_wt);
	igraph_integer_t edgefrom, edgeto;
	for(i=0; i<no_of_edges; i++){
		igraph_edge(graph, (igraph_integer_t) i, &edgefrom, &edgeto);
		
		if(mode==IGRAPH_ALL){
			sum_wt[edgefrom] += (long double) VECTOR(*weights)[i];
			sum_wt[edgeto] += (long double) VECTOR(*weights)[i];
		}
		else if(mode==IGRAPH_IN)
			sum_wt[edgeto] += (long double) VECTOR(*weights)[i];
		else if(mode==IGRAPH_OUT)
			sum_wt[edgefrom] += (long double) VECTOR(*weights)[i];				
	}
//	for(i=0; i<no_of_nodes; i++)

//		Rprintf("\nsumwt:%Lf",sum_wt[i]);

	long double temp;
	long int k;        
	if(type==1){
		for(i=0; i<no_of_nodes; i++){
			temp = (long double)VECTOR(*cores)[i];
                	temp = pow(temp,a);	
			sum_wt[i] = pow(sum_wt[i],b);
			temp = temp*sum_wt[i];
			//long double c = a + b;
			temp = pow(temp, (1/(a+b)));
			k = roundl(temp);
          //    	  Rprintf("\nK:%li",k);  
			VECTOR(*cores)[i] = k;
	//      	  Rprintf("\ncore:%li",(long int)VECTOR(*cores)[i]); 
			if(maxdeg < k){
				maxdeg=k;
			//	Rprintf("\nmax degree: %li",maxdeg);
			}
		}
	  } 
	if(type==2){ 
			for(i=0; i<no_of_nodes; i++){
			temp = (long double)VECTOR(*cores)[i];
                	temp = temp * l;	
			sum_wt[i] = sum_wt[i] * (1-l);
			temp = temp + sum_wt[i];			
			//long double c = a + b;
			k = temp;
              	  //Rprintf("\nK:%li",k);  
			VECTOR(*cores)[i] = k;
	      	//  Rprintf("\ncore:%li",(long int)VECTOR(*cores)[i]); 
			if(maxdeg < k){
				maxdeg=k;
			//	Rprintf("\nmax degree: %li",maxdeg);
			}
		}
	  } 
		
   }
 // Rprintf("\nmax degree: %li",maxdeg);
  bin=igraph_Calloc(maxdeg+1, long int);
  if (bin==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, bin);

  /* degree histogram */
  for (i=0; i<no_of_nodes; i++) {
    bin[ (long int)VECTOR(*cores)[i] ] += 1;
  }
  
  /* start pointers */
  j=0;
  for (i=0; i<=maxdeg; i++) {
    long int k=bin[i];
    bin[i] = j;
    j += k;
  }
  
  /* sort in vert (and corrupt bin) */
  for (i=0; i<no_of_nodes; i++) {
    pos[i] = bin[(long int)VECTOR(*cores)[i]];
    vert[pos[i]] = i;
    bin[(long int)VECTOR(*cores)[i]] += 1;
  }
  
  /* correct bin */
  for (i=maxdeg; i>0; i--) {
    bin[i] = bin[i-1];
  }
  bin[0]=0;

 // for (i=0; i<=maxdeg; i++) {
   // Rprintf("\nbin before loop %li is :%li",i,bin[i]);
 // }
  /* this is the main algorithm */
  IGRAPH_VECTOR_INIT_FINALLY(&neis, maxdeg);
  for (i=0; i<no_of_nodes; i++) {
    long int v=vert[i];
   // Rprintf("\nv is :%li",v);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, omode));
    for (j=0; j<igraph_vector_size(&neis); j++) {
      long int u=(long int) VECTOR(neis)[j];
      if (VECTOR(*cores)[u] > VECTOR(*cores)[v]) {
	long int du=(long int) VECTOR(*cores)[u];
	long int pu=pos[u];
	long int pw=bin[du];
	long int w=vert[pw];
	if (u != w) {
	  pos[u]=pw;
	  pos[w]=pu;
	  vert[pu]=w;
	  vert[pw]=u;
	}
	bin[du] += 1;
//	Rprintf("\nnew core of %li is :%li",u,(long int)VECTOR(*cores)[u]);
	VECTOR(*cores)[u] -= 1;
//	Rprintf("\nnew core of %li is :%li",u,(long int)VECTOR(*cores)[u]);
      }
    }
  }

//  for (i=0; i<=maxdeg; i++) {
 //   Rprintf("\nbin %li is :%li",i,bin[i]);
 // }
  
/*  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Rprintf ( "\nTime after completion: %s", asctime (timeinfo) );
*/

  end_time = clock();
  float time_diff = ((float) (end_time - start_time) / 1000000.0F ) * 1000; // CPU time taken in milliseconds
  Rprintf( "%f", time_diff );	
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);

  igraph_free(bin);
  igraph_free(pos);
  igraph_free(vert);
  igraph_free(sum_wt);
  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}
