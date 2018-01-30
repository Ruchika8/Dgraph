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

int dgraph_wtcoreness2(const igraph_t *graph, igraph_vector_t *cores, igraph_neimode_t mode, const igraph_vector_t* weights, 
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
  long int *bin, *vert, *pos, *tempbin,*visited,*neivisited;
  long double a=alpha, b=beta, l=lambda;
  long int maxdeg;
  long double minwt,sum=0,meanwt;
  long int i, j=0, z;
  long double tempdeg, tempwt;
  long int k;
  long double *sum_wt;// sum of adjacent nodes
  igraph_vector_t current_degree; //current degree of nodes
  igraph_vector_t neis;
  igraph_neimode_t omode;
  igraph_integer_t edgefrom, edgeto;
  IGRAPH_VECTOR_INIT_FINALLY(&current_degree, no_of_nodes);  

//igraph_bool_t isdirected=FALSE;
  
  //isdirected=igraph_is_directed(graph);
/*  if(igraph_is_directed(graph))
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
	
//  mode=omode=IGRAPH_ALL;
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
  visited=igraph_Calloc(no_of_nodes, long int);
  if (visited==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, visited);
  neivisited=igraph_Calloc(no_of_nodes, long int);
  if (neivisited==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, neivisited);

  //current_degree=igraph_Calloc(no_of_nodes, igraph_vector_t);
 // if (current_degree==0) {
  //  IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  //}	
  //IGRAPH_FINALLY(igraph_free, current_degree);
  
  IGRAPH_CHECK(igraph_degree(graph, &current_degree, igraph_vss_all(), mode, 
			     IGRAPH_LOOPS));
  /* maximum degree + degree of vertices */
  IGRAPH_CHECK(igraph_degree(graph, cores, igraph_vss_all(), mode, 
			     IGRAPH_LOOPS));
  maxdeg = (long int) igraph_vector_max(cores);
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
		
		VECTOR(*weights)[i] = ((long double)(VECTOR(*weights)[i]))/meanwt;	    
		//sum = sum/minwt;
		 //= roundl(sum);//confusion
	}
*/	minwt =(long double) igraph_vector_min(weights);
	for(i=0; i<no_of_edges; i++){
		
		sum = ((long double)(VECTOR(*weights)[i]))/minwt;	    
		//sum = sum/minwt;
		VECTOR(*weights)[i] = roundl(sum);
	}
	sum_wt=igraph_Calloc(no_of_nodes, long double);
	if (sum_wt==0) {
    		IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  	}
  	IGRAPH_FINALLY(igraph_free, sum_wt);
	
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
	long double temp;
	long int k;        
	if(type==1){
		for(i=0; i<no_of_nodes; i++){
			temp = (long double)VECTOR(*cores)[i];
                	temp = pow(temp,a);	
			tempwt = pow(sum_wt[i],b);
			temp = temp * tempwt;
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
			tempwt = sum_wt[i] * (1-l);
			temp = temp + tempwt;
			//long double c = a + b;
			k = temp;
          //    	  Rprintf("\nK:%li",k);  
			VECTOR(*cores)[i] = k;
	//      	  Rprintf("\ncore:%li",(long int)VECTOR(*cores)[i]); 
			if(maxdeg < k){
				maxdeg=k;
			//	Rprintf("\nmax degree: %li",maxdeg);
			}
		}
	  } 
  }
  
  bin=igraph_Calloc(maxdeg+1, long int);
  if (bin==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, bin);

  tempbin=igraph_Calloc(maxdeg+2, long int);
  if (tempbin==0) {
    IGRAPH_ERROR("Cannot calculate k-cores", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tempbin);

  /* degree histogram */
  for (i=0; i<no_of_nodes; i++) {
    bin[ (long int)VECTOR(*cores)[i] ] += 1;
  }
  
  /* start pointers */
  long int start=0;
  for (i=0; i<=maxdeg; i++) {
    long int num=bin[i];
    bin[i] = start;    
    start += num;
    if(num!=0)
    	tempbin[i+1]=1;
  }
  
  /* sort in vert (and corrupt bin) */
  for (i=0; i<no_of_nodes; i++) {
    visited[i]=0;
    pos[i] = bin[(long int)VECTOR(*cores)[i]];
    vert[pos[i]] = i;
    bin[(long int)VECTOR(*cores)[i]] += 1;
  }
  
  /* correct bin */
  for (i=maxdeg; i>0; i--) {
    bin[i] = bin[i-1];
  }
  bin[0]=0;
/*Rprintf("\n");
for (i=0; i<maxdeg; i++) {
	Rprintf("	bin of %li:%li ",i,bin[i]);
}

Rprintf("\n");
for (i=0; i<maxdeg; i++) {
	Rprintf("	tempbin of %li:%li ",i,tempbin[i]);
}
for (i=0; i<no_of_nodes; i++) {
    Rprintf("	pos of %li:%li ",i+1,pos[i]+1);
    Rprintf("	vert of %li:%li ",i+1,vert[i]+1);
 }
  */

  /* this is the main algorithm */
  long int min_core = (long int) igraph_vector_min(cores);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, maxdeg);
  i=0;
 // long int maxnode=0;
  while (i<no_of_nodes) {
	
    long int v=vert[i];
    Rprintf("	in while v: %li",v);
    if(visited[v]==1){
		Rprintf("in visited");
		i=0;
		while(visited[vert[i]]!=0 && i<no_of_nodes)
			i++;
		continue;		
    }
    Rprintf("\n after while");
    visited[v]=1;
    //maxnode++;
	for(long int c=0; c<no_of_nodes; c++){
		neivisited[c]=0;
	}
	      Rprintf("\n core of %li is: %li ",v+1,(long int) VECTOR(*cores)[v]);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, omode));
    for (j=0; j<igraph_vector_size(&neis); j++) {
      long int u=(long int) VECTOR(neis)[j];
	Rprintf("u:%li",u);
      if (VECTOR(*cores)[u] > VECTOR(*cores)[v] && visited[u]!=1 && neivisited[u]!=1) {
	neivisited[u]=1;
	VECTOR(current_degree)[u] -= 1;
	for(z=0; z<no_of_edges; z++){
		igraph_edge(graph, (igraph_integer_t) z, &edgefrom, &edgeto);
		if(mode==IGRAPH_ALL){
			
			if((edgefrom==v && edgeto==u) || (edgefrom==u && edgeto==v))		
			{//Rprintf("\nin");
				sum_wt[u] -= (long double) VECTOR(*weights)[z];}		}
		else if(mode==IGRAPH_IN){
			if(edgefrom==v && edgeto==u)
				sum_wt[u] -= (long double) VECTOR(*weights)[z];
		}
		else if (mode==IGRAPH_OUT){
			if(edgefrom==u && edgeto==v)
				sum_wt[u] -= (long double) VECTOR(*weights)[z];
		}				
	}
	Rprintf("\nsumwt of %li :%Lf",u+1,sum_wt[u]);
	tempdeg = (long double)VECTOR(current_degree)[u];
	if(type==1){
		tempdeg = pow(tempdeg,a);
		tempwt = pow(sum_wt[u],b);
		tempdeg = tempdeg * tempwt;
		tempdeg = pow(tempdeg, (1/(a+b)));
		k = roundl(tempdeg);
	}
	if(type==2){
		//Rprintf("\nsumwt:%Lf",sum_wt[u]);
		//Rprintf("\ndeg:%Lf",tempdeg);

		tempdeg = tempdeg * l;
		tempwt = sum_wt[u] * (1-l);
		tempdeg = tempdeg + tempwt;
		k = tempdeg;
        }

	Rprintf("\nnew core:%li",k);
	long int du=(long int) VECTOR(*cores)[u];
	Rprintf("     old core:%li",du);
	//long int difference = du - k;
       // Rprintf("\ndiff:%li",difference);
//	if (difference < min_core)
//		difference = min_core;		
//for(long int it=0; it<difference; it++){
	long int count=0;
	long int it=0;
	for(it=du;it>k;it--){
		if(tempbin[it]==1)
			count++;
	}
	
       du=(long int) VECTOR(*cores)[u];	
       long int tempdu=du, x=du, binval=bin[du];
       Rprintf("\ncount :%li",count);
	for(it=0; it<=count; it++){
//		while(bin[(long int) VECTOR(*cores)[u]]!=bin[k]){
//		if(bin[(long int) VECTOR(*cores)[u]]!=bin[k]){
		long int pu=pos[u];
		Rprintf("	pu:%li",pu);
		long int pw=bin[x];
		long int w=vert[pw];
		if (u != w) {
	  		pos[u]=pw;
	  		pos[w]=pu;
	  		vert[pu]=w;
	  		vert[pw]=u;
		}
		while(bin[x]==bin[x-1] && x!=k){
			bin[x]++;
			x--;
		}
		if(x!=k)
			bin[x]++;
		x--;

	}
	if(bin[tempdu]==bin[tempdu+1] && tempdu<maxdeg){
		tempbin[tempdu+1]=0;	
	}	
	tempdu++;
	while(binval==bin[tempdu] && tempdu<=maxdeg){
		bin[tempdu]++;
		tempdu++;	
	}
	//tempbin[k]=1;
	if(k<maxdeg)
        	tempbin[k+1]=1;
        VECTOR(*cores)[u]=k;
/*

	for(it=0;it<=no_of_nodes;it++)
	{
		Rprintf("	Vert %li",vert[it]);
	}
Rprintf("\n");
	for(it=0;it<=no_of_nodes;it++)
	{
		Rprintf("	Pos %li",pos[it]);
	}
Rprintf("\n");
	for (it=0; it<maxdeg; it++) {
	Rprintf("	bin of %li:%li ",it,bin[it]);
}

Rprintf("\n");
for (it=0; it<maxdeg; it++) {
	Rprintf("	tempbin of %li:%li ",it,tempbin[it]);
}
*/

      }//endif(VECTOR(*cores)[u] > VECTOR(*cores)[v])
    }//endforj
   i++;
  }//endfori
	Rprintf("after for");
/*
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Rprintf ( "\nTime after completion: %s", asctime (timeinfo) );
  */

  end_time = clock();
  float time_diff = ((float) (end_time - start_time) / 1000000.0F ) * 1000; // CPU time taken in milliseconds
  Rprintf( "\n %f", time_diff );

  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(1);
Rprintf( "\n first");
  igraph_vector_destroy(&current_degree);
  IGRAPH_FINALLY_CLEAN(1);
Rprintf( "\n second" );
  
  igraph_free(bin);
Rprintf( "\n thirdd" );
  igraph_free(pos);
Rprintf( "\n fifth" );
  igraph_free(vert);
Rprintf( "\n sixth" );
  igraph_free(tempbin);
Rprintf( "\n fourth" );
  igraph_free(visited);
Rprintf( "\n seventh" );
  igraph_free(sum_wt);
Rprintf( "\n eight" );
  igraph_free(neivisited);
  IGRAPH_FINALLY_CLEAN(7);
  return 0;
}

