/****s* FHI-aims/metis_nearest_wrapper
 * NAME
 *   metis_nearest_wrapper
 * SYNOPSIS
 *   void metis_nearest_wrapper(int *xadj,int *adjncy,int *grid_partition,int *n_grid_points, \
 *			   int *n_grid_batches)
 * PURPOSE
 *   Partition a graph that contains the nearest neighbour information of the grid points.
 * INPUTS
 *   o xadj -- row pointers for the adjacency graph
 *   o adjncy -- the adjacency graph
 *   o n_grid_points -- number of points in the grid
 *   o n_grid_bathces -- desired number of desired
 * OUTPUT
 *   o grid_partition -- partitioning of the grid
 * AUTHOR
 *   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
 * SEE ALSO
 *   Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
 *   Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
 *   "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
 *   Computer Physics Communications (2008), submitted.
 * COPYRIGHT
 *   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
 *   e.V. Please note that any use of the "FHI-aims-Software" is subject to
 *   the terms and conditions of the respective license agreement."
 * HISTORY
 *   Release version, FHI-aims (2008).
 * SOURCE
 */
#include <stdio.h>
#ifdef Add_
void metis_nearest_wrapper_(int *xadj,int *adjncy,int *grid_partition,int *n_grid_points, \
			    int *n_grid_batches)
#else
     void metis_nearest_wrapper(int *xadj,int *adjncy,int *grid_partition,int *n_grid_points, \
			   int *n_grid_batches)
#endif
{

  int nparts;
  int *vwgt;
  int *adjwgt;
  int wgtflag;
  int numflag;
  int edgecut;
  int options[] = {0,0,0,0,0};

  /* initialize parameters for METIS */
  vwgt = NULL;
  adjwgt = NULL;
  wgtflag = 0;
  nparts = *n_grid_batches;
  numflag = 1;
  
  METIS_PartGraphKway(n_grid_points, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, \
		      &nparts, options, &edgecut, grid_partition);

  return;
}

/*******/
