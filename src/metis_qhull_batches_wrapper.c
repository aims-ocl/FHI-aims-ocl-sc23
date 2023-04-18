/****s* FHI-aims/metis_qhull_batches_wrapper
 * NAME
 *   metis_qhull_batches_wrapper
 * SYNOPSIS
 *   void metis_qhull_batches_wrapper(double *coms_x, double *coms_y, double *coms_z, \
 *			 int* batch_weights, int *batch_partition, int *n_grid_batches, int* n_tasks)
 * PURPOSE
 *   Distribute the grid batches to MPI-tasks by first buiding a Delaunay triangulation of the
 *   center of mass points of the batches and then partitioning the graph with Metis.
 * INPUTS
 *   o coms_x -- x-coordinates of the centers of mass for the grid batches
 *   o coms_y -- y-coordinates of the centers of mass for the grid batches
 *   o coms_z -- z-coordinates of the centers of mass for the grid batches
 *   o batch_weights -- weights of each of the batches (i.e. size of the batch)
 *   o n_grid_batches -- number of batches in the grid
 *   o n_tasks -- number of tasks to distribute the batches over
 * OUTPUT
 *   o batch_partition -- partitioning of the batches to n_tasks bins
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
#include "qhull_a.h"

#ifdef Add_
void metis_qhull_batches_wrapper_(double *coms_x, double *coms_y, double *coms_z, \
				  int* batch_weights, int *batch_partition, \
				  int *n_grid_batches, int *n_tasks)
#else
void metis_qhull_batches_wrapper(double *coms_x, double *coms_y, double *coms_z, \
				 int* batch_weights, int *batch_partition, \
				 int *n_grid_batches, int* n_tasks)
#endif
{

  int i_point;
  int i;
  int etype;
  int numtet;
  int numflag;
  int edgecut;
  int nparts;
  int *tetlist;
  int *nxadj;
  int *nadjncy;
  int *vwgt;
  int *adjwgt;
  int wgtflag;
  int options[5];

  vertexT *vertex, **vertexp;
  setT *vertices;

  int dim;	            /* dimension of points */
  int numpoints;            /* number of points */
  coordT *points;           /* array of coordinates for each point */
  boolT ismalloc;           /* True if qhull should free points in qh_freeqhull() or reallocation */
  char flags[]= "qhull d Qt Qbb"; /* option flags for qhull, see qh_opt.htm */
  FILE *outfile= stdout;    /* output from qh_produce_output()
			       use NULL to skip qh_produce_output() */
  FILE *errfile= stderr;    /* error messages from qhull code */
  int exitcode;             /* 0 if no error from qhull */
  facetT *facet;	    /* set by FORALLfacets */
  int curlong, totlong;	    /* memory remaining after qh_memfreeshort */

  /* initialize dim, numpoints, points[], ismalloc here */
  ismalloc = False;
  dim = 3;
  numpoints = *n_grid_batches;

  points = malloc( dim*numpoints*sizeof(coordT) );

  i = 0;
  for (i_point=0; i_point < numpoints; i_point++)
    {
      points[i++] = coms_x[i_point];
      points[i++] = coms_y[i_point];
      points[i++] = coms_z[i_point];
    }

  exitcode= qh_new_qhull (dim, numpoints, points, ismalloc, \
			  flags, outfile, errfile); 

  /* count the resulting tetras in the Delaunay triangulation */
  numtet = 0;
  if (!exitcode){
    FORALLfacets {
      if (!facet->upperdelaunay) {
	numtet = numtet + 1;
      }
    }
  }

  /* load the tetras to a structure for METIS */
  tetlist = malloc( 4*numtet*sizeof(int) );
  i = 0;
  FORALLfacets {
    if (!facet->upperdelaunay) {
      FOREACHvertex_(facet->vertices) {
	tetlist[i++] = qh_pointid(vertex->point);
      }
    }
  }

  /* free qhull memory */
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort (&curlong, &totlong);
  if (curlong || totlong) 
    fprintf (errfile, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
  
  /* initialize parameters for METIS */
  etype = 2;
  nparts = *n_tasks;
  numflag = 0;

  nxadj = malloc( (numpoints+1)*sizeof(int) );
  nadjncy = malloc( (20*numpoints)*sizeof(int) );

  /*  
  i = 0;
  printf("Before MeshToNodal %d %d \n", numtet, numpoints);
  for (i_point=0; i_point < numpoints; i_point++)
    {
      printf("%d %d %d %d \n", tetlist[i++], tetlist[i++], tetlist[i++], tetlist[i++]);
    } 
  */

  METIS_MeshToNodal( &numtet, &numpoints, tetlist, &etype, &numflag, \
		     nxadj, nadjncy );


  vwgt = malloc( numpoints*sizeof(int) );
  for (i_point=0; i_point < numpoints; i_point++)
    {
      vwgt[i_point] = batch_weights[i_point];
    }

  adjwgt = NULL;
  wgtflag = 2;
  options[0] = 0;
  
  METIS_PartGraphKway( &numpoints, nxadj, nadjncy, vwgt, adjwgt, &wgtflag, \
		       &numflag, &nparts, options, &edgecut, batch_partition );

  /*
  for (i_point=0; i_point < numpoints; i_point++)
    {
      batch_partition[i_point] = batch_partition[i_point] + 1;
    }
  */

  free(vwgt);
  free(nxadj);
  free(nadjncy);
  free(points);
  free(tetlist);

  return;

}

/******/ 
