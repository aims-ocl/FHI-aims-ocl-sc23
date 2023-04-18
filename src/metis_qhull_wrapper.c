/****s* FHI-aims/metis_qhull_wrapper
 * NAME
 *   metis_qhull_wrapper
 * SYNOPSIS
 *   void metis_qhull_wrapper(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
 *  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
 * PURPOSE
 *   Build a Delaunay triangulation of the grid points using qhull and partition the graph with Metis.
 * INPUTS
 *   o grid_points_x -- x-coordinates of the grid points
 *   o grid_points_y -- y-coordinates of the grid points
 *   o grid_points_z -- z-coordinates of the grid points
 *   o n_grid_bathces -- desired number of batches
 *   o n_grid_points -- number of points in the grid
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
#include "qhull_a.h"

#ifdef Add_
void metis_qhull_wrapper_(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
#else
void metis_qhull_wrapper(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
#endif
{

  int i_point;
  int i;
  int etype;
  int numtet;
  int numflag;
  int edgecut;
  int *epart;
  int nparts;
  int *tetlist;
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
  numpoints = *n_grid_points;

  points = malloc( dim*numpoints*sizeof(coordT) );

  i = 0;
  for (i_point=0; i_point < numpoints; i_point++)
    {
      points[i++] = grid_points_x[i_point];
      points[i++] = grid_points_y[i_point];
      points[i++] = grid_points_z[i_point];
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
  epart = malloc( numtet*sizeof(int) );
  etype = 2;
  nparts = *n_grid_batches;
  numflag = 0;
  
  METIS_PartMeshNodal(&numtet, &numpoints, tetlist, \
		      &etype, &numflag, &nparts, &edgecut, epart, grid_partition);
  
  for (i_point=0; i_point < numpoints; i_point++)
    {
      grid_partition[i_point] = grid_partition[i_point] + 1;
    }

  free(points);
  free(epart);
  free(tetlist);

  return;

}

/******/
