//****s* FHI-aims/metis_tetgen_wrapper
// NAME
//   metis_tetgen_wrapper
// SYNOPSIS
//   void metis_tetgen_wrapper(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
//  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
// PURPOSE
//   Build a Delaunay triangulation of the grid points using TetGen and partition the graph with Metis.
// INPUTS
//   o grid_points_x -- x-coordinates of the grid points
//   o grid_points_y -- y-coordinates of the grid points
//   o grid_points_z -- z-coordinates of the grid points
//   o n_grid_bathces -- desired number of batches
//   o n_grid_points -- number of points in the grid
// OUTPUT
//   o grid_partition -- partitioning of the grid
// AUTHOR
//   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
// SEE ALSO
//   Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
//   Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
//   "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
//   Computer Physics Communications (2008), submitted.
// COPYRIGHT
//   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
//   e.V. Please note that any use of the "FHI-aims-Software" is subject to
//   the terms and conditions of the respective license agreement."
// HISTORY
//   Release version, FHI-aims (2008).
// SOURCE
#include <iostream>
#include "tetgen.h"
using namespace std;

#ifdef Add_
extern "C" void metis_tetgen_wrapper_(double *, double *,double *, int *, int *, int *);
#else
extern "C" void metis_tetgen_wrapper(double *, double *,double *, int *, int *, int *);
#endif
extern "C" void METIS_PartMeshNodal(int *, int *, int *, int *, int *, int *, int *, int*, int *);

#ifdef Add_
void metis_tetgen_wrapper_(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
#else
void metis_tetgen_wrapper(double *grid_points_x, double *grid_points_y, double *grid_points_z, \
  		    int *grid_partition, int *n_grid_points, int *n_grid_batches)
#endif
{
  tetgenio in, out;

  int i_point;
  int i;
  int etype;
  int numflag;
  int edgecut;
  int *epart;
  int nparts;

  in.numberofpoints = *n_grid_points;
  in.pointlist = new REAL[in.numberofpoints * 3];
  in.firstnumber = 1;

  //  cout << in.numberofpoints << endl;

  i = 0;
  for (i_point=0; i_point < in.numberofpoints; i_point++)
    {
      in.pointlist[i++] = grid_points_x[i_point];
      in.pointlist[i++] = grid_points_y[i_point];
      in.pointlist[i++] = grid_points_z[i_point];
    }

  //  tetrahedralize("Q", &in, &out);
  tetrahedralize("Q", &in, &out);
  //  out.save_nodes("testout");
  //out.save_elements("testout");
  //out.save_faces("testout"); 

  //  cout << out.numberofpoints << endl;
  *n_grid_points = out.numberofpoints;

  epart = new int[out.numberoftetrahedra];
  etype = 2;
  nparts = *n_grid_batches;
  numflag = 1;

  METIS_PartMeshNodal(&out.numberoftetrahedra, &out.numberofpoints, out.tetrahedronlist, \
		      &etype, &numflag, &nparts, &edgecut, epart, grid_partition);

  //  for (i_point=0; i_point < out.numberofpoints; i_point++)
  //  {
  //    cout << grid_partition[i_point] << endl;
  //  }

  // delete[] in.pointlist;
  delete[] epart;

  return;

}
/******/
