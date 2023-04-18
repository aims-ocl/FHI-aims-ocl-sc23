#ifndef _AIMS_ATOM_CENTERED_GRID_BUNDLE_H_
#define _AIMS_ATOM_CENTERED_GRID_BUNDLE_H_

#include <vector>
#include <map>

#include <aimsAtomCenteredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>

using namespace std;

// re-use grids which have been triangulated earlier
#define USE_GRID_MAP

class aimsAtomCenteredGridBundle
{
  public:
    aimsAtomCenteredGridBundle();
    aimsAtomCenteredGridBundle(int _mpiAtom0, int _mpiAtom1, int _nAtoms);
    void setCurrentAtomIndex(int idx);
    int getCurrentAtomIndex();
    const double *getCenter();
    void setCenter(double *xyz); // set the center for the current atom
    int getAtomicNumber();
    void setAtomicNumber(int z); // set the atomic number for the current atom
    void addPoint(double *xyz, double *dat); // add point data to the current atom
    void addPoint(vector <double> &singlePointData); // add point data to the current atom
    void calculateGrids() { calculateGrids(false); };
    void calculateGrids(bool verbose); // perform the actual sorting and triangulation of the grids
    aimsAtomCenteredGrid* getCurrentAtom();
    void freeAtom();
    map<int, aimsAtomCenteredGrid> atoms;
    void setArrayNames(vector<string> names);

  private:
    map<int, vtkSmartPointer<vtkCellArray> > grids;
    int currentAtom; // index of the current atom
    int mpiAtom0; // lower and
    int mpiAtom1; // upper indices of the atoms
    int nAtoms; // total number of atoms
    vector<string> arrayNames;
};

#endif
