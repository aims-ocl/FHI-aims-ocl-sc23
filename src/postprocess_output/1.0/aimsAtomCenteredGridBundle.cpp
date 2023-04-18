#include <assert.h>
#include <aimsAtomCenteredGridBundle.h>

aimsAtomCenteredGridBundle::aimsAtomCenteredGridBundle()
{
  aimsAtomCenteredGridBundle(-1, -1, -1);
}

aimsAtomCenteredGridBundle::aimsAtomCenteredGridBundle(int _mpiAtom0, int _mpiAtom1, int _nAtoms)
{
  currentAtom = -1;
  mpiAtom0 = _mpiAtom0;
  mpiAtom1 = _mpiAtom1;
  nAtoms = _nAtoms;
  //
  for (int i=0; i<nAtoms; ++i)
  {
    atoms.insert( pair<int,aimsAtomCenteredGrid>(i, aimsAtomCenteredGrid() ) );
  }
}

int aimsAtomCenteredGridBundle::getCurrentAtomIndex()
{
  return currentAtom;
}

void aimsAtomCenteredGridBundle::setCurrentAtomIndex(int idx)
{
  assert( (idx >= 0)&&((idx < nAtoms)) );
  currentAtom = idx;
}

const double *aimsAtomCenteredGridBundle::getCenter()
{
  return atoms[currentAtom].getCenter();
}

void aimsAtomCenteredGridBundle::setCenter(double *xyz)
{
  atoms[currentAtom].setCenter(xyz);
}

int aimsAtomCenteredGridBundle::getAtomicNumber()
{
  return atoms[currentAtom].getAtomicNumber();
}

void aimsAtomCenteredGridBundle::setAtomicNumber(int z)
{
  atoms[currentAtom].setAtomicNumber(z);
}

void aimsAtomCenteredGridBundle::addPoint(double *xyz, double *dat)
{
  assert( (currentAtom >= mpiAtom0)&&((currentAtom <= mpiAtom1)) );
  atoms[currentAtom].addPoint(xyz, dat);
}

void aimsAtomCenteredGridBundle::addPoint(vector <double> &singlePointData)
{
  assert( (currentAtom >= mpiAtom0)&&((currentAtom <= mpiAtom1)) );
  atoms[currentAtom].addPoint(singlePointData);
}

void aimsAtomCenteredGridBundle::calculateGrids(bool verbose)
{
  for (int i=mpiAtom0; i<=mpiAtom1; ++i)
  {
    atoms[i].setArrayNames(arrayNames);
    atoms[i].sortPoints();
    atoms[i].preparePolyData();
#ifdef USE_GRID_MAP
    int z;
    z = atoms[i].getAtomicNumber();
    map<int, vtkSmartPointer<vtkCellArray> >::iterator it;
    it = grids.find(z);
    //
    if (it == grids.end()) // grid has not been encountered yet
    {
      if (verbose)
        cout << "Triangulating grid for atom #" << i << " ..." << endl;
      atoms[i].triangulate();
      grids[z] = atoms[i].getTriangulation();
    }
    else // grid is already known
    {
      if (verbose)
        cout << "Reusing known grid for atom #" << i << " ..." << endl;
      atoms[i].setTriangulation( grids[z] );
    }
#else
    if (verbose)
      cout << "Triangulating grid for atom #" << i << " ..." << endl;
    atoms[i].triangulate();
#endif
  }
}

aimsAtomCenteredGrid* aimsAtomCenteredGridBundle::getCurrentAtom()
{
  assert( (currentAtom >= mpiAtom0)&&((currentAtom <= mpiAtom1)) );
  aimsAtomCenteredGrid* atomPointer;
  atomPointer = &(atoms[currentAtom]);
  return atomPointer;
}

void aimsAtomCenteredGridBundle::freeAtom()
{
  atoms.erase(currentAtom);
}

void aimsAtomCenteredGridBundle::setArrayNames(vector<string> names)
{
  arrayNames = names;
}
