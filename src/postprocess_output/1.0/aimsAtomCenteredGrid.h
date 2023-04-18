#ifndef _AIMS_ATOM_CENTERED_GRID_H_
#define _AIMS_ATOM_CENTERED_GRID_H_

#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vector>

// Use <libqhullvtk>, which wraps a qhull delaunay triangulation for the use with VTK.
// It is not fully functional, yet, and therefore disabled by default.
#ifdef HAVE_QHULLVTK
#include <qhullvtk.h>
#endif

// Use the VTK polydata cleaner before feeding into
// the VTK Delaunay triangulation.  Should not be used.
//#define USE_CLEANER

using namespace std;

class aimsAtomCenteredGrid
{
  public:
    int    atomicNumber;
    bool   atomicNumberSet;
    double nucleusCoordinate[3];
    bool   nucleusCoordinateSet;

    vtkSmartPointer<vtkPoints>           points;
    vtkSmartPointer<vtkCellArray>        vertices;
    vtkSmartPointer<vtkPolyData>         polydata;
    vtkSmartPointer<vtkDelaunay3D>       delaunay3D;
    vtkSmartPointer<vtkUnstructuredGrid> ugrid;

    vtkSmartPointer<vtkDoubleArray>      scalar;
    vector< vtkSmartPointer<vtkDoubleArray> > vtkDataArrays;

    //void init();
    aimsAtomCenteredGrid();
    void setCenter(double *xyz);
    void setAtomicNumber(int z);
    const double *getCenter();
    int getAtomicNumber();
    int getNumberOfGridPoints();

    void addPoint(double *xyz, double *dat);
    void addPoint(vector <double> &singlePointData);

    void sortPoints();

    // connect points, vertices, scalar data to a polydata object
    void preparePolyData();
    void triangulate();
    void dumpPolyData(int atomId);
    void dumpUGrid(int atomId);

    vector<double> query(double *point, vtkSmartPointer<vtkCell> *previousCell, vtkIdType *previousCellId);

    vtkSmartPointer<vtkCellArray> getTriangulation();
    void setTriangulation(vtkSmartPointer<vtkCellArray> cells);
    vtkSmartPointer<vtkPoints> getPoints();
    void setPoints(vtkSmartPointer<vtkPoints> points);
    void setArrayNames(vector<string> names);

  private:
    vector< vector<double> > gridPointData;
    vector<string> arrayNames;
    void triangulateVtk();

#ifdef HAVE_QHULLVTK
    void triangulateLibQhull();
#endif

    friend bool pointComparator(const vector<double> &pointDataA, const vector<double> &pointDataB);
};
// helper function for the spatial sorting of points
bool pointComparator(const vector<double> &coordA, const vector<double> &coordB);

#endif // _AIMS_ATOM_CENTERED_GRID_H_
