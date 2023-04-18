#ifndef _AIMS_RECTILINEAR_GRID_H_
#define _AIMS_RECTILINEAR_GRID_H_

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>

#include <aimsAtomCenteredGrid.h>

#include <vector>

class aimsRectilinearGrid
{
  private:
    double *cartX;
    double *cartY;
    double *cartZ;
    int    nx;
    int    ny;
    int    nz;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    bool verbose;
    bool debug;

    vtkSmartPointer<vtkRectilinearGrid> rectGrid;
    vector< vtkSmartPointer<vtkDoubleArray> > vtkDataArrays;

    void initialize(
        int nX, int nY, int nZ,
        double xMin, double xMax,
        double yMin, double yMax,
        double zMin, double zMax);

    void destroy();


  public:
    aimsRectilinearGrid();
    aimsRectilinearGrid(
        int &nX, int &nY, int &nZ,
        double &xMin, double &xMax,
        double &yMin, double &yMax,
        double &zMin, double &zMax);
    ~aimsRectilinearGrid();

    void initializeDataArrays(vector<string> dataNames);

    void addFromAimsAtomCenteredGrid(aimsAtomCenteredGrid *atom);
    double getDataValue(ssize_t &dataIdx, vtkIdType &index);
    void setDataArray(double *dataArray, vtkIdType dataSize, string dataName, ssize_t dataIdx);

    void dumpVtr(string &fileName);
    void dumpVti(string &fileName);

    void setVerbosity(bool _verbose, bool _debug);
};

#endif
