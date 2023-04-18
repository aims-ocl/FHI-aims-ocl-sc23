#include <assert.h>

#include <iostream>
#include <string>
#include <sstream>

#include <aimsAtomCenteredGrid.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMath.h>

#include <algorithm>

using namespace std;

//void aimsAtomCenteredGrid::init()
aimsAtomCenteredGrid::aimsAtomCenteredGrid()
{
  points     = vtkSmartPointer<vtkPoints>::New();
  vertices   = vtkSmartPointer<vtkCellArray>::New();
  scalar     = vtkSmartPointer<vtkDoubleArray>::New();
  polydata   = vtkSmartPointer<vtkPolyData>::New();
  delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();
  ugrid      = vtkSmartPointer<vtkUnstructuredGrid>::New();
  //
  // set some default values
  atomicNumber = -1;
  atomicNumberSet=false;
  //
  for (int i=0; i<3; i++)
    nucleusCoordinate[i]=0.0;
  nucleusCoordinateSet=false;
  //
  scalar->SetName("pot");
}

void aimsAtomCenteredGrid::setCenter(double *xyz)
{
  for (int i=0; i<3; i++)
    nucleusCoordinate[i]=xyz[i];
  nucleusCoordinateSet=true;
}

void aimsAtomCenteredGrid::setAtomicNumber(int z)
{
  atomicNumber=z;
  atomicNumberSet=true;
}

const double * aimsAtomCenteredGrid::getCenter()
{
  return &nucleusCoordinate[0];
}

int aimsAtomCenteredGrid::getAtomicNumber()
{
  return atomicNumber;
}

int aimsAtomCenteredGrid::getNumberOfGridPoints()
{
  return (points->GetNumberOfPoints());
}


// add a point, specified by the point's coordinates x,y,z followed by ir, ia, followed by the scalar data
void aimsAtomCenteredGrid::addPoint(vector<double> &singlePointData)
{
  assert( singlePointData.size() > 5 );
  gridPointData.push_back(singlePointData);
}
// legacy wrapper, don't use this any more
void aimsAtomCenteredGrid::addPoint(double *xyz, double *dat)
{
  vector<double> singlePointData;
  for (int i=0; i<3; i++)
    singlePointData.push_back(xyz[i]);
  singlePointData.push_back(dat[0]);
  addPoint(singlePointData);
}


// sort the points
void aimsAtomCenteredGrid::sortPoints()
{
//  cout << "aimsAtomCenteredGrid::sortPoints()" << endl;
  sort(gridPointData.begin(), gridPointData.end(), pointComparator);
}
//
bool pointComparator(const vector<double> &pointDataA, const vector<double> &pointDataB)
{
  assert( pointDataA.size() >= 5 );
  assert( pointDataB.size() >= 5 );
  //
  if ( fabs(pointDataA[3] - pointDataB[3]) < 1.e-12 ) // both points are on the same radial shell
  {
    return (pointDataA[4] < pointDataB[4]); // decide based on the angular (Lebedev) index
  }
  else
  {
    return (pointDataA[3] < pointDataB[3]); // decide based on the radial index
  }
  //
}


// connect points, vertices, scalar data to a polydata object
void aimsAtomCenteredGrid::preparePolyData()
{
//  cout << "aimsAtomCenteredGrid::preparePolyData()" << endl;
  //
  // number of points in the data set
  ssize_t nPoints = gridPointData.size();
  assert( nPoints >= 1 );
  // number of data values per point -- the first 5 are x,y,z,ir,ia
  ssize_t nValues = (gridPointData[0].size() - 5);
  assert( nValues >= 0);
  //
  if ((ssize_t)arrayNames.size() != nValues)
  {
    arrayNames.clear();
    for (ssize_t i=0; i < nValues; i++ )
    {
      stringstream is;
      is << i;
      arrayNames.push_back( string("field-")+is.str() );
    }
  }
  //
  vtkDataArrays.clear();
  for (ssize_t i=0; i < nValues; i++ )
  {
    vtkDataArrays.push_back( vtkSmartPointer<vtkDoubleArray>::New() );
    vtkDataArrays.back()->SetName( arrayNames[i].c_str() );
  }
  //
  for (vector< vector<double> >::iterator it=gridPointData.begin(); it != gridPointData.end(); ++it)
  {
    vector<double> &pDat = *it;
    //
    assert( (ssize_t)pDat.size() == (nValues + 5) );
    //
    vtkIdType pid[1];
    pid[0] = points->InsertNextPoint(pDat[0], pDat[1], pDat[2]);
    vertices->InsertNextCell(1, pid);
    // pDat[3] and pDat[4] hold the radial and angular indices, data values start from 5

    //    scalar->InsertNextValue(pDat[5]);

    for (ssize_t i=0; i < (ssize_t)nValues; i++ )
    {
      vtkDataArrays[i]->InsertNextValue( pDat[i+5] );
    }
  }
  //
  gridPointData.clear();
  //
  polydata->SetPoints(points);
  polydata->SetVerts(vertices);
//  polydata->GetPointData()->AddArray(scalar);
  for (ssize_t i=0; i < (ssize_t)nValues; i++ )
  {
    polydata->GetPointData()->AddArray( vtkDataArrays[i] );
  }
}

void aimsAtomCenteredGrid::triangulate()
{
#ifdef HAVE_QHULLVTK
  triangulateLibQhull();
#else
  triangulateVtk();
#endif
}

// Triangulate with the VTK-built-in Delaunay triangulation
// is somewhat slower than the qHull triangulation below
void aimsAtomCenteredGrid::triangulateVtk()
{
  delaunay3D->SetTolerance(1.e-6);
#ifdef USE_CLEANER
  vtkSmartPointer<vtkCleanPolyData> cleaner =
      vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInput(polydata);
  cleaner->Update();
  delaunay3D->SetInputConnection(cleaner->GetOutputPort());
#else
  delaunay3D->SetInput(polydata);
#endif
  delaunay3D->Update();
  //
  ugrid = delaunay3D->GetOutput();
  //
//  cout << "atom::triangulateVtk(): " << ugrid->GetNumberOfPoints() << " points, "
//      << ugrid->GetNumberOfCells() << " tetrahedra" << endl;
}


#ifdef HAVE_QHULLVTK
// perform the Delaunay triangulation using functions from the QHull library
void aimsAtomCenteredGrid::triangulateLibQhull()
{
  ugrid->SetPoints( polydata->GetPoints() );
  qDelaunay3D( ugrid );
  ugrid->GetPointData()->PassData( polydata->GetPointData() );
//  cout << "atom::triangulateLibQhull(): " << ugrid->GetNumberOfPoints() << " points, "
//      << ugrid->GetNumberOfCells() << " tetrahedra" << endl;
}
#endif


vtkSmartPointer<vtkCellArray> aimsAtomCenteredGrid::getTriangulation()
{
//  cout << "aimsAtomCenteredGrid::getTriangulation()" << endl << flush;
  vtkSmartPointer<vtkCellArray> triangulation = ugrid->GetCells();
  return triangulation;
}

void aimsAtomCenteredGrid::setTriangulation(vtkSmartPointer<vtkCellArray> cells)
{
//  cout << "aimsAtomCenteredGrid::setTriangulation()" << endl << flush;
  ugrid->SetPoints( polydata->GetPoints() );
  ugrid->GetPointData()->PassData( polydata->GetPointData() );
  ugrid->SetCells(VTK_TETRA, cells);
  return;
}

vtkSmartPointer<vtkPoints> aimsAtomCenteredGrid::getPoints()
{
//  cout << "aimsAtomCenteredGrid::getPoints()" << endl << flush;
  vtkSmartPointer<vtkPoints> points = ugrid->GetPoints();
  return points;
}

void aimsAtomCenteredGrid::setPoints(vtkSmartPointer<vtkPoints> points)
{
//  cout << "aimsAtomCenteredGrid::setPoints()" << endl << flush;
  ugrid->SetPoints(points);
  return;
}

// dump the polydata object of an atom, i.e. the grid points (w/o grid)
void aimsAtomCenteredGrid::dumpPolyData(int atomId)
{
  string outFileName;
  stringstream ss;
  ss << atomId;
  outFileName.assign("atomPolyData");
  outFileName.append( ss.str() );
  outFileName.append(".vtk");

  vtkSmartPointer<vtkPolyDataWriter> pDataWriter =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  pDataWriter->SetFileName( outFileName.c_str() );
  pDataWriter->SetFileTypeToBinary();
  pDataWriter->SetInput(polydata);
  pDataWriter->Write();
}


// dump the Delaunay-triangulated grid of an atom
void aimsAtomCenteredGrid::dumpUGrid(int atomId)
{
  string outFileName;
  stringstream ss;
  ss << atomId;
  outFileName.assign("atomUGrid");
  outFileName.append( ss.str() );
  outFileName.append(".vtk");

  vtkSmartPointer<vtkUnstructuredGridWriter> ugridWriter =
      vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  ugridWriter->SetFileName( outFileName.c_str() );
  ugridWriter->SetFileTypeToBinary();
  ugridWriter->SetInput(ugrid);
  ugridWriter->Write();
}


// compute and return the interpolated value for "point"
//double aimsAtomCenteredGrid::query(double *point, vtkSmartPointer<vtkCell> *previousCell, vtkIdType *previousCellId)
vector<double> aimsAtomCenteredGrid::query(double *point, vtkSmartPointer<vtkCell> *previousCell, vtkIdType *previousCellId)
{
  vtkSmartPointer<vtkCell> currentCell=NULL;
  vtkIdType currentCellId=-1;
  vtkSmartPointer<vtkIdList> cellPointIds; // list for the points in the current cell
  double weights[16];  // interpolation weights
  double tol2=1.e-6;   // tolerance, some small number
  int    subId=0;      // not so sure ...
  double pcoords[3];   // cell-related coordinates, unused
  double value;

  vector<double> data;
  data.clear();

  // find the cell in the unstructured grid where the
  // current point in the rectilinear grid is located
  currentCellId=ugrid->FindCell( point, *previousCell, *previousCellId,
      tol2, subId, pcoords, weights );

  if(currentCellId>=0)  // did we find the cell?
  {
    currentCell=ugrid->GetCell(currentCellId);
    //
    // Tetrahedral grid -> is the cell defined by 4 points?
    // This should always be the case!
    if(currentCell->GetNumberOfPoints()==4)
    {
      cellPointIds=currentCell->GetPointIds();
      //
      for (int j=0; j < (int)vtkDataArrays.size(); j++)
      {
        value=0.0;
        //
        for (vtkIdType i=0; i<4; i++)
        {
          // do the interpolation based on the weights we got from FindCell
          value += (weights[i] * vtkDataArrays[j]->GetValue( cellPointIds->GetId(i) ));
        }
        //
        data.push_back(value);
      }
      //
      // saving the information about where to start reduces the runtime
      *previousCell  =currentCell;
      *previousCellId=currentCellId;
    }

  }
  else
  {
    *previousCellId=-1;
  }

  return data;
}


void aimsAtomCenteredGrid::setArrayNames(vector<string> names)
{
  arrayNames = names;
}
