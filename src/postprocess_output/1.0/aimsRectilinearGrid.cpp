#include <iostream>
#include <string>
#include <sstream>

#include <aimsRectilinearGrid.h>
#include <simpleTimer.h>
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
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>



using namespace std;

void aimsRectilinearGrid::initialize(
    int nX, int nY, int nZ,
    double xMin, double xMax,
    double yMin, double yMax,
    double zMin, double zMax)
{
//  simpleTimer T("aimsRectilinearGrid::initialize()");

  nx=nX;
  ny=nY;
  nz=nZ;
  xmin=xMin;
  xmax=xMax;
  ymin=yMin;
  ymax=yMax;
  zmin=zMin;
  zmax=zMax;

//  cout << "  Creating rectilinear grid with ..." << endl;
//  cout << "    (nx, ny, nz)=(" << nx   << ", " << ny   << ", " << nz << ")" << endl;
//  cout << "    (xmin, xmax)=(" << xmin << ", " << xmax               << ")" << endl;
//  cout << "    (ymin, ymax)=(" << ymin << ", " << ymax               << ")" << endl;
//  cout << "    (zmin, zmax)=(" << zmin << ", " << zmax               << ")" << endl;

  cartX = new double[nx];
  cartY = new double[ny];
  cartZ = new double[nz];
  {
    cartX[0]=xmin;
    double dx=(xmax-xmin)/double(nx-1);
    for (int ix=1; ix<nx; ix++)
      cartX[ix]=cartX[ix-1]+dx;
    cartY[0]=ymin;
    double dy=(ymax-ymin)/double(ny-1);
    for (int iy=1; iy<ny; iy++)
      cartY[iy]=cartY[iy-1]+dy;
    cartZ[0]=zmin;
    double dz=(zmax-zmin)/double(nz-1);
    for (int iz=1; iz<nz; iz++)
      cartZ[iz]=cartZ[iz-1]+dz;
  }

  rectGrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  rectGrid->SetDimensions(nx,ny,nz);

  vtkSmartPointer<vtkDoubleArray> xArray =
      vtkSmartPointer<vtkDoubleArray>::New();
  for (int ix=0; ix<nx; ix++)
    xArray->InsertNextValue( cartX[ix] );
  vtkSmartPointer<vtkDoubleArray> yArray =
      vtkSmartPointer<vtkDoubleArray>::New();
  for (int iy=0; iy<ny; iy++)
    yArray->InsertNextValue( cartY[iy] );
  vtkSmartPointer<vtkDoubleArray> zArray =
      vtkSmartPointer<vtkDoubleArray>::New();
  for (int iz=0; iz<nz; iz++)
    zArray->InsertNextValue( cartZ[iz] );

  rectGrid->SetXCoordinates(xArray);
  rectGrid->SetYCoordinates(yArray);
  rectGrid->SetZCoordinates(zArray);

  verbose = false;
  debug = false;
} // end initialize()


void aimsRectilinearGrid::initializeDataArrays(vector<string> dataNames)
{
  int n = dataNames.size();
  //
  for (int j=0; j<n; j++)
  {
    vtkDataArrays.push_back( vtkSmartPointer<vtkDoubleArray>::New() );

    char *s;
    s = new char[32];
    sprintf(s, "array-%d", j);
    vtkDataArrays.back()->SetName( s );
    delete [] s;

//    vtkDataArrays.back()->SetName( dataNames[j].c_str() );

    {
      // blank the data array
      vtkIdType linDim=nx*ny*nz;
      vtkDataArrays.back()->SetNumberOfTuples(linDim);
      for (vtkIdType i=0; i<linDim; i++)
        vtkDataArrays.back()->SetValue(i, 0.0);
    }
  }
}


void aimsRectilinearGrid::destroy()
{
  delete [] cartX;
  delete [] cartY;
  delete [] cartZ;
//  rectData->Reset();
//  rectData->Squeeze();
  for (int i=0; i < (int)vtkDataArrays.size(); i++)
  {
    vtkDataArrays[i]->Reset();
    vtkDataArrays[i]->Squeeze();
  }
  vtkDataArrays.clear();
};

// constructor with some defaults
aimsRectilinearGrid::aimsRectilinearGrid()
{
  initialize(128, 128, 128, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);
};

aimsRectilinearGrid::aimsRectilinearGrid(
    int &nX, int &nY, int &nZ,
    double &xMin, double &xMax,
    double &yMin, double &yMax,
    double &zMin, double &zMax)
{
  initialize(nX, nY, nZ, xMin, xMax,
      yMin, yMax, zMin, zMax);
};

aimsRectilinearGrid::~aimsRectilinearGrid()
{
  destroy();
};

void aimsRectilinearGrid::addFromAimsAtomCenteredGrid( aimsAtomCenteredGrid * atom )
{
//  simpleTimer T("aimsRectilinearGrid::addFromAimsAtomCenteredGrid()");

  vtkSmartPointer<vtkCell> previousCell=NULL;
  vtkIdType previousCellId=-1;
  double point[3];         // coordinates of the current point
  double value=0.0;        // dummy variable to hold and update the rectGrid value
  vtkIdType rectGridIdx=0; // count linearly through the rectGrid data array
  vector<double> data;

  int nhit = 0;
  int nmiss = 0;

  for (int iz=0; iz<nz; iz++)
  {
    for (int iy=0; iy<ny; iy++)
    {
      for (int ix=0; ix<nx; ix++)
      {
        point[0]=cartX[ix];
        point[1]=cartY[iy];
        point[2]=cartZ[iz];

        data.clear();

        data = atom->query(point, &previousCell, &previousCellId);

        if(previousCellId >= 0) // query could be satisfied
        {
          //assert(data.size() == vtkDataArrays.size());
          //
          for (int i=0; i < (int)vtkDataArrays.size(); i++)
          {
            value = vtkDataArrays[i]->GetValue(rectGridIdx);
            value += data[i];
            vtkDataArrays[i]->SetValue(rectGridIdx, value);
          }
          ++nhit;
        }
        else
        {
          ++nmiss;
        }

        rectGridIdx++;
      }
    }
  } // end for iz ...

  if (debug)
    cout << "  nhit, nmiss = " << nhit << ", " << nmiss << endl;
} // end function


// access data
double aimsRectilinearGrid::getDataValue(ssize_t &dataIdx, vtkIdType &index)
{
  return ( vtkDataArrays[dataIdx]->GetValue(index) );
}


// prepare the data array
void aimsRectilinearGrid::setDataArray(double *dataArray, vtkIdType dataSize, string dataName, ssize_t dataIdx)
{
  assert( dataSize==(nx*ny*nz) );
  vtkDataArrays[dataIdx]->Reset();
  vtkDataArrays[dataIdx]->SetNumberOfTuples(dataSize);

  //vtkDataArrays[dataIdx]->SetArray(dataArray, dataSize, 1);
  //... does *not* work (segfault)!  Use instead:
  for (vtkIdType i=0; i<dataSize; i++)
    vtkDataArrays[dataIdx]->SetValue(i, dataArray[i]);

  vtkDataArrays[dataIdx]->SetName( dataName.c_str() );
}


// write data as vtr file
void aimsRectilinearGrid::dumpVtr(string &fileName)
{
  for (int i=0; i < (int)vtkDataArrays.size(); i++)
  {
    rectGrid->GetPointData()->AddArray( vtkDataArrays[i] );
//    cout << "###" << vtkDataArrays[i]->GetName() << endl;
  }
  vtkSmartPointer<vtkXMLRectilinearGridWriter> writer =
      vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(rectGrid);
  //writer->SetDataModeToAscii();
  writer->SetDataModeToBinary();
  writer->Write();
};



// write data as vti file
void aimsRectilinearGrid::dumpVti(string &fileName)
{
  double dx, dy, dz;
  {
    double eps;
    eps = 1.e-6;

    if (nx>1)
      dx = (xmax - xmin)/double(nx-1);
    else
      dx = -1.;

    if (ny>1)
      dy = (ymax - ymin)/double(ny-1);
    else
      dy = -1.;

    if (nz>1)
      dz = (zmax - zmin)/double(nz-1);
    else
      dz = -1.;

    // TODO   check if deltas are equal
  }

  vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
  // Specify the size of the image data
  imageData->SetDimensions( nx, ny, nz );
  imageData->SetOrigin( xmin, ymin, zmin );
  imageData->SetSpacing( dx, dy, dz );

  // TODO  implement support for more than one scalar field
  imageData->SetNumberOfScalarComponents(1);
  imageData->SetScalarTypeToDouble();
  int* dims = imageData->GetDimensions();

  vtkIdType rectGridIdx=0;
  for (int z = 0; z < dims[2]; z++)
  {
    for (int y = 0; y < dims[1]; y++)
    {
      for (int x = 0; x < dims[0]; x++)
      {
        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
        // TODO  implement support for more than one scalar field
        pixel[0] = vtkDataArrays[0]->GetValue(rectGridIdx);
        rectGridIdx++;
      }
    }
  }

  vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(imageData);
  writer->Write();
};



void aimsRectilinearGrid::setVerbosity(bool _verbose, bool _debug)
{
  verbose = _verbose;
  debug = _debug;
};
