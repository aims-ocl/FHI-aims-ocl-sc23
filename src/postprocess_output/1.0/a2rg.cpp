/**
 *
 *  a2rg
 *
 *  purpose:
 *      map FHI-aims atom-centered grid data to a rectilinear grid
 *
 *  input:
 *      data on atom-centered grids from FHI-aims in binary or ascii format
 *
 *  output:
 *      data on a VTK rectilinear grid (.vtr file)
 *      data as 3D VTK image data file (.vti file)
 *      atom coordinates in VTK polydata format (.vtp file)
 *
 */



// C include files
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/types.h>

// MPI include file
#include <mpi.h>

// C++ include files
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>

// C++/VTK include files
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMath.h>

// include files from the package
#include <aimsAtomCenteredGrid.h>
#include <aimsAtomCenteredGridBundle.h>
#include <aimsRectilinearGrid.h>
#include <chemicalSymbols.h>
#include <simpleTimer.h>

using namespace std;

// coordinate correction factor
#define BOHRFACTOR 0.52917721


void usage(char *exeName)
{
  cerr << "Usage: " << exeName << " [OPTION]... FHI-aims-data-file[s]" << endl;
  cerr << endl;
  cerr << "  Options" << endl;
  cerr << "    --help     print this help message" << endl;
  cerr << "    --ascii    read data files in ascii mode, default: binary" << endl;
  cerr << "    --verbose  print information on what is being done"  << endl;
  cerr << "    --debug    print debug information"  << endl;
  cerr << "    --timing   print timing output" << endl;
  cerr << endl;
  cerr << "  Options with arguments" << endl;
  cerr << "    --nx       <number of grid points in x direction, default: 128>" << endl;
  cerr << "    --ny       <number of grid points in y direction, default: 128>" << endl;
  cerr << "    --nz       <number of grid points in z direction, default: 128>" << endl;
  cerr << "    --xmin     <coordinate minimum in x direction, default: auto>" << endl;
  cerr << "    --xmax     <coordinate maximum in x direction, default: auto>" << endl;
  cerr << "    --ymin     <coordinate minimum in y direction, default: auto>" << endl;
  cerr << "    --ymax     <coordinate maximum in y direction, default: auto>" << endl;
  cerr << "    --zmin     <coordinate minimum in z direction, default: auto>" << endl;
  cerr << "    --zmax     <coordinate maximum in z direction, default: auto>" << endl;
  cerr << "    --vtrout   <VTK rectilinear grid output file name (.vtr), default: disabled>" << endl;
  cerr << "    --vtiout   <VTK 3D image data output file name (.vti), default: disabled." << endl;
  cerr << "                (.vti files are required only to perform volume rendering.)>" << endl;
  cerr << "    --vtpout   <VTK atom coordinates output file name (.vtp), default: disabled>" << endl;
  cerr << endl;
}



int main(int argc, char **argv)
{
  // MPI initialization
  int mpiSize;
  int mpiRank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);


  bool timing  = false;
  bool verbose = false;
  bool debug   = false;
  bool binary  =  true;

  int reclen  =     32;
  int nAtoms  =     -1;
  int nx      =    128;
  int ny      =    128;
  int nz      =    128;
  double xmin =  1.e99;
  double ymin =  1.e99;
  double zmin =  1.e99;
  double xmax = -1.e99;
  double ymax = -1.e99;
  double zmax = -1.e99;

  // a "small number" for random coordinate modifications
  double  eps =  2.5e-5;

  string vtrFileName;
  vtrFileName.assign("disable");
  string vtiFileName;
  vtiFileName.assign("disable");
  string atomsFileName;
  atomsFileName.assign("disable");

  // names of the data sets
  vector<string> fieldIds;

  // parse the argument vector, handle command line parameters
  while (true)
  {
    static struct option long_options[] =
        {
            {"help",   no_argument,       0, 'h'},
            {"ascii",  no_argument,       0, 'b'},
            {"debug",  no_argument,       0, 'd'},
            {"verbose",no_argument,       0, 'v'},
            // number of points in each direction
            {"nx",     required_argument, 0, 'l'},
            {"ny",     required_argument, 0, 'm'},
            {"nz",     required_argument, 0, 'n'},
            // coordinate limits for the grid
            {"xmin",   required_argument, 0, 'p'},
            {"xmax",   required_argument, 0, 'q'},
            {"ymin",   required_argument, 0, 'r'},
            {"ymax",   required_argument, 0, 's'},
            {"zmin",   required_argument, 0, 't'},
            {"zmax",   required_argument, 0, 'u'},
            // input/output file names
            {"vtrout", required_argument, 0, 'o'},
            {"vtpout", required_argument, 0, 'a'},
            {"vtiout", required_argument, 0, 'i'},
            // end of structure
            {0, 0, 0, 0}
        };

    int option_index = 0;
    int c;
    c = getopt_long(argc, argv, "hbdvl:m:n:p:q:r:s:t:u:o:a:i:",
        long_options, &option_index);

    if (c == -1) break;

    switch (c)
    {
      case 'h':
        usage(argv[0]);
        return 1;
        break;
      case 'b':
        binary = false;
        break;
      case 'l':
        nx = atoi(optarg);
        break;
      case 'm':
        ny = atoi(optarg);
        break;
      case 'n':
        nz = atoi(optarg);
        break;
      case 'p':
        xmin = atof(optarg);
        break;
      case 'q':
        xmax = atof(optarg);
        break;
      case 'r':
        ymin = atof(optarg);
        break;
      case 's':
        ymax = atof(optarg);
        break;
      case 't':
        zmin = atof(optarg);
        break;
      case 'u':
        zmax = atof(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      case 'd':
        debug = true;
        break;
      case 'o':
        vtrFileName.assign(optarg);
        break;
      case 'i':
        vtiFileName.assign(optarg);
        break;
      case 'a':
        atomsFileName.assign(optarg);
        break;
      default:
        break;
    }
  }

  if (    (vtrFileName.compare("disable")   == 0)
       || (vtiFileName.compare("disable")   == 0)
       || (atomsFileName.compare("disable") == 0) )
  {
    cout << "No output file name given -- nothing to do.  Exiting." << endl;
    return 1;
  }

  if (timing)
    simpleTimer T("main()::TOTAL", mpiRank);

  // read record length from the first data file
  if (binary)
  {
    char *inFileName = argv[optind];

    if (verbose)
      cout << "Reading record length from " << inFileName << " ..." << endl;

    FILE *inFileP = fopen (inFileName, "rb");
    assert(inFileP != NULL);
    char line[256];
    while( fread(line, reclen, 1, inFileP) )
    {
      if (! strncmp("RECL", line, 4))
      {
        memcpy(&reclen, &(line[5]), sizeof(reclen));
        break;
      }
    }
    fclose(inFileP);

    if (debug)
      cout << "  reclen = " << reclen << endl;
  }


  // read nAtoms from the first data file
  {
    char *inFileName = argv[optind];

    if (verbose)
      cout << "Reading nAtoms from " << inFileName << " ..." << endl;

    FILE *inFileP = fopen (inFileName, (binary ? "rb" : "r"));
    assert(inFileP != NULL);
    char line[256];
    while(    (binary && (fread(line, reclen, 1, inFileP)) )
           || (fgets( line, sizeof(line), inFileP ) != NULL) )
    {
      if (! strncmp("ATOMS", line, 5))
      {
        if (binary)
        {
          memcpy(&nAtoms, &(line[6]), sizeof(nAtoms));
        }
        else
        {
          char tmpBuf[64];
          sscanf(line, "%s %d", tmpBuf, &nAtoms);
        }
        break;
      }
    }
    fclose(inFileP);

    if (debug)
      cout << "  nAtoms = " << nAtoms << endl;
  }


  // calculate MPI offsets
  int mpiNAtoms=nAtoms/mpiSize;
  int mpiAtom0=mpiRank*mpiNAtoms;
  int mpiAtom1;
  if (mpiRank<(mpiSize-1))
  {
    mpiAtom1=(mpiRank+1)*mpiNAtoms-1;
  }
  else
  {
    mpiAtom1=nAtoms-1;
  }
  mpiNAtoms=mpiAtom1-mpiAtom0;

  aimsAtomCenteredGridBundle myAtoms(mpiAtom0, mpiAtom1, nAtoms);

  // contract the coordinate extents if there's only a single grid point
  if (nx==1) xmax = xmin;
  if (ny==1) ymax = ymin;
  if (nz==1) zmax = zmin;

  double coordMin[3];    // autodetection: minima of the coordinates
  double coordMax[3];    //                maxima of the coordinates
  for (int j=0; j<3; j++)
  {
    coordMin[j] = +1.e99;  // useful initialization for
    coordMax[j] = -1.e99;  // finding the true min and max
  }

  int currentAtom;
  // (1) read the data files
  {
    if (timing)
      simpleTimer T("main()::DATA READ-IN", mpiRank);

    if (verbose)
      cout << "Reading data ..." << endl;

    // Read all data files given as command line arguments:
    //   <optind> is set by getopt_long() to the first entry in argv[]
    //   after all regular options which is in our case the list of
    //   files to process.
    for (int i=optind; i<argc; i++)
    {
      char *inFileName = argv[i];

      if (verbose)
        cout << inFileName << endl;

      FILE *inFileP = fopen (inFileName, (binary ? "rb" : "r"));
      assert(inFileP != NULL);
      currentAtom=-1;

      char line[256];
      char tmpBuf[64];
      double xyz[3];       // buffer for the coordinate triple

      // TEMPORARY FIX:  Bug in AIMS-output: The FIELDS entry is present in the dataset twice.
//      bool DebugFirstFields = true;

      while(    (binary && (fread(line, reclen, 1, inFileP)) )
             || (fgets( line, sizeof(line), inFileP ) != NULL) )
      {
        if (! strncmp("ATOMS", line, 5))
        {
          // do nothing
        }
        else if (! strncmp("RECL", line, 4))
        {
          // do nothing
        }
        else if (! strncmp("FIELDS", line, 6))
        {
          string _line;
          _line = &(line[7]);
          istringstream iss( _line );
          copy(istream_iterator<string>(iss),
              istream_iterator<string>(),
              back_inserter<vector<string> >(fieldIds));

          {
            vector<string> fieldIdsTmp;
            fieldIdsTmp = fieldIds;

            sort(fieldIdsTmp.begin(), fieldIdsTmp.end());

            bool duplicate;
            duplicate = false;

            for(int i = 1; i < (int)fieldIdsTmp.size(); i++)
            {
              if (fieldIdsTmp[i] == fieldIdsTmp[i-1])
              {
                duplicate = true;
                break;
              }
            }

            if (duplicate)
            {
              if (debug)
                cout << "Duplicate fieldId encountered, using default fieldIds!" << endl;
              char *s;
              s = new char[16];
              for(int i = 0; i < (int)fieldIds.size(); i++)
              {
                sprintf(s, "Array-%d", i);
                fieldIds[i] = s;
              }
              delete [] s;
            }
          }


          if (debug)
          {
            cout << "  fieldIds = ";
            vector<string>::iterator it;
            for(it = fieldIds.begin(); it < fieldIds.end(); it++)
            {
              cout << *it << " ";
            }
            cout << endl;
          }
        }
        else if (! strncmp("CENTER", line, 6))
        {
          // required for a sticks-and-balls model (output to .vtp file)
          if (binary)
            memcpy(&xyz[0], &(line[7]), 3*sizeof(double));
          else
            sscanf(line, "%s %lg %lg %lg", tmpBuf, &xyz[0], &xyz[1],  &xyz[2]);
          //
#ifdef BOHRFACTOR
          for (int j=0; j<3; j++)
            xyz[j] *= BOHRFACTOR;
#endif
          //
          for (int j=0; j<3; j++)
          {
            // move the coordinate randomly by an epsilon,
            // otherwise VTK data structures are completely flat in certain cases
            double mov=vtkMath::Random(-eps, eps);
            xyz[j] = xyz[j] + mov;
          }
          //
          myAtoms.setCenter(xyz);

          if (debug)
            cout << "  centerXYZ = " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;

        }
        else if (! strncmp("SPECIES", line, 7))
        {
          // required for a sticks-and-balls model (output to .vtp file)
          char elemBuf[3];
          if (binary)
          {
            // TODO : TEST!
            memcpy(&elemBuf, &(line[8]), 2*sizeof(char));
            if ( (elemBuf[1]>=65) && (elemBuf[1]<=90) )
              elemBuf[2] = '\0';
            else
              elemBuf[1] = '\0';
          }
          else
          {
            sscanf(line, "%s %s", tmpBuf, elemBuf);
          }
          myAtoms.setAtomicNumber( symbol_to_z(string(elemBuf)) );
          //
          if (debug)
            cout << "  species = " << elemBuf << endl;
        }
        else if (! strncmp("ATOM", line, 4))
        {
          // switch to a given ATOM ID
          if (binary)
            memcpy(&currentAtom, &(line[5]), sizeof(currentAtom));
          else
            sscanf(line, "%s %d", tmpBuf, &currentAtom);
          currentAtom=currentAtom-1;   // mind Fortran vs C counting
          myAtoms.setCurrentAtomIndex(currentAtom);
          //
          if (debug)
            cout << "  currentAtom = " << currentAtom << endl;
        }
        else
        {
          // add data only iff the current MPI process is responsible
          if ((currentAtom >= mpiAtom0) && (currentAtom <= mpiAtom1))
          {
            int ir; // radial index
            int ia; // angular index
            double *data = new double[fieldIds.size()];
            //
            if (binary)
            {
              int offset;
              offset = 0;
              memcpy(&xyz[0],    &(line[offset]), 3*sizeof(double));
              offset += 3*sizeof(double);
              memcpy(&ir,        &(line[offset]), sizeof(ir));
              offset += sizeof(ir);
              memcpy(&ia,        &(line[offset]), sizeof(ia));
              offset += sizeof(ia);
              // TODO : add variable number of scalars
              memcpy(&(data[0]), &(line[offset]), fieldIds.size()*sizeof(double));
            }
            else
            {
              string _line;
              _line = &(line[7]);
              istringstream iss( _line );
              vector<string> fieldValues;
              copy(istream_iterator<string>(iss),
                   istream_iterator<string>(),
                   back_inserter<vector<string> >(fieldValues));
              //
              assert( fieldValues.size() == (fieldIds.size() + 5) );
              //
              xyz[0] = atof( fieldValues[0].c_str() );
              xyz[1] = atof( fieldValues[1].c_str() );
              xyz[2] = atof( fieldValues[2].c_str() );
              ir     = atoi( fieldValues[3].c_str() );
              ia     = atoi( fieldValues[4].c_str() );
              //
              for (size_t j=0; j<fieldIds.size(); j++)
                data[j] = atof( fieldValues[j+5].c_str() );
            }
            //
#ifdef BOHRFACTOR
            for (int j=0; j<3; j++)
              xyz[j] *= BOHRFACTOR;
#endif
            //
            for (int j=0; j<3; j++)
            {
              // move the coordinate randomly by an epsilon,
              // otherwise the Delaunay triangulation has problems
              double mov=vtkMath::Random(-eps, eps);
              xyz[j] = xyz[j] + mov;
            }
            //
            for (int j=0; j<3; j++)
            {
              if (coordMin[j]>xyz[j])  coordMin[j]=xyz[j];
            }
            for (int j=0; j<3; j++)
            {
              if (coordMax[j]<xyz[j])  coordMax[j]=xyz[j];
            }
            //
            vector <double> tmpPointData;
            tmpPointData.clear();
            tmpPointData.push_back(xyz[0]);
            tmpPointData.push_back(xyz[1]);
            tmpPointData.push_back(xyz[2]);
            tmpPointData.push_back(double(ir));
            tmpPointData.push_back(double(ia));
            for (size_t j=0; j<fieldIds.size(); j++)
              tmpPointData.push_back(data[j]);
            //
            myAtoms.addPoint(tmpPointData);
            //
            if (false && debug)
            {
              cout << "  pointData = ";
              vector<double>::iterator it;
              for(it = tmpPointData.begin(); it < tmpPointData.end(); it++)
              {
                cout << *it << " ";
              }
              cout << endl;
            }
            //
            delete[] data;
          }
        }
      }
      fclose (inFileP);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  } // data read-in loop



  // sync the minima and maxima between all processes
  MPI_Allreduce(MPI_IN_PLACE, coordMin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, coordMax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  // Impose the min and max values for the box:
  // the "unrealistically initialized" values may have been altered via the command line
  if (xmin < 1.e99) coordMin[0]=xmin;
  if (ymin < 1.e99) coordMin[1]=ymin;
  if (zmin < 1.e99) coordMin[2]=zmin;
  if (xmax >-1.e99) coordMax[0]=xmax;
  if (ymax >-1.e99) coordMax[1]=ymax;
  if (zmax >-1.e99) coordMax[2]=zmax;



  // calculate unstructured grids from polydata
  {
    if (timing)
      simpleTimer T("main()::DELAUNAY3D", mpiRank);
    myAtoms.setArrayNames(fieldIds);
    myAtoms.calculateGrids(verbose);
    MPI_Barrier(MPI_COMM_WORLD);
  }


  // dump polydata and the triangulated unstructured grids for debugging
  if (false)
  {
    myAtoms.atoms[1].dumpPolyData(1);
    myAtoms.atoms[2].dumpPolyData(2);
    myAtoms.atoms[1].dumpUGrid(1);
    myAtoms.atoms[2].dumpUGrid(2);
  }


  // (3) map to rectilinear grids in parallel
  {
    aimsRectilinearGrid *R;
    R = new aimsRectilinearGrid( nx, ny, nz,
        coordMin[0], coordMax[0],
        coordMin[1], coordMax[1],
        coordMin[2], coordMax[2]);
    R->setVerbosity(verbose, debug);
    R->initializeDataArrays( fieldIds );


    // map each atom grid to the rectilinear grid
//    if (false)
    {
      if (timing)
        simpleTimer T("main()::MAPPING OF DATA TO RECT GRID", mpiRank);
      //
      for (currentAtom=mpiAtom0; currentAtom<=mpiAtom1; currentAtom++)
      {
        if (verbose)
          cout << "Mapping atom-centered grid #" << currentAtom << " ..." << endl;
        myAtoms.setCurrentAtomIndex(currentAtom);
        R->addFromAimsAtomCenteredGrid( myAtoms.getCurrentAtom() );
        //
        if (mpiRank > 0) // mpiRank 0 still needs the data, see below
          myAtoms.freeAtom();
      }
      //
      MPI_Barrier(MPI_COMM_WORLD);
    }


    // global summation
    {
      if (verbose)
        cout << "Performing global data summation ..." << endl;
      if (timing)
        simpleTimer T("main()::MPI DATA SUMMATION", mpiRank);
      //
      for (ssize_t i=0; i < (ssize_t)fieldIds.size(); i++)
      {
        vtkIdType indexXYZ;
        indexXYZ=0;
        double *bufferXY;
        double *bufferXYZ;
        bufferXY = new double[nx*ny];
        bufferXYZ = new double[nx*ny*nz];
        for (int iz=0; iz<nz; iz++)
        {
          vtkIdType indexXY;
          indexXY=0;
          for (int iy=0; iy<ny; iy++)
          {
            for (int ix=0; ix<nx; ix++)
            {
              bufferXY[indexXY] = R->getDataValue(i, indexXYZ);
              indexXYZ++;
              indexXY++;
            }
          }
          MPI_Reduce(bufferXY, &(bufferXYZ[iz*(nx*ny)]), (nx*ny), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        } // end for iz ...
        //
        MPI_Barrier(MPI_COMM_WORLD);
        //
        if (mpiRank == 0)
        {
          R->setDataArray(bufferXYZ, (nx*ny*nz), fieldIds[i], i);
        }
        //
        delete [] bufferXY;
        delete [] bufferXYZ;
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }

    if (mpiRank == 0)
    {
      if (verbose)
        cout << "Writing data ..." << endl;

      if (timing)
        simpleTimer T("main()::FINAL IO", mpiRank);

      // output data for a sticks-and-balls model as .vtp file
      if ( atomsFileName.compare("disable")!=0 )
      {
        vtkSmartPointer<vtkPoints>    points       = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> vertices     = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkIntArray>  atomicNumber = vtkSmartPointer<vtkIntArray>::New();
        vtkIdType pid[1];
        //
        atomicNumber->SetNumberOfComponents(1);
        atomicNumber->SetName("element");
        //
        for (int i=0; i<nAtoms; i++)
        {
          myAtoms.setCurrentAtomIndex(i);
          //
          pid[0] = points->InsertNextPoint( myAtoms.getCenter() );
          vertices->InsertNextCell( 1, pid );
          atomicNumber->InsertNextValue( myAtoms.getAtomicNumber() );
          //
          myAtoms.freeAtom();
        }
        //
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkXMLPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        //
        polydata->SetPoints(points);
        polydata->SetVerts(vertices);
        polydata->GetPointData()->AddArray(atomicNumber);
        pwriter->SetFileName( atomsFileName.c_str() );
        pwriter->SetInput(polydata);
        pwriter->SetDataModeToBinary();
        pwriter->Write();
      }

      // write data to .vtr file
      if ( vtrFileName.compare("disable")!=0 )
        R->dumpVtr( vtrFileName );

      // write data to .vti file (for volume rendering)
      if ( vtiFileName.compare("disable")!=0 )
        R->dumpVti( vtiFileName );

/*
      // debug output for GNUPLOT
      if (nx==1)
      {
        indexXYZ=0;
        FILE *gpfile;
        char fileName[32];
        sprintf(fileName, "gnuplot.dat");
        gpfile=fopen(fileName, "w");
        //
        for (int iz=0; iz<nz; iz++)
        {
          for (int iy=0; iy<ny; iy++)
          {
            {
              fprintf(gpfile, "%e\n", bufferXYZ[indexXYZ]);
              indexXYZ++;
            }
          }
          fprintf(gpfile, "\n");
        }
        fclose(gpfile);
      }
*/

    } // end of data output on rank 0

    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
