# ---
#  Makefile for the FHI-aims data postprocessing tool 'a2rg'.
#
#  Requires VTK (5.8.0) and MPI.  On RZG Linux clusters,
#  type "module load impi/4.0.3 vtk/5.8" before compiling and
#  running this code.
# ---


### --- executable, source file, compiler settings ---
PROG := a2rg
#SUFX := .exe
SRCS := a2rg.cpp aimsAtomCenteredGrid.cpp aimsAtomCenteredGridBundle.cpp aimsRectilinearGrid.cpp simpleTimer.cpp chemicalSymbols.cpp

# Intel MPI compiler wrapper for g++ (available on all clusters, e.g. AIMS, after 'module load impi/4.0.3'):
CCC := mpigxx

# *** uncomment the following two lines for the local OpenMPI and VTK on the VIZ cluster ***
#CCC := mpic++
#VTK_HOME := /u/system/vtk/5.8.0

CCFLAG := -O3 -msse2 -Wall -Wno-deprecated
CCFLAG += -DMPICH_IGNORE_CXX_SEEK
OBJS := ${SRCS:.cpp=.o}
HDRS := ${SRCS:.cpp=.h}
INCL := -I.

### --- use libqhullvtk to run the delaunay triangulation (experimental) ---
QHVTK := false

### --- VTK ---
INCL  += -I$(VTK_HOME)/include/vtk
LIBD  += -L$(VTK_HOME)/lib/vtk
RPATH += -Wl,-rpath,$(VTK_HOME)/lib/vtk
#
LIBS += -lCosmo -lvtkexpat -lvtkHybrid -lvtkNetCDF -lvtkverdict
LIBS += -lMapReduceMPI -lvtkFiltering -lvtkImaging -lvtkParallel -lvtkViews
LIBS += -lVPIC -lvtkfreetype -lvtkInfovis -lvtkpng -lvtkVolumeRendering
LIBS += -lvtkalglib -lvtkftgl -lvtkIO -lvtkproj4 -lvtkWidgets
LIBS += -lvtkCharts -lvtkGenericFiltering -lvtkjpeg -lvtkRendering -lvtkzlib
LIBS += -lvtkCommon -lvtkGeovis -lvtklibxml2 -lvtksqlite
LIBS += -lvtkDICOMParser -lvtkGraphics -lvtkmetaio -lvtksys
LIBS += -lvtkexoIIc -lvtkhdf5 -lvtkNetCDF_cxx -lvtktiff
#
LIBS += -lm

### --- libqhullvtk, currently not used ---
ifeq (${QHVTK},true)
    CCFLAG += -DHAVE_QHULLVTK
    QHVTK_HOME := /afs/ipp/home/k/khr/RZG/Delaunay/delaunay-svn/trunk/qhullvtk
    INCL += -I$(QHVTK_HOME)
    LIBD += -L$(QHVTK_HOME)
    LIBS += -lqhullvtk
    RPATH += -Wl,-rpath,$(QHVTK_HOME)
endif


### --- COMPILATION ---

all: default

default: $(PROG)$(SUFX)

$(OBJS): %.o: %.cpp
	$(CCC) $(CCFLAG) $(INCL) -o $@ -c $<

$(PROG)$(SUFX): Makefile $(OBJS)
	$(CCC) $(CCFLAG) $(LIBD) $(LIBS) $(RPATH) -o $(PROG)$(SUFX) $(OBJS)

.PHONY: clean realclean

realclean: clean
	rm -f $(PROG)$(SUFX)

clean:
	rm -f *.o
