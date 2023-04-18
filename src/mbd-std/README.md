# `mbdvdw` — Many-body dispersion method

This Fortran package implements the Many-body dispersion method including analytical forces and stress, self-consistency and reciprocal space/Ewald summation for periodic systems.

The code is used in two electronic-structure codes: [Quantum Espresso](http://www.quantum-espresso.org) and [FHI-aims](https://aimsclub.fhi-berlin.mpg.de).

The included Python module wraps the non-self-consistent part of the code and enables to use the package as a standalone program.

### Installation

The installation requires Python 2 or 3 with Numpy, [cffi](https://cffi.readthedocs.io/en/latest/) and [mpi4py](https://mpi4py.readthedocs.io/). All these are included in the [Anaconda](https://anaconda.org) Python.

To build the Python extension, simply run `make`. The build may require adjustments via environment variables. The default run is equivalent to

```bash
make FC=gfortran MPIFC=mpifort PYTHON=python3 LAPACK='-framework Accelerate' FFLAGS=''
```

and should work on macOS with the scientific build stack installed.

To test whether the build was successful, run

```bash
python3 pymbdvdw.py
```

### Usage

See `pymbdvdw.py`.
