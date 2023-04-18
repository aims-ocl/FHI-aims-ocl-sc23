"""
Implements Python 2/3 interface to the FHI-aims C bindings. See
python_interface.f90 for the Fortran side of the interface.

COPYRIGHT

Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note
that any use of the "FHI-aims-Software" is subject to the terms and
conditions of the respective license agreement.
"""
from __future__ import print_function
import sys
import os
from pathlib import Path
import imp
import numpy as np
import traceback
from collections import namedtuple
import code
try:
    from c_aims import ffi, lib as aims
except ImportError:  # enable importing the module outside of cffi context
    class FFI:
        def def_extern(self):
            def decorator(f):
                return f
            return decorator
    ffi = FFI()
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
except ImportError:  # support Pythons withou mpi4py
    comm = None
    rank = 0


# Python counterparts of the C structs in python_interface.h
GridPoint = namedtuple('GridPoint', 'coords index_atom index_radial index_angular')
BatchOfPoints = namedtuple('BatchOfPoints', 'size points')


# Python counterpart to AimsContext_t
class AimsContext:
    """
    Contains data from FHI-aims. The data components are accessed as attributes
    and cannot be assigned to, only modified in-place. All data are pure Python
    objects. The Numpy arrays are linked to the same memory locations as the
    corresponding Fortran arrays and hence are in Fortran byte order.
    Furthermore, the object defines convenience methods that wrap some MPI
    synchronisation.

    Examples:
    >>> print(ctx.coords)
    >>> ctx.rho[:] = 0.d0
    >>> print(len(ctx.batches[0].points))
    >>> points, rho = ctx.gather_all_grids(['rho'])
    """
    def __init__(self, dict_):
        self.__dict__['_dict'] = dict_
        self.__dict__['__name__'] = 'AimsContext'

    def __repr__(self):
        return "<AimsContext '{}'>" \
            .format(', '.join(k for k in self.keys() if k[0] != '_'))

    def __contains__(self, key):
        return key in self._dict

    def __getattr__(self, key):
        try:
            return self._dict[key]
        except KeyError as e:
            raise AttributeError(
                "'AimsContext' object has no attribute '{}'".format(key)
            ) from e

    def __setattr__(self, key, val):
        raise TypeError('Cannot assign to attributes, modify values in-place')

    def keys(self):
        """View of all names of items in the context."""
        return self._dict.keys()

    def items(self):
        """View of all items in the context."""
        return self._dict.items()

    @property
    def rank(self):
        """MPI rank of the process."""
        return rank

    def evaluate_chi_0(self, r_grid, r_prime, u):
        u = np.array([u.real, u.imag])
        chi_0 = np.empty((r_grid.shape[0],))
        aims.c_evaluate_chi_0(
            ffi.cast('double *', r_grid.ctypes.data),
            r_grid.shape[0],
            ffi.cast('double *', r_prime.ctypes.data),
            ffi.cast('double *', u.ctypes.data),
            ffi.cast('double *', chi_0.ctypes.data),
        )
        return chi_0

    def gather_all_grids(self, arraynames=None):
        """
        Gather grid points and associated quantities from all MPI processes.

        Requires mpi4py. `arraynames` is a list of quantities that are to be
        gathered together with the grid points. Only quantities available in the
        context object can be specified. On the root process, returns a list
        of arrays where the first one lists coordinates of the grid points
        (3, n_pts) and the subsequent correspond to `arraynames` (dim, n_pts).
        On other processes, returns a list of None's of an equal length.

        Examples:
        >>> points, rho, partition_tab = ctx.gather_all_grids(['rho', 'partition_tab'])
        """
        if not comm:
            raise RuntimeError("'gather_all_grids' not available without mpi4py")
        # combine points from all batches
        arrays = [np.array([pt.coords for batch in self.batches for pt in batch.points]).T]
        arraynames = arraynames or []
        n_spin = self.rho.shape[0]
        for name in arraynames:
            if name == 'rho':
                arrays.append(self.rho)
            elif name == 'partition_tab':
                arrays.append(self.partition_tab)
            elif name == 'rho_gradient':
                arrays.append(self.rho_gradient.reshape(3*n_spin, -1))
            elif name == 'kinetic_density':
                arrays.append(self.kinetic_density)
        # stack all arrays to a single one that will be passed by MPI
        total_array = np.vstack(arrays)
        # receive buffer to hold numbers of points on all processes
        all_n_points = np.empty(comm.Get_size(), dtype=int) if rank == 0 else None
        comm.Gather(np.array(total_array.shape[1]), all_n_points)  # gather numbers of points
        n_rows = total_array.shape[0]
        # receive buffer to hold data from all processes
        all_total_array = np.empty((n_rows, all_n_points.sum()), order='F') \
            if rank == 0 else None
        # Since the numbers of points differ accross processes, we have to use
        # MPI_Gatherv instead of MPI_Gather. The receive buffer for Gatherv is
        # a tuple (buffer, counts, displacements, type), where counts are the
        # sizes of the data across processes and displacements are their
        # intended locations in the receive buffer.
        comm.Gatherv(total_array, (
            all_total_array,
            tuple(n_rows*all_n_points),
            (0, *(n_rows*all_n_points).cumsum()[:-1]),
            MPI.DOUBLE
        ) if rank == 0 else None)  # receive buffer only needed on root
        if rank == 0:
            # unpack the single array on root
            ret_arrays = [all_total_array[0:3, :]]
            i_row = 3
            for name in arraynames:
                if name == 'rho':
                    ret_arrays.append(all_total_array[i_row:i_row+n_spin, :])
                    i_row += n_spin
                elif name == 'partition_tab':
                    ret_arrays.append(all_total_array[i_row, :])
                    i_row += 1
                elif name == 'rho_gradient':
                    ret_arrays.append(
                        all_total_array[i_row:i_row+n_spin*3, :]
                        .reshape(3, n_spin, -1)
                    )
                    i_row += n_spin*3
                if name == 'kinetic_density':
                    ret_arrays.append(all_total_array[i_row:i_row+n_spin, :])
                    i_row += n_spin
            assert i_row == n_rows
            return ret_arrays
        else:  # return None's on other processes
            return (1+len(arraynames))*[None]


@ffi.def_extern()
def call_python_cffi(c_ctx, filename, event, rank):
    """
    This is the entry point from Fortran. `c_ctx` is the AimsContext_t C
    struct, `filename` is the Python script or 'REPL', `event` is either
    'parse' or 'run' and `rank` is the MPI id of the caller.
    """
    try:
        filename = ffi.string(filename).decode()
        event = ffi.string(event).decode()
        if event == 'parse':
            retcode = do_parse(filename, rank)
        elif event == 'run':
            retcode = do_run(c_ctx, filename, rank)
    except:
        traceback.print_exc(file=sys.stdout)
        retcode = err_codes['unknown']
    sys.stdout.flush()  # move Python buffer to system buffer
    try:
        # Make sure system buffers are flushed. This is not supported on at
        # least some MPI platforms.
        os.fsync(sys.stdout.fileno())
    except OSError:
        pass
    return retcode


hooks = {}

err_codes = {
    'forming_ctx': 1, 'REPL': 11, 'REPL_exit': 12, 'file_not_exists': 2,
    'load_user_script': 3, 'no_run': 4, 'user_script': 10,
    'user_script_parse': 5, 'unknown': 20
}


def do_parse(filename, rank):
    """
    Load the file as a Python module and execute it's `parse` attribute without
    any arguments when present. Returns -1 on normal exit, other codes on error.
    """
    def rootprint(*args, **kwargs):
        if rank == 0:
            print(*args, **kwargs)
    if filename == 'REPL':
        return -1
    if not Path(filename).is_file():
        rootprint(' *** Error: Python file {} does not exist'.format(filename))
        return err_codes['file_not_exists']
    if rank == 0:
        print('  Trying to import Python file "{}":'.format(filename))
        print('===================================')
        with open(filename) as f:
            print(f.read().strip())
        print('===================================')
    try:
        user_script = load_module(filename)
    except:
        if rank == 0:
            traceback.print_exc(file=sys.stdout)
        rootprint(' *** Error: Cannot load user script')
        return err_codes['load_user_script']
    if not hasattr(user_script, 'run'):
        rootprint(' *** Error: User scripts need to provide function "run(ctx)"')
        return err_codes['no_run']
    if hasattr(user_script, 'parse'):
        try:
            user_script.parse()
        except:
            if rank == 0:
                traceback.print_exc(file=sys.stdout)
            rootprint(' *** Error: User script failed during parsing {}'.format(rank))
            return err_codes['user_script_parse']
    hooks[filename] = user_script
    return -1


def load_module(pathname):
    """Load an arbitrary file as a Python module."""
    path = Path(pathname)
    modulename = path.stem
    module = imp.new_module(modulename)
    exec(compile(path.open().read(), path.name, 'exec'), module.__dict__)
    return module


banner = """\
There is a local variable `ctx`. See `help(ctx)` for details. Press CTRL+D to
continue the aims run. Type `exit(1)` to abort aims."""


def do_run(c_ctx, filename, rank):
    """
    Transorm the C context into a Python context. If `filename` is 'REPL', an
    interactive Python console is launched with the context as a local
    variable. If `filename` is a file, it is loaded as a module and it's `run`
    attribute is executed with the context as a single argument. Returns -1 on
    normal exit, other codes on error.
    """
    def rootprint(*args, **kwargs):
        if rank == 0:
            print(*args, **kwargs)
    try:
        ctx = {}
        ctx['_c_ctx'] = c_ctx
        ctx['_ffi'] = ffi
        ctx['_aims'] = aims
        ctx['coords'] = get_ndarray(c_ctx.coords, (3, c_ctx.n_atoms))
        ctx['elements'] = [
            ffi.string(c_ctx.elements[
                (c_ctx.species[i]-1)*2:(c_ctx.species[i]-1)*2+2
            ]).strip().decode() for i in range(c_ctx.n_atoms)
        ]
        ctx['rho'] = get_ndarray(c_ctx.rho, (c_ctx.n_spin, c_ctx.n_full_points))
        ctx['rho_gradient'] = \
            get_ndarray(c_ctx.rho_gradient, (3, c_ctx.n_spin, c_ctx.n_full_points))
        ctx['kinetic_density'] = \
            get_ndarray(c_ctx.rho_gradient, (c_ctx.n_spin, c_ctx.n_full_points))
        ctx['partition_tab'] = \
            get_ndarray(c_ctx.partition_tab, (c_ctx.n_full_points,))
        ctx['hirshfeld_volume'] = \
            get_ndarray(c_ctx.hirshfeld_volume, (c_ctx.n_atoms,))
        ctx['freeintegral'] = get_ndarray(c_ctx.freeintegral, (c_ctx.n_atoms,))
        ctx['hirshfeldw'] = get_ndarray(c_ctx.hirshfeldw, (c_ctx.n_atoms,))
        ctx['batches'] = [
            BatchOfPoints(
                size=batch.size,
                points=[GridPoint(
                    coords=get_ndarray(point.coords, (3,)),
                    index_atom=point.index_atom,
                    index_radial=point.index_radial,
                    index_angular=point.index_radial
                ) for point in citer(batch.points, batch.size)]
            ) for batch in citer(c_ctx.batches, c_ctx.n_my_batches)
        ]
        try:
            ctx['KS_eigenvector'] = get_ndarray(
                c_ctx.KS_eigenvector, (
                    c_ctx.n_basis, c_ctx.n_states, c_ctx.n_spin,
                    c_ctx.n_k_points_task
                )
            ) if c_ctx.real_eigenvectors else get_ndarray(
                c_ctx.KS_eigenvector_complex, (
                    c_ctx.n_basis, c_ctx.n_states, c_ctx.n_spin,
                    c_ctx.n_k_points_task
                ),
                dtype='complex'
            )
        except TypeError:
            ctx['KS_eigenvector'] = None
        ctx['KS_eigenvalue'] = get_ndarray(c_ctx.KS_eigenvalue, (
            c_ctx.n_states, c_ctx.n_spin, c_ctx.n_k_points
        ))
        ctx['occ_numbers'] = get_ndarray(c_ctx.occ_numbers, (
            c_ctx.n_states, c_ctx.n_spin, c_ctx.n_k_points
        ))
        ctx['hamiltonian'] = get_ndarray(c_ctx.hamiltonian, (
            c_ctx.n_hamiltonian_matrix_size, c_ctx.n_spin
        ))
        ctx = AimsContext(ctx)
    except:
        traceback.print_exc(file=sys.stdout)
        return err_codes['forming_ctx']
    if filename == 'REPL':
        try:
            code.interact(local={'ctx': ctx}, banner=banner)
        except SystemExit as e:
            if e.args[0] != 0:
                return err_codes['REPL_exit']
        except:  # not sure whether this can happen
            traceback.print_exc(file=sys.stdout)
            return err_codes['REPL']
        return -1
    try:
        hooks[filename].run(ctx)
    except:
        if rank == 0:
            traceback.print_exc(file=sys.stdout)
        rootprint(' *** Error: User script failed on process {}'.format(rank))
        return err_codes['user_script']
    return -1


def citer(obj, siz):
    """Wrap a C pointer as an interator."""
    for i in range(siz):
        yield obj[i]


def get_ndarray(pointer, shape, dtype='float'):
    """Wrap a C pointer as a Numpy array."""
    return np.ndarray(
        shape=shape,
        buffer=ffi.buffer(pointer, np.prod(shape)*8),
        order='F',
        dtype=dtype
    )
