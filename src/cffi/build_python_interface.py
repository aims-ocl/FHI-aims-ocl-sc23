import cffi
from pathlib import Path

ffi = cffi.FFI()

srcdir = Path(__file__).resolve().parent

with (srcdir/'python_interface.h').open() as f:
    ffi.embedding_api(f.read())
with (srcdir/'../aims_c_api.h').open() as f:
    ffi.cdef(f.read())
ffi.set_source(
    'c_aims',
    '#include "{}/python_interface.h"\n'.format(srcdir) +
    '#include "{}/aims_c_api.h"'.format(srcdir/'..')
)
with (srcdir/'python_interface.py').open() as f:
    ffi.embedding_init_code(f.read())
ffi.emit_c_code('python_interface.c')
