#!/usr/bin/env python

import struct
import numpy as np
import scipy.sparse as sp

def read_elsi_to_csr(filename):
    mat = open(filename,"rb")
    data = mat.read()
    mat.close()

    # Get header
    start = 0
    end = 64
    header = struct.unpack("i"*16,data[start:end])

    # Number of basis functions (matrix size)
    n_basis = header[3]

    # Total number of non-zero elements
    nnz = header[5]

    # Get column pointer
    start = end
    end = start+n_basis*4
    col_ptr = struct.unpack("i"*n_basis,data[start:end])
    col_ptr += (nnz+1,)
    col_ptr = np.array(col_ptr)

    # Get row index
    start = end
    end = start+nnz*4
    row_idx = struct.unpack("i"*nnz,data[start:end])
    row_idx = np.array(row_idx)

    # Get non-zero value
    start = end
    end = start+nnz*8
    nnz_val = struct.unpack("d"*nnz,data[start:end])
    nnz_val = np.array(nnz_val)

    # Change convention
    for i_val in range(nnz):
        row_idx[i_val] -= 1

    for i_col in range(n_basis+1):
        col_ptr[i_col] -= 1

    return sp.csr_matrix((nnz_val,row_idx,col_ptr),shape=(n_basis,n_basis))

def read_elsi_to_csr_complex(filename):
    mat = open(filename,"rb")
    data = mat.read()
    mat.close()

    # Get header
    start = 0
    end = 64
    header = struct.unpack("i"*16,data[start:end])

    # Number of basis functions (matrix size)
    n_basis = header[3]

    # Total number of non-zero elements
    nnz = header[5]

    # Get column pointer
    start = end
    end = start+n_basis*4
    col_ptr = struct.unpack("i"*n_basis,data[start:end])
    col_ptr += (nnz+1,)
    col_ptr = np.array(col_ptr)

    # Get row index
    start = end
    end = start+nnz*4
    row_idx = struct.unpack("i"*nnz,data[start:end])
    row_idx = np.array(row_idx)

    # Get non-zero value
    start = end
    end = start+nnz*16
    nnz_val = struct.unpack("d"*nnz*2,data[start:end])
    nnz_val_real = np.array(nnz_val[0::2])
    nnz_val_imag = np.array(nnz_val[1::2])
    nnz_val = nnz_val_real + 1j*nnz_val_imag

    # Change convention
    for i_val in range(nnz):
        row_idx[i_val] -= 1

    for i_col in range(n_basis+1):
        col_ptr[i_col] -= 1

    return sp.csr_matrix((nnz_val,row_idx,col_ptr),shape=(n_basis,n_basis))
