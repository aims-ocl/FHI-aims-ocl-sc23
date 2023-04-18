#!/usr/bin/env python

import sys
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Read into csr
csr = io.mmread(sys.argv[1]).tocsr()

# Find sparsity pattern
spp = np.zeros(csr.shape,dtype=np.int8)

i_val = 0

for i_col in range(csr.shape[1]):
    this_nnz = csr.indptr[i_col+1]-csr.indptr[i_col]

    for j_val in range(i_val,i_val+this_nnz-1):
        if abs(csr.data[j_val]) > 0:
            tmp = 1
        else:
            tmp = 0

        spp[csr.indices[j_val],i_col] = tmp

    i_val += this_nnz

# Plot
plt.matshow(spp,cmap=cm.gray_r)
plt.show()
