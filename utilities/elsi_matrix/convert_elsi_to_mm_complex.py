#!/usr/bin/env python

import sys
import scipy.io as io
from read_elsi import read_elsi_to_csr_complex

io.mmwrite(sys.argv[2],read_elsi_to_csr_complex(sys.argv[1]))
