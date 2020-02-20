import sys, math, gzip
import numpy as np
import pandas as pd

from time import time

from libc.math cimport exp, fabs
from libc.stdint cimport int32_t, int8_t

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cpdef calc_colvec(const double [::1] rows, const double [::1] cols, int32_t window_size):

    cdef:
        int i, half_window_size, len_rows

    half_window_size = int(window_size / 2)

    len_rows = len(rows)

    assert len(rows) == len(cols)

    outvec = np.zeros(len_rows)
    cdef double[::1] outvec_view
    outvec_view = outvec

    outvec_view[0] = rows[i]

    for i in range(1, len_rows):
        outvec_view[i] = outvec_view[i - 1] - cols[i - 1] + rows[i]

    return outvec
