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
cpdef calc_covar(haplos, double [::1] autocovars, double theta):

    cdef:
        int mid_x, mid1, mid2, i, j, k, len_haps, width_haps
        double rowsum
        int8_t[:, :] haps_view
        double[::1] colvec_view
        double[::1] rowvec_view
        int32_t[::1] new_breakpoints

    haps_view = haplos
    len_haps = len(haplos)
    width_haps = haplos.shape[1]
    thetas = (1-theta)*(1-theta)
    colvec = np.zeros(len_haps)
    colvec_view = colvec

    rowvec = np.zeros(len_haps)
    rowvec_view = rowvec

    for i in range(0, len_haps):
        for j in range(i + 1, len_haps):

            ai = autocovars[i]
            aj = autocovars[j]

            n11, n01, n10 = 0, 0, 0
            for k in range(width_haps):
                if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                    n11 += 1
                elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                    n01 += 1
                elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                    n10 += 1

            f11 = n11/len_haps
            f1 = (n11+n10)/len_haps
            f2 = (n11+n01)/len_haps
            D = f11 - f1*f2
            Ds2 = (thetas*D)
            Ds2 = Ds2 * Ds2
            rsq = Ds2 / (ai * aj)

            rowvec_view[i] += rsq
            colvec_view[j] += rsq

        if i > 1:
            rowvec_view[i] = rowvec_view[i] + rowvec_view[i-1] - colvec_view[i-1]

    for i in range(0, len_haps):
        if i > 0:
            rowvec_view[i] = rowvec_view[i] / (i * (len_haps - i))
        else:
            rowvec_view[i] = rowvec_view[i] / len_haps

    breakpoint = np.argmin(rowvec_view)

    return breakpoint

