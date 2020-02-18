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
cpdef calc_colvec(haplos, double [::1] autocovars, double theta, int32_t window_size):

    cdef:
        int i, j, j1, j, k, pos1, pos2, len_haps, n11, n10, n01, len_g1_int, half_window_size
        double len_g1, f11, f1, f2, Ds2, D, nind, thetas, aj1, aj, rsq

    half_window_size = int(window_size / 2)
    haps = haplos

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    len_haps = len(haps)

    assert len_haps >= 2 * half_window_size

    thetas = (1-theta)*(1-theta)

    cdef int8_t[:, :] haps_view
    cdef double[::1] outvec_view

    haps_view = haps
    outvec = np.zeros(len_haps)
    outvec_view = outvec

    for i in range(0, half_window_size):

        for j in range(i):

            j = i - j

            if j < 0:
                continue

            ai = autocovars[i]
            aj = autocovars[j]
            n11, n01, n10 = 0, 0, 0
            for k in range(len_g1_int):
                if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                    n11 += 1
                elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                    n01 += 1
                elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                    n10 += 1

            f11 = n11/len_g1
            f1 = (n11+n10)/len_g1
            f2 = (n11+n01)/len_g1
            D = f11 - f1*f2
            Ds2 = (thetas*D)
            Ds2 = Ds2 * Ds2
            rsq = Ds2 / (ai * aj)

            outvec_view[i] += rsq


    for i in range(half_window_size, len_haps):

        for j in range(half_window_size):

            j = i - j

            ai = autocovars[i]
            aj = autocovars[j]

            n11, n01, n10 = 0, 0, 0
            for k in range(len_g1_int):
                if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                    n11 += 1
                elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                    n01 += 1
                elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                    n10 += 1

            f11 = n11/len_g1
            f1 = (n11+n10)/len_g1
            f2 = (n11+n01)/len_g1
            D = f11 - f1*f2
            Ds2 = (thetas*D)
            Ds2 = Ds2 * Ds2
            rsq = Ds2 / (ai * aj)

            outvec_view[i] += rsq


    return outvec



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cpdef calc_rowvec(haplos, double [::1] autocovars, double theta, int32_t window_size):

    cdef:
        int i, j, j1, j, k, pos1, pos2, len_haps, n11, n10, n01, len_g1_int, half_window_size
        double len_g1, f11, f1, f2, Ds2, D, nind, thetas, aj1, aj, rsq

    half_window_size = int(window_size / 2)
    haps = haplos

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    len_haps = len(haps)

    assert len_haps >= 2 * half_window_size

    thetas = (1-theta)*(1-theta)

    cdef int8_t[:, :] haps_view
    cdef double[::1] outvec_view

    haps_view = haps
    outvec = np.zeros(len_haps)
    outvec_view = outvec

    for i in range(0, len_haps - half_window_size):

        for j in range(half_window_size):

            j = i + j

            ai = autocovars[i]
            aj = autocovars[j]
            n11, n01, n10 = 0, 0, 0
            for k in range(len_g1_int):
                if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                    n11 += 1
                elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                    n01 += 1
                elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                    n10 += 1

            f11 = n11/len_g1
            f1 = (n11+n10)/len_g1
            f2 = (n11+n01)/len_g1
            D = f11 - f1*f2
            Ds2 = (thetas*D)
            Ds2 = Ds2 * Ds2
            rsq = Ds2 / (ai * aj)

            outvec_view[i] += rsq


    for i in range(len_haps - half_window_size, len_haps):

        for j in range(half_window_size):

            j = i + j

            if j > len_haps:
                continue

            ai = autocovars[i]
            aj = autocovars[j]

            n11, n01, n10 = 0, 0, 0
            for k in range(len_g1_int):
                if haps_view[i][k] == 1 and haps_view[j][k] == 1:
                    n11 += 1
                elif haps_view[i][k] == 0 and haps_view[j][k] == 1:
                    n01 += 1
                elif haps_view[i][k] == 1 and haps_view[j][k] == 0:
                    n10 += 1

            f11 = n11/len_g1
            f1 = (n11+n10)/len_g1
            f2 = (n11+n01)/len_g1
            D = f11 - f1*f2
            Ds2 = (thetas*D)
            Ds2 = Ds2 * Ds2
            rsq = Ds2 / (ai * aj)

            outvec_view[i] += rsq

    return outvec
