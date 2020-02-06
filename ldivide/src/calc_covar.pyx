#!/usr/bin/env python3
#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

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
cpdef calc_covar(haplos, int32_t [::1] autocovars, double theta, int32_t window_size):

    cdef:
        int i, j, j1, j2, k, pos1, pos2, len_haps, n11, n10, n01, len_g1_int, half_window_size
        double len_g1, f11, f1, f2, Ds2, D, nind, thetas, aj1, aj2

    half_window_size = int(window_size / 2)
    haps = haplos.values

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    len_haps = len(haps)

    assert len_haps >= 2 * half_window_size

    thetas = (1-theta)*(1-theta)

    cdef long[:] allpos_view
    cdef long[:] g1_view, g2_view
    cdef long[:, :] haps_view
    cdef double[:] outvec_view

    allpos_view = allpos.values
    haps_view = haps
    outvec = np.zeros(len_haps)

    for i in range(0, half_window_size):

            for j in range(i):

                j1 = i - j
                j2 = i + j

                if j1 < 0:
                    break

                aj2 = autocovars[j1]
                n11, n01, n10 = 0, 0, 0
                for k in range(len_g1_int):
                    if haps_view[j1][k] == 1 and haps_view[j2][k] == 1:
                        n11 += 1
                    elif haps_view[j1][k] == 0 and haps_view[j2][k] == 1:
                        n01 += 1
                    elif haps_view[j1][k] == 1 and haps_view[j2][k] == 0:
                        n10 += 1

                f11 = n11/len_g1
                f1 = (n11+n10)/len_g1
                f2 = (n11+n01)/len_g1
                D = f11 - f1*f2
                Ds2 = (thetas*D)
                Ds2 = Ds2 * Ds2

                outvec_view[i] += Ds2 / (aj1 * aj2)


    for i in range(half_window_size, len_haps - half_window_size):

            for j in range(half_window_size):

                j1 = i - j
                j2 = i + j

                aj1 = autocovars[j1]
                aj2 = autocovars[j2]

                n11, n01, n10 = 0, 0, 0
                for k in range(len_g1_int):
                    if haps_view[j1][k] == 1 and haps_view[j2][k] == 1:
                        n11 += 1
                    elif haps_view[j1][k] == 0 and haps_view[j2][k] == 1:
                        n01 += 1
                    elif haps_view[j1][k] == 1 and haps_view[j2][k] == 0:
                        n10 += 1

                f11 = n11/len_g1
                f1 = (n11+n10)/len_g1
                f2 = (n11+n01)/len_g1
                D = f11 - f1*f2
                Ds2 = (thetas*D)
                Ds2 = Ds2 * Ds2

                outvec_view[i] += Ds2 / (aj1 * aj2)


    for i in range(len_haps - half_window_size, len_haps):

            for j in range(i):

                j1 = i - j
                j2 = i + j

                if j2 >= len_haps:
                    break

                aj1 = autocovars[j1]
                aj2 = autocovars[j2]

                n11, n01, n10 = 0, 0, 0
                for k in range(len_g1_int):
                    if haps_view[j1][k] == 1 and haps_view[j2][k] == 1:
                        n11 += 1
                    elif haps_view[j1][k] == 0 and haps_view[j2][k] == 1:
                        n01 += 1
                    elif haps_view[j1][k] == 1 and haps_view[j2][k] == 0:
                        n10 += 1

                f11 = n11/len_g1
                f1 = (n11+n10)/len_g1
                f2 = (n11+n01)/len_g1
                D = f11 - f1*f2
                Ds2 = (thetas*D)
                Ds2 = Ds2 * Ds2

                outvec_view[i] += Ds2 / (aj1 * aj2)
