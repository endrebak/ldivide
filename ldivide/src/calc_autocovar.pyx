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
cpdef calc_autocovar(int8_t [::1] haplos, double theta):

    cdef:
        int i, j, k, pos1, pos2, len_allpos, n11, n10, n01, len_g1_int
        double len_g1, f11, f1, f2, Ds2, D, nind, thetas, thetas2

    haps = haplos.values

    records = []

    len_g1_int = haps.shape[1]
    len_g1 = float(haps.shape[1])
    len_allpos = len(allpos)
    thetas = (1-theta)*(1-theta)
    thetas2 = (theta/2.0)*(1-theta/2.0)

    cdef long[:] allpos_view
    cdef long[:] g1_view, g2_view
    cdef long[:, :] haps_view

    allpos_view = allpos.values
    haps_view = haps

    for i in range(len_allpos):

            n11, n01, n10 = 0, 0, 0
            for k in range(len_g1_int):
                if haps_view[i][k] == 1 and haps_view[i][k] == 1:
                    n11 += 1

                f11 = n11/len_g1
                D = f11
                Ds2 = thetas*D

                Ds2 = Ds2 + thetas2

                print(Ds2)
                j = j+1
