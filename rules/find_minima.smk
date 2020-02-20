import math
import bisect
from time import time

from ldivide.src.calc_square import calc_rowvec, calc_colvec, calc_final

import scipy.ndimage.filters as filters
import scipy.signal as sig

import numpy as np

def apply_filter_get_minima(vector, width):

    start = time()
    a = sig.get_window('hanning', 2 * (width + 1))

    ga = filters.convolve1d(vector, a/a.sum())

    minima_a = sig.argrelextrema(ga, np.less)[0]
    end = time()
    print(f"(Took {end - start} seconds.)")

    return minima_a


rule find_minima:
    input:
        "{prefix}/1kg/vector/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/minima/{population}/{chromosome}.pq"
    run:
        f = input[0]
        vector = pd.read_parquet(f).squeeze()
        vector.name = "Value"

        n_snps_between_breakpoints = 5000

        n_bpoints = int(math.ceil( len(vector) / n_snps_between_breakpoints - 1 ))

        """First do linear search to find max filter width, then do binary search to find filter width that
gets you close to desired number of breakpoints."""

        # nb_minima = int(1e6)
        search_val = 100
        last_search_val = 50
        while True:
            print(f"Trying a search width of {search_val}", end="\n")
            minima = apply_filter_get_minima(vector, search_val)
            if len(minima) <= n_bpoints:
                print(f"Found {len(minima)} minima, done with initial search.")
                break

            print(f"Found {len(minima)} minima, searching for something less than {n_bpoints} in initial search.")
            last_search_val = search_val
            search_val *= 2

        if not len(minima) == n_bpoints:
            last_nb_minima = -1000
            lo = last_search_val
            hi = search_val

            while lo <= hi:
                mid = (lo + hi) // 2

                print("---" * 10)
                print(f"Trying a search width of {mid}", end="\n")
                minima = apply_filter_get_minima(vector, mid)
                nb_minima = len(minima)

                print(f"Found {nb_minima} minima, searching for something in ({n_bpoints - 10}, {n_bpoints + 10}).")

                if abs(nb_minima - n_bpoints) <= 10:
                    break

                if nb_minima > n_bpoints:
                    lo = mid + 1
                elif n_bpoints > nb_minima:
                    hi = mid - 1
                else:
                    break

        pd.DataFrame({"Minima": minima}).to_parquet(output[0])


rule compute_row:
    input:
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt",
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/rowvectors/{population}/{chromosome}.pq"
    run:
        theta = float(open(input.theta).readline().strip())
        haplos = pd.read_parquet(input.haplos)
        autocovars = pd.read_parquet(input.autocorr).squeeze().values

        result = calc_rowvec(haplos.values, autocovars, theta, 5000)
        result = pd.Series(result)
        result = result.to_frame()
        result.columns = result.columns.astype(str)
        result.to_parquet(output[0])


rule compute_col:
    input:
        haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
        theta = "{prefix}/theta/{chromosome}/{population}.txt",
        autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/colvectors/{population}/{chromosome}.pq"
    run:
        theta = float(open(input.theta).readline().strip())
        haplos = pd.read_parquet(input.haplos)
        autocovars = pd.read_parquet(input.autocorr).squeeze().values

        result = calc_colvec(haplos.values, autocovars, theta, 5000)
        result = pd.Series(result)
        print(result)
        result = result.to_frame()
        result.columns = result.columns.astype(str)
        print(result)
        result.to_parquet(output[0])


rule compute_square:
    input:
        col = "{prefix}/1kg/colvectors/{population}/{chromosome}.pq",
        row = "{prefix}/1kg/rowvectors/{population}/{chromosome}.pq",
    output:
        "{prefix}/1kg/square/{population}/{chromosome}.pq"
    run:
        col = pd.read_parquet(input.col).squeeze()
        row = pd.read_parquet(input.row).squeeze()

        print(row.isnull().sum())
        # print(col.isnull().sum())
        final = calc_final(row.values, col.values, 5000)

        df = pd.Series(final).to_frame()
        print(df)
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])


rule normalize_square:
    input:
        "{prefix}/1kg/square/{population}/{chromosome}.pq"
    output:
        "{prefix}/1kg/normalized_square/{population}/{chromosome}.pq"
    run:
        s = pd.read_parquet(input[0]).squeeze()

        l = len(s)
        n = np.ones(l) * ((2499 * 2498)/2)
        flank_n = (np.arange(2, 2501) * (np.arange(2, 2501) - 1))/2
        n[:2499] = flank_n
        n[len(n)-2499:] = flank_n[::-1]

        print(s)
        s = s / n
        print(s)
        df = s.to_frame()
        df.columns = df.columns.astype(str)
        df.to_parquet(output[0])





# rule find_local_minima:
#     input:
#         minima = "{prefix}/1kg/minima/{population}/{chromosome}.pq",
#         haplos = "{prefix}/1kg/{population}/{chromosome}.pq",
#         theta = "{prefix}/theta/{chromosome}/{population}.txt",
#         autocorr = "{prefix}/1kg/autocorr/{population}/{chromosome}.pq",
#     output:
#         "{prefix}/1kg/local_minima/{population}/{chromosome}.pq"
#     run:
#         theta = float(open(input.theta).readline().strip())
#         haplos = pd.read_parquet(input.haplos)
#         minima = pd.read_parquet(input.minima)
#         autocovars = pd.read_parquet(input.autocorr).squeeze().values

#         a = minima.squeeze().values
#         diff_last = len(autocovars) - a[-1]
#         shifted = a + np.r_[(np.diff(a)/2).astype(int), diff_last]
#         midpoints = np.r_[0, shifted]

#         new_breakpoints = []

#         for i, (m1, m2) in enumerate(zip(midpoints, midpoints[1:])):
#             print("-----" * 5)
#             print(np.round(i/len(midpoints), 3))
#             haps = haplos.iloc[m1:m2].values
#             auto = autocovars[m1:m2]

#             new_breakpoint = find_local_minima(haps, auto, theta) + m1
#             new_breakpoints.append(new_breakpoint)
#             print("Distance from border:", min(abs(new_breakpoint - m2), abs(new_breakpoint - m1)))
#             print("Distance from old breakpoint:", abs(new_breakpoint - a[i]))

#         new_breakpoints = np.array(new_breakpoints)

#         pd.DataFrame({"Breakpoints": new_breakpoints}).to_parquet(output[0])
