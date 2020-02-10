import math
import bisect
from time import time

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
        "{prefix}/1kg/vector/{population}/{chromosome}.tsv.gz"
    output:
        "{prefix}/1kg/minima/{population}/{chromosome}.tsv.gz"
    run:
        f = input[0]
        vector = pd.read_table(f, sep="\t", header=None, names=["Value"], squeeze=True)

        n_snps_between_breakpoints = 5000

        n_bpoints = int(math.ceil( len(vector) / n_snps_between_breakpoints - 1 ))

        """Increase width until less than bp, then do linear searches around that point."""

        # nb_minima = int(1e6)
        search_val = 100
        last_search_val = 50
        while True:
            print(f"Trying a search width of {search_val}", end="\n")
            minima = apply_filter_get_minima(vector, search_val)
            print(f"Found {len(minima)} minima, searching for something greater than {n_bpoints} in initial search.")
            if len(minima) <= n_bpoints:
                break
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





