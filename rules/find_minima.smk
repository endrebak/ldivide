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
        diag = "{prefix}/1kg/vector/{population}/{window_size}/{chromosome}.pq",
        other = "{prefix}/1kg/normalized_square/{population}/{window_size}/{chromosome}.pq"
    output:
        statistics = "{prefix}/1kg/minima/{population}/{window_size}/{chromosome}.tsv",
        minima_diag = "{prefix}/1kg/minima/{population}/{window_size}/diag/{chromosome}.tsv",
        minima_other = "{prefix}/1kg/minima/{population}/{window_size}/other/{chromosome}.tsv"
    run:
        f = input.diag
        f2 = input.other

        v = pd.read_parquet(f).squeeze()
        v2 = pd.read_parquet(f2).squeeze()

        rowdicts = []
        minima_diag = []
        minima_other = []

        for width in range(0, 30001, 1000):
            r = apply_filter_get_minima(v, width)
            r2 = apply_filter_get_minima(v2, width)

            df = pd.Series(r).to_frame()
            df.insert(df.shape[1], "Type", "diagonal")
            df.insert(df.shape[1], "Width", width)

            df2 = pd.Series(r2).to_frame()
            df2.insert(df2.shape[1], "Type", "other")
            df2.insert(df2.shape[1], "Width", width)
            print(f"----{width}----")
            print(len(r))
            print(len(r2))
            lr = len(r)
            lr2 = len(r2)
            rowdicts.append({"Chromosome": wildcards.chromosome, "Width": width, "Minima": lr, "Type": "diagonal"})
            rowdicts.append({"Chromosome": wildcards.chromosome, "Width": width, "Minima": lr2, "Type": "other"})

            minima_diag.append(df)
            minima_other.append(df2)

        df = pd.DataFrame(rowdicts)

        df.to_csv(output.statistics, sep="\t", index=False)

        diag = pd.concat(minima_diag)
        other = pd.concat(minima_other)

        diag.to_csv(output.minima_diag, sep="\t", index=False)
        other.to_csv(output.minima_other, sep="\t", index=False)


rule find_local_minima:
    input:
        minima = "{prefix}/1kg/minima/{population}/{window_size}/{chromosome}.pq",
        vector = "{prefix}/1kg/normalized_square/{population}/{window_size}/{chromosome}.pq"
    output:
        "{prefix}/1kg/local_minima/{population}/{window_size}/{chromosome}.tsv"
    run:
        vector = pd.read_parquet(input.vector).squeeze()
        v = vector.values
        minima = pd.read_parquet(input.minima)

        a = minima.squeeze().values
        diff_last = len(vector) - a[-1]
        shifted = a + np.r_[(np.diff(a)/2).astype(int), diff_last]
        midpoints = np.r_[1, shifted] # not zero-idx because first snp is always zero

        new_breakpoints = []

        for i, (m1, m2) in enumerate(zip(midpoints, midpoints[1:])):

            new_vector = vector[m1:m2]
            new_breakpoint_value = new_vector.min()
            new_breakpoint = new_vector.idxmin()
            # equal_value = new_vector[new_vector == new_breakpoint_value]
            # print(equal_value)

            midpoint_distance = min(abs(new_breakpoint - m2), abs(new_breakpoint - m1))
            old_breakpoint_distance = abs(new_breakpoint - a[i])

            old_value = vector.values[a[i]]
            new_value = vector.values[new_breakpoint]

            new_breakpoints.append({"Chromosome": wildcards.chromosome,
                                    "MidpointStart": m1,
                                    "MidpointEnd": m2,
                                    "NewBreakpoint": new_breakpoint,
                                    "OldBreakpoint": a[i],
                                    "BorderDistance": min(abs(new_breakpoint - m2), abs(new_breakpoint - m1)),
                                    "OldBreakpointDistance": abs(new_breakpoint - a[i]),
                                    "OldValue": old_value,
                                    "NewValue": new_value,
                                    "ValueDifference": abs(old_value - new_value),
            })

        new_breakpoints = pd.DataFrame.from_dict(new_breakpoints)

        new_breakpoints.to_csv(output[0], sep="\t", index=False)


# rule collect_minima:
#     input:
#         expand("{{prefix}}/1kg/local_minima/{{population}}/{chromosome}.tsv", chromosome=chromosomes)
#     output:
#         "{prefix}/1kg/local_minima/{population}/genome.tsv"
#     run:
#         df = pd.concat([pd.read_table(f) for f in input])
#         print("df.median()")
#         print(df.median())
#         print(("df[df.BorderDistance > df.OldBreakpointDistance].median()"))
#         print((df.BorderDistance > df.OldBreakpointDistance).sum())
#         print(df[df.BorderDistance > df.OldBreakpointDistance].median())
#         print(("df[df.BorderDistance < df.OldBreakpointDistance].median()"))
#         print((df.BorderDistance < df.OldBreakpointDistance).sum())
#         print(df[df.BorderDistance < df.OldBreakpointDistance].median())
#         print("df.mean()")
#         print(df.mean())
#         print(("df[df.BorderDistance > df.OldBreakpointDistance].mean()"))
#         print((df.BorderDistance > df.OldBreakpointDistance).sum())
#         print(df[df.BorderDistance > df.OldBreakpointDistance].mean())
#         print(("df[df.BorderDistance < df.OldBreakpointDistance].mean()"))
#         print((df.BorderDistance < df.OldBreakpointDistance).sum())
#         print(df[df.BorderDistance < df.OldBreakpointDistance].mean())
#         # print("df.mean()")
#         # print(df.mean())
#         # print("df.min()")
#         # print(df.min())
#         # print("df.max()")
#         # print(df.max())
#         df.to_csv(output[0], sep="\t", index=False)

