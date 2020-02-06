"""Just a sanity-check for the code in calc_covar.pyx"""


half_window_size = 3
len_haps = 10

# half_window_size += 1

for i in range(0, half_window_size):

        for j in range(i):

            j1 = i - j
            j2 = i + j

            if j1 < 0:
                break

            print(i, j1, j2)

for i in range(half_window_size, len_haps - half_window_size):

        for j in range(half_window_size):

            j1 = i - j
            j2 = i + j

            print(i, j1, j2)

for i in range(len_haps - half_window_size, len_haps):

        for j in range(i):

            j1 = i - j
            j2 = i + j

            if j2 >= len_haps:
                break

            print(i, j1, j2)
