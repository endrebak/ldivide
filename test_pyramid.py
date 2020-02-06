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

# 1 1 1
# 2 2 2
# 2 1 3
# 3 3 3
# 3 2 4
# 3 1 5
# 4 4 4
# 4 3 5
# 4 2 6
# 5 5 5
# 5 4 6
# 5 3 7
# 6 6 6
# 6 5 7
# 6 4 8
# 7 7 7
# 7 6 8
# 7 5 9
# 8 8 8
# 8 7 9
# 9 9 9
