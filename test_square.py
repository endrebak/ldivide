mid1, mid2 = 4, 10

for i in range(mid1, mid2):
    print(f"--- {i} ---")
    for j in range(i + 1, mid2):
        print((i, j))

