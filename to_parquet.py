import pandas as pd
from sys import stdin, argv
import numpy as np

outfile = argv[1]
samples_file = argv[2]

n = 0
with open(samples_file) as f:
    for _ in f:
        n += 1

# times two for variants, +1 for index
n = (n * 2) + 1

dtypes = {k: np.int8 for k in range(n)}
dtypes[0] = np.int64

df = pd.read_table(stdin, header=None, dtype=dtypes, index_col=0, sep="\s+")
# print(df)
valid_haplos = ((df == 1) | (df == 0)).all(axis=1)
# print()
# print(valid_haplos)
df = df[valid_haplos]
# print(df)
df.columns = df.columns.astype(str)
df.index.name = "Position"

df.to_parquet(outfile)
