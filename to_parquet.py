import pandas as pd
from sys import stdin, argv
import numpy as np

outfile = argv[1]

df = pd.read_table(stdin, header=None, dtype=np.int8, sep="\s+", index_col=0)
invalid_haplos = (df != 1) | (df != 0)
df = df[~invalid_haplos.any(axis=1)]
df.columns = df.columns.astype(str)

df.to_parquet(outfile)
