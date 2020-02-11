import pandas as pd
from sys import stdin, argv

outfile = argv[1]

df = pd.read_table(stdin, header=None, dtype=bool, sep="\s+")
df.columns = df.columns.astype(str)

df.to_parquet(outfile)
